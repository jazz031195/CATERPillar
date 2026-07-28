#include "openglwindow.h"
#include <QSurfaceFormat>
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QtMath>
#include <QWheelEvent>
#include <QMouseEvent>
#include <random>
#include <future>
#include <thread>
#include <algorithm>
#include <mutex>
#include <vector>
#include <QPainter>


OpenGLWindow::OpenGLWindow(QWindow *parent)
    : QOpenGLWindow(NoPartialUpdate, parent),
      orbitTheta(0.0f), 
      orbitPhi(0.0f), 
      cameraDistance(100.0f), // Try bumping this to 100.0f if the axons are huge
      zoomFactor(1.0f),      // <--- CRITICAL: Must be 1.0, not 0
      totalBatchedVertices(0)
{
    // By default QSurfaceFormat has depthBufferSize == 0, which means the
    // window is created WITHOUT a depth buffer. glEnable(GL_DEPTH_TEST) is
    // then silently a no-op, and spheres get drawn in draw-call order
    // instead of true depth order (back objects overpaint front ones).
    // We must request a depth buffer before the underlying native window
    // is created.
    QSurfaceFormat format;
    format.setDepthBufferSize(24);
    format.setStencilBufferSize(8);
    format.setSamples(4); // optional: anti-aliasing
    setFormat(format);

    connect(&timer, &QTimer::timeout, this, static_cast<void(QWindow::*)()>(&OpenGLWindow::update));
    timer.start(16); // Refresh every ~16 ms (60 FPS)
}
OpenGLWindow::~OpenGLWindow() {}

template <typename T>
constexpr const T& clamp(const T& v, const T& lo, const T& hi)
{
    return (v < lo) ? lo : (hi < v) ? hi : v;
}

namespace {
    inline QColor lerp(const QColor& a, const QColor& b, double t) {
        t = clamp(t, 0.0, 1.0);
        return QColor(
            static_cast<int>(a.red()   + t * (b.red()   - a.red())),
            static_cast<int>(a.green() + t * (b.green() - a.green())),
            static_cast<int>(a.blue()  + t * (b.blue()  - a.blue()))
        );
    }
}

void OpenGLWindow::buildBatchedGeometry()
{
    geometryReady = false;
    batchedVertices.clear();
    batchedColors.clear();
    batchedRadii.clear(); // Clear the new array

    size_t numSpheres = spherePositions.size();
    if (numSpheres == 0) return;
    
    batchedVertices.reserve(numSpheres * 3);
    batchedColors.reserve(numSpheres * 3);
    batchedRadii.reserve(numSpheres); // Reserve space

    for (size_t i = 0; i < numSpheres; ++i) {

        if (isSphereInFrustum(spherePositions[i], sphereRadii[i])) {
            const QVector3D& pos = spherePositions[i];
            const QColor& col = axonColors[i];

            batchedVertices.push_back(pos.x());
            batchedVertices.push_back(pos.y());
            batchedVertices.push_back(pos.z());

            batchedColors.push_back(col.redF());
            batchedColors.push_back(col.greenF());
            batchedColors.push_back(col.blueF());
            
            batchedRadii.push_back(sphereRadii[i]); // Save the specific radius
        }
    }
    
    totalBatchedVertices = batchedVertices.size() / 3;
    geometryReady = true;
}

void OpenGLWindow::setSpheres(const std::vector<std::vector<double>>& x,
                              const std::vector<std::vector<double>>& y,
                              const std::vector<std::vector<double>>& z,
                              const std::vector<std::vector<double>>& radius,
                              const std::vector<int>& groupIds)
{
    // --- 1. Reset and Pre-calculate Colors (Fast) ---
    spherePositions.clear();
    sphereRadii.clear();
    axonColors.clear();

    if (x.empty()) return;

    std::vector<QColor> colors(x.size());
    std::vector<size_t> idxAxon, idxG1, idxG2, idxG3, idxBV;
    for (size_t i = 0; i < groupIds.size(); ++i) {
        switch (static_cast<SphereGroup>(groupIds[i])) {
            case SphereGroup::Axon:   idxAxon.push_back(i); break;
            case SphereGroup::Glial1: idxG1.push_back(i);   break;
            case SphereGroup::Glial2: idxG2.push_back(i);   break;
            case SphereGroup::Glial3: idxG3.push_back(i);   break;
            case SphereGroup::Blood:  idxBV.push_back(i);   break;
        }
    }

    auto assignGradient = [&](const std::vector<size_t>& idx, const QColor& c0, const QColor& c1) {
        if (idx.empty()) return;
        for (size_t k = 0; k < idx.size(); ++k) {
            double t = (idx.size() > 1) ? static_cast<double>(k) / (idx.size() - 1) : 1.0;
            colors[idx[k]] = lerp(c0, c1, t);
        }
    };

    assignGradient(idxBV, QColor(255, 200, 200), QColor(170, 0, 0));
    assignGradient(idxG1, QColor(200, 255, 200), QColor(0, 150, 0));
    assignGradient(idxG2, QColor(220, 255, 220), QColor(0, 110, 0));
    assignGradient(idxG3, QColor(200, 220, 255), QColor(0, 60, 170));
    if (!idxAxon.empty()) {
        std::vector<QColor> tmp(idxAxon.size());
        generateRandomColors(idxAxon.size(), tmp);
        for (size_t k = 0; k < idxAxon.size(); ++k) colors[idxAxon[k]] = tmp[k];
    }

    if (x.size() != y.size() || x.size() != z.size() || x.size() != radius.size() || x.size() != groupIds.size()) {
        qDebug() << "Input size mismatch";
        return;
    }

    for (size_t i = 0; i < x.size(); ++i) {
        if (x[i].size() != y[i].size() || x[i].size() != z[i].size() || x[i].size() != radius[i].size()) {
            qDebug() << "Bundle size mismatch at" << i;
            return;
        }
    }


    // --- 2. Pass 1: Sum and Count (Fast single-thread) ---
    double totalSumX = 0, totalSumY = 0, totalSumZ = 0;
    long long totalCount = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        totalCount += x[i].size();
        for (const auto& val : x[i]) totalSumX += val;
        for (const auto& val : y[i]) totalSumY += val;
        for (const auto& val : z[i]) totalSumZ += val;
    }

    if (totalCount == 0) return;

    avgX = totalSumX / totalCount;
    avgY = totalSumY / totalCount;
    avgZ = totalSumZ / totalCount;
    QVector3D avg(avgX, avgY, avgZ);


    // --- 3. Pass 2: Multithreaded Processing ---
    // Reserve memory so vectors don't reallocate
    spherePositions.resize(totalCount);
    sphereRadii.resize(totalCount);
    axonColors.resize(totalCount);

    // Verify total count matches reserved size
    size_t verifyTotal = 0;
    for (size_t i = 0; i < x.size(); ++i) verifyTotal += x[i].size();
    Q_ASSERT(verifyTotal == (size_t)totalCount);
    Q_ASSERT(spherePositions.size() == (size_t)totalCount);

    unsigned int numThreads = std::thread::hardware_concurrency();
    if (numThreads == 0) numThreads = 2; 

    std::vector<std::future<void>> futures;
    size_t bundlesPerThread = (x.size() + numThreads - 1) / numThreads;

    // We need to track where each thread starts writing in the flat vectors
    std::vector<size_t> threadOffsets(numThreads, 0);
    size_t currentOffset = 0;
    for (unsigned int t = 0; t < numThreads; ++t) {
        threadOffsets[t] = currentOffset;
        size_t start = t * bundlesPerThread;
        size_t end = std::min(start + bundlesPerThread, x.size());
        for (size_t i = start; i < end; ++i) currentOffset += x[i].size();
    }

    for (unsigned int t = 0; t < numThreads; ++t) {
        size_t startBundle = t * bundlesPerThread;
        size_t endBundle = std::min(startBundle + bundlesPerThread, x.size());
        size_t writePos = threadOffsets[t];

        if (startBundle >= endBundle) break;

        futures.push_back(std::async(std::launch::async, [=, &x, &y, &z, &radius, &colors, &avg]() {
            size_t localWrite = writePos;
            for (size_t i = startBundle; i < endBundle; ++i) {
                for (size_t j = 0; j < x[i].size(); ++j) {
                    this->spherePositions[localWrite] = QVector3D(x[i][j], y[i][j], z[i][j]) - avg;
                    this->sphereRadii[localWrite] = static_cast<float>(radius[i][j]);
                    this->axonColors[localWrite] = colors[i];
                    localWrite++;
                }
            }
        }));
    }

    for (auto& f : futures) f.get();

    // Finalize
    initialspherePositions = spherePositions;
    initialsphereRadii = sphereRadii;
    initialaxonColors = axonColors;

    // Ensure sphere template exists even if GL init timing was off
    if (sphereVertices.empty()) {

        generateSphereVBO(30, 30, 1.0f);
    }
    buildBatchedGeometry();

    update();
}

void OpenGLWindow::setVoxelBounds(const QVector3D& minCorner, const QVector3D& maxCorner)
{
    voxelMinCorner = minCorner;
    voxelMaxCorner = maxCorner;
    hasVoxelBounds = true;
    update();
}

void OpenGLWindow::resetCamera()
{
    orbitTheta = 0.0f;
    orbitPhi = 0.0f;
    cameraDistance = 100.0f; // Pull the camera back to a safe distance
    zoomFactor = 1.0f;       // Reset the scroll wheel multiplier
    update();
}

void OpenGLWindow::initializeGL()
{
    initializeOpenGLFunctions();
    glEnable(GL_DEPTH_TEST);
    
    // THIS ALLOWS THE SHADER TO CONTROL POINT SIZES
    glEnable(GL_PROGRAM_POINT_SIZE); 

    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);

    shaderProgram = new QOpenGLShaderProgram(this);

    const char *vertexShaderSource =
        "#version 330 core\n"
        "layout(location = 0) in vec3 position;\n"
        "layout(location = 1) in vec3 color;\n"
        "layout(location = 2) in float radius;\n"
        "out vec3 fragColor;\n"
        "uniform mat4 mvp;\n"            // Model-View-Projection matrix
        "uniform float viewportHeight;\n" // Viewport height in pixels
        "uniform float tanHalfFovY;\n"    // tan(fovY / 2), matches projectionMatrix's fovY
        "void main() {\n"
        "    gl_Position = mvp * vec4(position, 1.0);\n"
        "    fragColor = color;\n"
        // Physically-based point-sprite size: the on-screen diameter of a sphere of
        // world radius `radius` at (view-space) distance gl_Position.w, so that two
        // spheres touching in world space (center distance == r1+r2) still visibly
        // touch on screen regardless of how different their radii are (e.g. a large
        // soma vs. a thin process) -- an ad hoc size formula here previously caused
        // touching spheres of very different radii to render with a visible gap.
        "    float calculatedSize = (radius * viewportHeight) / (gl_Position.w * tanHalfFovY);\n"
        "    gl_PointSize = clamp(calculatedSize, 1.0, 500.0);\n"
        "}\n";

    const char *fragmentShaderSource =
        "#version 330 core\n"
        "in vec3 fragColor;\n"
        "out vec4 finalColor;\n"
        "void main() {\n"
        //   This bit of math turns the flat square point into a perfect circle
        "    vec2 coord = gl_PointCoord - vec2(0.5);\n"
        "    if(length(coord) > 0.5) discard;\n"
        "    finalColor = vec4(fragColor, 1.0);\n"
        "}\n";

    shaderProgram->addShaderFromSourceCode(QOpenGLShader::Vertex, vertexShaderSource);
    shaderProgram->addShaderFromSourceCode(QOpenGLShader::Fragment, fragmentShaderSource);
    shaderProgram->link();

    // Simple flat-color line shader, used only to draw the voxel wireframe box.
    const char *lineVertexShaderSource =
        "#version 330 core\n"
        "layout(location = 0) in vec3 position;\n"
        "uniform mat4 mvp;\n"
        "void main() {\n"
        "    gl_Position = mvp * vec4(position, 1.0);\n"
        "}\n";

    const char *lineFragmentShaderSource =
        "#version 330 core\n"
        "uniform vec3 lineColor;\n"
        "out vec4 finalColor;\n"
        "void main() {\n"
        "    finalColor = vec4(lineColor, 1.0);\n"
        "}\n";

    lineShaderProgram = new QOpenGLShaderProgram(this);
    lineShaderProgram->addShaderFromSourceCode(QOpenGLShader::Vertex, lineVertexShaderSource);
    lineShaderProgram->addShaderFromSourceCode(QOpenGLShader::Fragment, lineFragmentShaderSource);
    lineShaderProgram->link();
}

void OpenGLWindow::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h);
    projectionMatrix.setToIdentity();
    // CHANGED: 1.0f and 10000.0f to restore Z-buffer precision
    projectionMatrix.perspective(fovYDegrees, float(w) / float(h), 10.0f, 10000.0f);
    viewportHeightPx = float(h);
}
QVector3D rotateAround(const QVector3D& position, float deltaTheta, float deltaPhi) {
    float x = position.x();
    float y = position.y();
    float z = position.z();  // Z remains the same after rotation around the z-axis
    float newX = x; 
    float newY = y; 
    float newZ= z;

    if (deltaTheta == 0.0f && deltaPhi == 0.0f) {
        return position;
    }

    if (deltaTheta != 0.0f){
        // Calculate the new x and y coordinates after the rotation around z axis
        newX = x * cos(deltaTheta) - y * sin(deltaTheta);
        newY = x * sin(deltaTheta) + y * cos(deltaTheta);

        x = newX;
        y = newY;
    }
    if (deltaPhi != 0.0f){
        // Calculate the new x and z coordinates after the rotation around x axis
        newX = x * cos(deltaPhi) - z * sin(deltaPhi);
        newZ = x * sin(deltaPhi) + z * cos(deltaPhi);
    }

    // Return the new position vector
    return QVector3D(newX, newY, newZ);
}
void OpenGLWindow::paintGL()
{
    // 1. Safety checks and window title update
    if (!geometryReady || totalBatchedVertices == 0) return;
    
    this->setTitle(QString("Cells: %1 | Cam Zoom: %2")
                   .arg(totalBatchedVertices)
                   .arg(cameraDistance * zoomFactor));

    // 2. Clear the screen and enable depth testing
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    
    // Crucial: allow the vertex shader to dictate point sizes
    glEnable(GL_PROGRAM_POINT_SIZE); 

    // 3. Calculate Camera Position (Orbit Math)
    QVector3D targetCenter(0.0f, 0.0f, 0.0f);
    float currentRadius = cameraDistance * zoomFactor;
    float radTheta = qDegreesToRadians(orbitTheta);
    float radPhi = qDegreesToRadians(orbitPhi);

    cameraPosition = targetCenter + QVector3D(
        currentRadius * cos(radPhi) * sin(radTheta),
        currentRadius * sin(radPhi),
        currentRadius * cos(radPhi) * cos(radTheta)
    );

    // 4. Calculate Matrices
    QMatrix4x4 viewMatrix;
    viewMatrix.lookAt(cameraPosition, targetCenter, QVector3D(0, 1, 0));
    
    // Combine Projection and View into a single Model-View-Projection matrix
    QMatrix4x4 mvpMatrix = projectionMatrix * viewMatrix;

    // 5. Render using the Shader
    if (!shaderProgram || !shaderProgram->isLinked()) {
        qWarning() << "Shader program is not ready!";
        return;
    }

    shaderProgram->bind();

    // -- Set Uniforms (Global variables for this draw call) --
    shaderProgram->setUniformValue("mvp", mvpMatrix);
    shaderProgram->setUniformValue("viewportHeight", viewportHeightPx);
    shaderProgram->setUniformValue("tanHalfFovY", static_cast<float>(qTan(qDegreesToRadians(fovYDegrees / 2.0f))));

    // -- Set Attributes (Per-point data variables) --
    // Location 0: Positions (X, Y, Z)
    shaderProgram->enableAttributeArray(0);
    shaderProgram->setAttributeArray(0, GL_FLOAT, batchedVertices.data(), 3);

    // Location 1: Colors (R, G, B)
    shaderProgram->enableAttributeArray(1);
    shaderProgram->setAttributeArray(1, GL_FLOAT, batchedColors.data(), 3);

    // Location 2: Radii (Float)
    shaderProgram->enableAttributeArray(2);
    shaderProgram->setAttributeArray(2, GL_FLOAT, batchedRadii.data(), 1);

    // 6. Execute the Draw Call! (1 Vertex = 1 Sphere)
    glDrawArrays(GL_POINTS, 0, totalBatchedVertices);

    // 7. Cleanup state for the next frame
    shaderProgram->disableAttributeArray(0);
    shaderProgram->disableAttributeArray(1);
    shaderProgram->disableAttributeArray(2);
    shaderProgram->release();

    // 8. Draw the voxel boundary wireframe, so tips near a wall can be judged
    // against the actual simulation box instead of guessed at.
    if (hasVoxelBounds) {
        drawVoxelWireframe();
    }
}

void OpenGLWindow::drawVoxelWireframe()
{
    if (!lineShaderProgram || !lineShaderProgram->isLinked()) {
        return;
    }

    // Recenter the box corners by the same average offset applied to the
    // sphere positions in setSpheres(), so the wireframe lines up with them.
    QVector3D avg(avgX, avgY, avgZ);
    QVector3D lo = voxelMinCorner - avg;
    QVector3D hi = voxelMaxCorner - avg;

    QVector3D corners[8];
    for (int c = 0; c < 8; ++c) {
        corners[c] = QVector3D(
            (c & 1) ? hi.x() : lo.x(),
            (c & 2) ? hi.y() : lo.y(),
            (c & 4) ? hi.z() : lo.z()
        );
    }

    // Two corners are connected by an edge iff they differ in exactly one bit.
    std::vector<GLfloat> lineVertices;
    lineVertices.reserve(12 * 2 * 3);
    for (int c = 0; c < 8; ++c) {
        for (int bit = 1; bit <= 4; bit <<= 1) {
            if (c & bit) continue; // only emit each edge once
            const QVector3D &a = corners[c];
            const QVector3D &b = corners[c | bit];
            lineVertices.push_back(a.x()); lineVertices.push_back(a.y()); lineVertices.push_back(a.z());
            lineVertices.push_back(b.x()); lineVertices.push_back(b.y()); lineVertices.push_back(b.z());
        }
    }

    QMatrix4x4 viewMatrix;
    viewMatrix.lookAt(cameraPosition, QVector3D(0.0f, 0.0f, 0.0f), QVector3D(0, 1, 0));
    QMatrix4x4 mvpMatrix = projectionMatrix * viewMatrix;

    lineShaderProgram->bind();
    lineShaderProgram->setUniformValue("mvp", mvpMatrix);
    lineShaderProgram->setUniformValue("lineColor", QVector3D(1.0f, 1.0f, 0.0f)); // yellow

    lineShaderProgram->enableAttributeArray(0);
    lineShaderProgram->setAttributeArray(0, GL_FLOAT, lineVertices.data(), 3);

    glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(lineVertices.size() / 3));

    lineShaderProgram->disableAttributeArray(0);
    lineShaderProgram->release();
}

bool OpenGLWindow::isSphereInFrustum(const QVector3D& pos, float radius) {
    // 1. Get the direction vector from the camera to the sphere
    QVector3D toSphere = pos - cameraPosition;
    
    // 2. Simple distance check: if the sphere is too far behind the camera, 
    // or too far away from the center, skip it.
    // (This is a simplified distance-based culling)
    float dist = toSphere.length();
    
    // If it's outside our render distance, don't draw
    if (dist > (cameraDistance * zoomFactor * 10.0f)) return false;
    
    return true;
}
void OpenGLWindow::drawSphere(const QVector3D& position, float radius, const QColor& color)
{
    // Translate the model to the sphere's position
    glPushMatrix();
    glTranslatef(position.x(), position.y(), position.z());

    // Scale the sphere (optional if your sphere's vertices are not already scaled)
    glScalef(radius, radius, radius);

    // Set the sphere's color
    glColor3f(color.redF(), color.greenF(), color.blueF());

    // Enable vertex array and normal array
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

    // Bind the VBO for the sphere's vertices
    glBindBuffer(GL_ARRAY_BUFFER, sphereVBO);
    glVertexPointer(3, GL_FLOAT, 0, 0);  // Point to the vertex data

    // Bind the VBO for the sphere's normals (if using lighting)
    // Assuming normals are stored in the same order as vertices
    glNormalPointer(GL_FLOAT, 0, &sphereNormals[0]);

    // Draw the sphere using the vertices from the VBO
    glDrawArrays(GL_TRIANGLE_STRIP, 0, sphereVertices.size() / 3);

    // Disable the arrays and unbind the VBO
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glPopMatrix();
}


// Mouse press event - store the initial position and track which button is pressed
void OpenGLWindow::mousePressEvent(QMouseEvent *event)
{
    lastMousePosition = event->pos();

    if (event->button() == Qt::LeftButton) {
        leftMousePressed = true;  // Left mouse button pressed for rotation
    } else if (event->button() == Qt::RightButton) {
        rightMousePressed = true;  // Right mouse button pressed for panning
    }
}
void OpenGLWindow::mouseMoveEvent(QMouseEvent *event)
{
    int deltaX = event->x() - lastMousePosition.x();
    int deltaY = event->y() - lastMousePosition.y();

    float sensitivity = 0.3f;  // Control speed of rotation

    if (leftMousePressed) {
        // Adjust spherical angles
        orbitTheta += deltaX * sensitivity;
        orbitPhi += deltaY * sensitivity;

        // Clamp phi to prevent flipping (avoid going over poles)
        orbitPhi = std::max(-89.0f, std::min(89.0f, orbitPhi));
    }

    lastMousePosition = event->pos();
    update();
}


// Mouse release event - reset the state of the mouse buttons
void OpenGLWindow::mouseReleaseEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton) {
        leftMousePressed = false;
        SphererotationX = 0.0f;
        SphererotationY = 0.0f;
    } else if (event->button() == Qt::RightButton) {
        rightMousePressed = false;
    }
}

void OpenGLWindow::wheelEvent(QWheelEvent *event)
{
    // If scrolling forward, zoom in (decrease radius)
    if (event->angleDelta().y() > 0) {
        zoomFactor *= 0.9f; 
    } else {
        zoomFactor *= 1.1f; 
    }
    
    // THIS LINE IS REQUIRED TO ACTUALLY DRAW THE ZOOM
    update(); 
}

void OpenGLWindow::generateRandomColors(int count, std::vector<QColor>& colors) {
    // Use a static RNG so successive calls don't reseed to the same sequence
    static thread_local std::mt19937 rng(std::random_device{}());

    // Ensure we overwrite exactly 'count' entries
    colors.clear();
    colors.reserve(count);

    // Avoid too-dark shades; force full alpha
    std::uniform_int_distribution<int> dist(32, 255); // min 32 to dodge near-black

    for (int i = 0; i < count; ++i) {
        colors.emplace_back(dist(rng), dist(rng), dist(rng), 255);
    }
}

void OpenGLWindow::generateSphereVBO(int slices, int stacks, float radius)
{
    sphereVertices.clear();
    sphereNormals.clear();

    for (int i = 0; i <= stacks; ++i) {
        // FIX: use (i) not (i-1) for lat0 to avoid the off-by-one
        double lat0 = M_PI * (-0.5 + double(i) / stacks);
        double z0   = qSin(lat0) * radius;
        double zr0  = qCos(lat0) * radius;

        double lat1 = M_PI * (-0.5 + double(i + 1) / stacks);
        double z1   = qSin(lat1) * radius;
        double zr1  = qCos(lat1) * radius;

        for (int j = 0; j <= slices; ++j) {
            double lng = 2 * M_PI * double(j) / slices; // FIX: j not (j-1)
            double x = qCos(lng);
            double y = qSin(lng);

            sphereVertices.push_back(x * zr0);
            sphereVertices.push_back(y * zr0);
            sphereVertices.push_back(z0);

            sphereVertices.push_back(x * zr1);
            sphereVertices.push_back(y * zr1);
            sphereVertices.push_back(z1);
        }
    }

    // GPU upload (only valid inside GL context)
    if (QOpenGLContext::currentContext() != nullptr) {
        if (sphereVBO) glDeleteBuffers(1, &sphereVBO);
        glGenBuffers(1, &sphereVBO);
        glBindBuffer(GL_ARRAY_BUFFER, sphereVBO);
        glBufferData(GL_ARRAY_BUFFER,
                     sphereVertices.size() * sizeof(float),
                     sphereVertices.data(),
                     GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    qDebug() << "Sphere template built:" << sphereVertices.size() << "floats";
}