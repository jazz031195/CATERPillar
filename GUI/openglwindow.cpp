#include "openglwindow.h"
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


    qDebug() << "buildBatchedGeometry called";
    qDebug() << "spherePositions.size():" << spherePositions.size();
    qDebug() << "sphereRadii.size():" << sphereRadii.size();
    qDebug() << "axonColors.size():" << axonColors.size();
    qDebug() << "sphereVertices.size():" << sphereVertices.size();
    geometryReady = false;


    batchedVertices.clear();
    batchedColors.clear();

    size_t numSpheres = spherePositions.size();
    if (numSpheres == 0 || sphereVertices.empty()) return;
    
    size_t floatsPerSphere = sphereVertices.size(); 
    batchedVertices.reserve(numSpheres * floatsPerSphere);
    batchedColors.reserve(numSpheres * floatsPerSphere);

    for (size_t i = 0; i < numSpheres; ++i) {
        const QVector3D& pos = spherePositions[i];
        float r = sphereRadii[i];
        const QColor& col = axonColors[i];

        size_t vertCount = (sphereVertices.size() / 3) * 3; // truncate to safe multiple of 3
        for (size_t v = 0; v + 2 < vertCount; v += 3) {
            batchedVertices.push_back((sphereVertices[v] * r) + pos.x());
            batchedVertices.push_back((sphereVertices[v+1] * r) + pos.y());
            batchedVertices.push_back((sphereVertices[v+2] * r) + pos.z());

            batchedColors.push_back(col.redF());
            batchedColors.push_back(col.greenF());
            batchedColors.push_back(col.blueF());
        }
    }
    totalBatchedVertices = batchedVertices.size() / 3;
    qDebug() << "SUCCESS: Generated" << totalBatchedVertices << "vertices.";
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
    std::vector<size_t> idxAxon, idxG1, idxG2, idxBV;
    for (size_t i = 0; i < groupIds.size(); ++i) {
        switch (static_cast<SphereGroup>(groupIds[i])) {
            case SphereGroup::Axon:   idxAxon.push_back(i); break;
            case SphereGroup::Glial1: idxG1.push_back(i);   break;
            case SphereGroup::Glial2: idxG2.push_back(i);   break;
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
    glEnable(GL_DEPTH_TEST);  // Enable depth test
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);  // Set a background color

    // Generate the VBO for the sphere geometry
    generateSphereVBO(30, 30, 1.0f);  // Example: 30 slices, 30 stacks, radius 1.0

}

void OpenGLWindow::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h);  // Update viewport to match new window size
    projectionMatrix.setToIdentity();
    projectionMatrix.perspective(45.0f, float(w) / float(h), 0.1f, 100000.0f);  // Perspective projection
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

    if (!geometryReady || totalBatchedVertices == 0) return;
    this->setTitle(QString("Vertices: %1 | Cam: %2")
                   .arg(totalBatchedVertices)
                   .arg(cameraDistance * zoomFactor));

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

    if (totalBatchedVertices == 0) return;

    QVector3D targetCenter(0.0f, 0.0f, 0.0f);
    float radius = cameraDistance * zoomFactor;
    float radTheta = qDegreesToRadians(orbitTheta);
    float radPhi = qDegreesToRadians(orbitPhi);

    cameraPosition = targetCenter + QVector3D(
        radius * cos(radPhi) * sin(radTheta),
        radius * sin(radPhi),
        radius * cos(radPhi) * cos(radTheta)
    );

    QMatrix4x4 modelViewMatrix;
    modelViewMatrix.setToIdentity();
    modelViewMatrix.lookAt(cameraPosition, targetCenter, QVector3D(0, 1, 0));

    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf(projectionMatrix.constData());
    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(modelViewMatrix.constData());

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);

    glVertexPointer(3, GL_FLOAT, 0, batchedVertices.data());
    glColorPointer(3, GL_FLOAT, 0, batchedColors.data());

    // Debug mode
    glPointSize(2.0f);
    glDrawArrays(GL_POINTS, 0, totalBatchedVertices);

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);


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