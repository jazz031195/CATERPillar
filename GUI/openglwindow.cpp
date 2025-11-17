#include "openglwindow.h"
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QtMath>
#include <QWheelEvent>
#include <QMouseEvent>
#include <random>

OpenGLWindow::OpenGLWindow(QWindow *parent)
    : QOpenGLWindow(NoPartialUpdate, parent),
      orbitTheta(0.0f), orbitPhi(0.0f), cameraDistance(50.0f)
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



void OpenGLWindow::setSpheres(const std::vector<std::vector<double>>& x,
                              const std::vector<std::vector<double>>& y,
                              const std::vector<std::vector<double>>& z,
                              const std::vector<std::vector<double>>& radius,
                              const std::vector<int>& groupIds)
{
    spherePositions.clear();
    sphereRadii.clear();
    axonColors.clear();
    initialspherePositions.clear();
    initialsphereRadii.clear();
    initialaxonColors.clear();

    // --- Build color per-bundle using groupIds ---
    Q_ASSERT(static_cast<int>(x.size()) == static_cast<int>(groupIds.size()));
    std::vector<QColor> colors(x.size());

    // Collect bundle indices per group to make smooth gradients within each group
    std::vector<size_t> idxAxon, idxG1, idxG2, idxBV;
    for (size_t i = 0; i < groupIds.size(); ++i) {
        switch (static_cast<SphereGroup>(groupIds[i])) {
            case SphereGroup::Axon:   idxAxon.push_back(i); break;
            case SphereGroup::Glial1: idxG1.push_back(i);   break;
            case SphereGroup::Glial2: idxG2.push_back(i);   break;
            case SphereGroup::Blood:  idxBV.push_back(i);   break;
        }
    }

    auto assignGradient = [&](const std::vector<size_t>& idx,
                              const QColor& c0, const QColor& c1)
    {
        if (idx.empty()) return;
        if (idx.size() == 1) {
            colors[idx[0]] = c1;
            return;
        }
        for (size_t k = 0; k < idx.size(); ++k) {
            double t = static_cast<double>(k) / static_cast<double>(idx.size() - 1);
            colors[idx[k]] = lerp(c0, c1, t);
        }
    };

    // Palettes:
    // Blood vessels: light red -> dark red
    assignGradient(idxBV, QColor(255, 200, 200), QColor(170, 0, 0));
    // Glial pop1: light green -> mid green
    assignGradient(idxG1, QColor(200, 255, 200), QColor(0, 150, 0));
    // Glial pop2: lighter to darker/different green to distinguish
    assignGradient(idxG2, QColor(220, 255, 220), QColor(0, 110, 0));
    // Axons: keep random or choose a neutral palette
    if (!idxAxon.empty()) {
        std::vector<QColor> tmp(idxAxon.size());
        generateRandomColors(idxAxon.size(), tmp); // your existing helper
        for (size_t k = 0; k < idxAxon.size(); ++k) colors[idxAxon[k]] = tmp[k];
    }

    // --- Compute averages and populate buffers ---
    double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
    long long count = 0;

    for (size_t i = 0; i < x.size(); ++i) {
        for (size_t j = 0; j < x[i].size(); ++j) {
            sumX += x[i][j];
            sumY += y[i][j];
            sumZ += z[i][j];
            ++count;
        }
    }
    if (count == 0) return; // nothing to draw, avoid div/0

    avgX = sumX / static_cast<double>(count);
    avgY = sumY / static_cast<double>(count);
    avgZ = sumZ / static_cast<double>(count);
    QVector3D avg(avgX, avgY, avgZ);

    maxX = std::numeric_limits<double>::lowest();
    maxY = std::numeric_limits<double>::lowest();
    maxZ = std::numeric_limits<double>::lowest();

    for (size_t i = 0; i < x.size(); ++i) {
        for (size_t j = 0; j < x[i].size(); ++j) {
            QVector3D pos(x[i][j], y[i][j], z[i][j]);
            spherePositions.push_back(pos - avg);
            sphereRadii.push_back(radius[i][j]);
            axonColors.push_back(colors[i]);  // color per bundle

            initialspherePositions.push_back(pos - avg);
            initialsphereRadii.push_back(radius[i][j]);
            initialaxonColors.push_back(colors[i]);

            maxX = std::max(maxX, x[i][j]);
            maxY = std::max(maxY, y[i][j]);
            maxZ = std::max(maxZ, z[i][j]);
        }
    }

    // Reset transforms/camera
    SphererotationX = 0.0f;
    SphererotationY = 0.0f;
    cameraPosition = QVector3D(maxX, maxY, maxZ + 50.0f);

    update();
}

void OpenGLWindow::resetCamera(){

    orbitTheta= 0.0f;
    orbitPhi = 0.0f;
    cameraDistance = maxZ+ 50.0f;
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
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Convert spherical coordinates to Cartesian (for camera movement)
    float radius = cameraDistance;  // Distance from center
    float radTheta = qDegreesToRadians(orbitTheta);
    float radPhi = qDegreesToRadians(orbitPhi);

    float x = radius * cos(radPhi) * sin(radTheta);
    float y = radius * sin(radPhi);
    float z = radius * cos(radPhi) * cos(radTheta);

    cameraPosition = QVector3D(x, y, z)*zoomFactor;  // Update camera position

    // Set up view matrix
    QMatrix4x4 modelViewMatrix;
    modelViewMatrix.setToIdentity();
    modelViewMatrix.lookAt(cameraPosition, QVector3D(0, 0, 0), QVector3D(0, 1, 0));

    QMatrix4x4 mvpMatrix = projectionMatrix * modelViewMatrix;
    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(mvpMatrix.constData());

    // Sort and draw spheres
    std::vector<std::pair<double, size_t>> distances;
    for (size_t i = 0; i < spherePositions.size(); ++i) {
        QVector3D diff = spherePositions[i] - cameraPosition;
        distances.push_back(std::make_pair(diff.lengthSquared(), i));
    }

    std::sort(distances.begin(), distances.end(), std::greater<std::pair<double, size_t>>());

    for (size_t i = 0; i < spherePositions.size(); ++i) {
        size_t index = distances[i].second;
        drawSphere(spherePositions[index], sphereRadii[index], axonColors[index]);
    }
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

// Zoom in or out depending on the scroll direction
void OpenGLWindow::wheelEvent(QWheelEvent *event)
{
    if (event->angleDelta().y() > 0) {
        zoomFactor *= 1.1f;  // Zoom in
    } else {
        zoomFactor /= 1.1f;  // Zoom out
    }

    update();  // Repaint with new zoom
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
        double lat0 = M_PI * (-0.5 + double(i - 1) / stacks);
        double z0 = qSin(lat0) * radius;
        double zr0 = qCos(lat0) * radius;

        double lat1 = M_PI * (-0.5 + double(i) / stacks);
        double z1 = qSin(lat1) * radius;
        double zr1 = qCos(lat1) * radius;

        for (int j = 0; j <= slices; ++j) {
            double lng = 2 * M_PI * double(j - 1) / slices;
            double x = qCos(lng);
            double y = qSin(lng);

            // Vertex and normal for the first point
            sphereVertices.push_back(x * zr0);
            sphereVertices.push_back(y * zr0);
            sphereVertices.push_back(z0);

            sphereNormals.push_back(x * zr0);
            sphereNormals.push_back(y * zr0);
            sphereNormals.push_back(z0);

            // Vertex and normal for the second point
            sphereVertices.push_back(x * zr1);
            sphereVertices.push_back(y * zr1);
            sphereVertices.push_back(z1);

            sphereNormals.push_back(x * zr1);
            sphereNormals.push_back(y * zr1);
            sphereNormals.push_back(z1);
        }
    }

    // Generate and bind VBO
    glGenBuffers(1, &sphereVBO);
    glBindBuffer(GL_ARRAY_BUFFER, sphereVBO);
    glBufferData(GL_ARRAY_BUFFER, sphereVertices.size() * sizeof(float), &sphereVertices[0], GL_STATIC_DRAW);
}