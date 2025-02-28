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

void OpenGLWindow::setSpheres(const std::vector<std::vector<double>>& x,
                              const std::vector<std::vector<double>>& y,
                              const std::vector<std::vector<double>>& z,
                              const std::vector<std::vector<double>>& radius)
{
    spherePositions.clear();
    sphereRadii.clear();
    axonColors.clear();

    std::vector<QColor> colors;
    generateRandomColors(x.size(), colors);

    double sumX = 0.0, sumY = 0.0 , sumZ = 0.0;
    int count = 0;

    // Populate sphere positions and calculate sums of x and y coordinates and max z
    for (size_t i = 0; i < x.size(); ++i) {
        for (size_t j = 0; j < x[i].size(); ++j) {

            sumX += x[i][j];
            sumY += y[i][j];
            sumZ += z[i][j];
            count++;
        }
    }

    maxX = std::numeric_limits<double>::lowest();
    maxY = std::numeric_limits<double>::lowest();
    maxZ = std::numeric_limits<double>::lowest();

    // Calculate average x and y
    avgX = sumX / count;
    avgY = sumY / count;
    avgZ = sumZ / count;
    QVector3D avg = QVector3D(avgX, avgY, avgZ);
    for (size_t i = 0; i < x.size(); ++i) {
        for (size_t j = 0; j < x[i].size(); ++j) {
            spherePositions.push_back(QVector3D(x[i][j], y[i][j], z[i][j]) - avg);
            sphereRadii.push_back(radius[i][j]);
            axonColors.push_back(colors[i]);

            initialspherePositions.push_back(QVector3D(x[i][j], y[i][j], z[i][j]) - avg);
            initialsphereRadii.push_back(radius[i][j]);
            initialaxonColors.push_back(colors[i]);
            maxX = std::max(maxX, x[i][j]);
            maxY = std::max(maxY, y[i][j]);
            maxZ = std::max(maxZ, z[i][j]);

        }
    }

    avgX = 0;
    avgY = 0;
    avgZ = 0;
    SphererotationX = 0.0f;
    SphererotationY = 0.0f;

    cameraPosition = QVector3D(maxX, maxY, maxZ+ 50.0f);
    update();

}

void OpenGLWindow::resetCamera(){

    orbitTheta= 0.0f;
    orbitPhi = 0.0f;
    cameraDistance = 50.0f;
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


// Generate random colors
void OpenGLWindow::generateRandomColors(int count, std::vector<QColor> &colors) {
    // Seed the random number generator with current time
    std::mt19937 rng(static_cast<unsigned int>(std::time(nullptr)));

    // Distribution for RGB values (0 to 255)
    std::uniform_int_distribution<int> dist(0, 255);

    // Generate random colors
    for (int i = 0; i < count; ++i) {
        int red = dist(rng);
        int green = dist(rng);
        int blue = dist(rng);
        colors.push_back(QColor(red, green, blue));
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