#ifndef OPENGLWINDOW_H
#define OPENGLWINDOW_H

#include <QOpenGLWindow>
#include <QOpenGLFunctions>
#include <QMatrix4x4>
#include <QVector3D>
#include <QTimer>
#include <QColor>
#include <vector>

class OpenGLWindow : public QOpenGLWindow, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    explicit OpenGLWindow(QWindow *parent = nullptr);
    ~OpenGLWindow();

    void setSpheres(const std::vector<std::vector<double>>& x,
                    const std::vector<std::vector<double>>& y,
                    const std::vector<std::vector<double>>& z,
                    const std::vector<std::vector<double>>& radius,
                    const std::vector<int>& groupIds);

    void resetCamera();
    enum class SphereGroup : int { Axon = 0, Glial1 = 1, Glial2 = 2, Blood = 3 };

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;

    // Event handlers
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;
    void generateSphereVBO(int slices, int stacks, float radius);

private:
    std::vector<QVector3D> spherePositions;
    std::vector<double> sphereRadii;
    std::vector<QColor> axonColors;  // Store colors for each axon

    std::vector<QVector3D> initialspherePositions;
    std::vector<double> initialsphereRadii;
    std::vector<QColor> initialaxonColors;  // Store colors for each axon

    double avgX;
    double avgY;
    double avgZ;
    double maxX;
    double maxY;
    double maxZ;

    float orbitTheta;
    float orbitPhi;
    float cameraDistance;
    


    QMatrix4x4 projectionMatrix;
    QTimer timer;
    float zoomFactor = 1.0f;  // Adjust this for zoom control
    float rotationX = 0.0f;   // X-axis rotation
    float rotationY = 0.0f;   // Y-axis rotation
    QVector3D cameraPosition;  // Camera position for panning

    float SphererotationX = 0.0f;   // X-axis rotation
    float SphererotationY = 0.0f;   // Y-axis rotation

    QPoint lastMousePosition;  // Last mouse position for tracking movement
    bool leftMousePressed = false;  // Flag to track if the left mouse button is pressed
    bool rightMousePressed = false;  // Flag to track if the right mouse button is pressed

    void drawSphere(const QVector3D &position, float radius, const QColor &color);
    void generateRandomColors(int count, std::vector<QColor> &colors);

    // sphere generation
    GLuint sphereVBO;
    std::vector<float> sphereVertices;
    std::vector<float> sphereNormals;
    
};

#endif // OPENGLWINDOW_H

