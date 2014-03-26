#pragma once

#include <QGraphicsScene>
struct DeformationPath;

class DeformScene : public QGraphicsScene
{
public:
    DeformScene();

	void addDeformationPath( DeformationPath * path );

protected:
    void drawBackground ( QPainter * painter, const QRectF & rect );
    void drawForeground ( QPainter * painter, const QRectF & rect );
};

// Utility: Trackball code
static inline float projectOnBall(float x, float y){
        const float size       = 1.0f;
        const float size2      = size*size;
        const float size_limit = size2*0.5;
        const float d = x*x + y*y;
        return d < size_limit ? sqrt(size2 - d) : size_limit/sqrt(d);
}
#include <qglviewer/camera.h>
static inline qglviewer::Quaternion deformedBallQuaternion(int prevX, int prevY, int x, int y, float cx, float cy,
        int viewportWidth, int viewportHeight, float rotationSensitivity = 1.0){
        // Points on the deformed ball
        float px = rotationSensitivity * (prevX  - cx) / viewportWidth;
        float py = rotationSensitivity * (cy - prevY)  / viewportHeight;
        float dx = rotationSensitivity * (x - cx)          / viewportWidth;
        float dy = rotationSensitivity * (cy - y)          / viewportHeight;

        const qglviewer::Vec p1(px, py, projectOnBall(px, py));
        const qglviewer::Vec p2(dx, dy, projectOnBall(dx, dy));
        const qglviewer::Vec axis = cross(p2,p1);
        const float angle = 2.0 * asin(sqrt(axis.squaredNorm() / p1.squaredNorm() / p2.squaredNorm()));
        return qglviewer::Quaternion(axis, angle);
}
