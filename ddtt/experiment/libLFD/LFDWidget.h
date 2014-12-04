#pragma once

#include <QWidget>
#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include "qglviewer/camera.h"

#include "LFD.h"
#include "GL/GLU.h"

namespace Ui {
class LFDWidget;
}

template <typename Scalar>
struct GLVertex{
	Scalar x, y, z;
	Scalar nx, ny, nz;
	GLVertex(Scalar X = 0, Scalar Y = 0, Scalar Z = 0, Scalar nX = 0, Scalar nY = 0, Scalar nZ = 0) : x(X), y(Y), z(Z), nx(nX), ny(nY), nz(nZ) {}
};
typedef GLVertex<double> GLVertexd;

template <typename Index>
struct GLFace{
	Index v1, v2, v3;
	GLFace(Index v1, Index v2, Index v3) :v1(v1), v2(v2), v3(v3){}
};
typedef GLFace<unsigned int> GLFaceui;

class MyGLWidget : public QOpenGLWidget, protected QOpenGLFunctions{
public:
	MyGLWidget(SurfaceMesh::SurfaceMeshModel *model, QWidget * parent) : model(model), QOpenGLWidget(parent){ isReady = false; isDone = false; }

	bool isReady;
	bool isDone;

protected:
    void initializeGL();
    void paintGL();
    SurfaceMesh::SurfaceMeshModel * model;
};

struct RenderOption{
	qglviewer::Vec cameraPos;
	bool isColor;
	bool isBackground;
	bool isLights;
	bool isShowCameras;
	RenderOption(const qglviewer::Vec & cameraPos = qglviewer::Vec(1, 1, 1)) : cameraPos(cameraPos) { isColor = isBackground = isLights = isShowCameras = true; }
};

class LFDWidget : public QWidget
{
    Q_OBJECT

public:
    explicit LFDWidget(SurfaceMesh::SurfaceMeshModel *model, QWidget *parent = 0);
    ~LFDWidget();

    Ui::LFDWidget *ui;
};
