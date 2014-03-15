#ifndef MYDRAWAREA_H
#define MYDRAWAREA_H

#include <qglviewer/qglviewer.h>
#include "SurfaceMeshModel.h"

extern QString curLabel;
extern QStringList AllLabels;
extern QVector<QColor> UniqueColors;
extern QVector<QString> labelNames;

static QString meshOps[] = {"rotateleft", "rotateright", "flip", "invertpart"};
enum MeshOperation{ ROTATE_LEFT, ROTATE_RIGHT, FLIP, INVERT_PART };
extern MeshOperation curOp;

class MyDrawArea : public QGLViewer
{
public:
    MyDrawArea(SurfaceMesh::SurfaceMeshModel * mesh, QString filename) : m(mesh), filename(filename){}
	~MyDrawArea();

	void draw();

	void postSelection(const QPoint& point);
	void mousePressEvent(QMouseEvent * event);

    SurfaceMesh::SurfaceMeshModel * m;
    QString filename;
};

#endif // MYDRAWAREA_H
