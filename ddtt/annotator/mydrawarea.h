#ifndef MYDRAWAREA_H
#define MYDRAWAREA_H

#include <qglviewer/qglviewer.h>
#include "RenderObjectExt.h"

#include "StructureGraph.h"

extern QString curLabel;
extern QStringList AllLabels;
extern QVector<QColor> LabelColors;
extern QVector<QString> labelNames;

static QString shapeOps[] = {"none", "rotateleft", "rotateright", "labelpart"};

enum ShapeOperation{ NONE_OP, ROTATE_LEFT, ROTATE_RIGHT,
					 LABEL_PART};

extern ShapeOperation curOp;

extern QString curLabel;
extern int curLabelIdx;

extern QVector<QColor> LabelColors;
extern QVector<QColor> ParentColors;
extern QMap < QString, int > labelIndices;

class MyDrawArea : public QGLViewer
{
	Q_OBJECT
public:
    MyDrawArea(QSharedPointer<Structure::Graph> mesh, QString filename);
	~MyDrawArea();

	void draw();

	void postSelection(const QPoint& point);
	void mousePressEvent(QMouseEvent * event);
	void focusInEvent(QFocusEvent * event);

    QSharedPointer<Structure::Graph> m;
    QString filename, basename;
	bool isDeleted;

	bool isDrawWireframe;
	bool isDoubleLight;
	QColor bg, fg;

	// Debug:
	QVector<RenderObject::Base*> debugItems;

signals:
	void gotFocus(MyDrawArea *);
};

extern MyDrawArea * lastSelected;

#endif // MYDRAWAREA_H
