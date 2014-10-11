#pragma once

#include "SurfaceMeshPlugins.h"
#include "RichParameterSet.h"
#include "interfaces/ModePluginDockWidget.h"

#include "deform-handle.h"

// ShapeOp
#include "Solver.h"

class deform: public SurfaceMeshModePlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "deform.plugin.starlab")
    Q_INTERFACES(ModePlugin)

    QIcon icon(){ return QIcon(":/images/deform.png"); }

public:
	deform() : widget(NULL), dockwidget(NULL), isDeformReady(false), solver(NULL), last_selected(-1) {}

    // Main functionality
    void create();
    void destroy(){}
    void decorate();

    // Selection
    void drawWithNames();
    bool postSelection( const QPoint& );
	int last_selected;

	bool keyPressEvent(QKeyEvent*);
	bool mouseMoveEvent(QMouseEvent*);
    bool mousePressEvent(QMouseEvent* event);

    QWidget * widget;
	ModePluginDockWidget * dockwidget;

	double worldRadius;

    // Always usable
    bool isApplicable() { return true; }

    // Deformation
	void create_handle(const Vector3 & p, size_t vid);
    QVector< QSharedPointer<DeformHandle> > handles;
	bool isDeformReady;
	bool isSolving;
	ShapeOp::Solver * solver;

public slots:
	void create_ROI();
	void apply_deformation();
	
signals:
};
