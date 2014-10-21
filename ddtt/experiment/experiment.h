#pragma once

#include "SurfaceMeshPlugins.h"
#include "interfaces/ModePluginDockWidget.h"

#include "StructureGraph.h"
#include "ShapeGraph.h"

class experiment: public SurfaceMeshModePlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "experiment.plugin.starlab")
    Q_INTERFACES(ModePlugin)

    QIcon icon(){ return QIcon(":/images/experiment.png"); }
	
public:
    QString name() { return "Deformation Correspondence"; }
    QString description() { return "Deformation Correspondence"; }

	bool isApplicable(){ return true; }
	
	void create();
    void destroy(){}
    void decorate();
	
	bool keyPressEvent(QKeyEvent*);
	bool mouseMoveEvent(QMouseEvent*);
	bool mousePressEvent(QMouseEvent* event);

	void doCorrespond();
	void doCorrespond2();

	void doCorrespondSearch();

    experiment() : widget(NULL), dockwidget(NULL) {}
    QWidget * widget;
	ModePluginDockWidget * dockwidget;

    QVector<Structure::ShapeGraph*> graphs;

public slots:
	void postCorrespond();
};
