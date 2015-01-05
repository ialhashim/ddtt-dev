#pragma once

#include "SurfaceMeshPlugins.h"
#include "interfaces/ModePluginDockWidget.h"

class SymRelationViz: public SurfaceMeshModePlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "SymRelationViz.plugin.starlab")
    Q_INTERFACES(ModePlugin)

    QIcon icon(){ return QIcon(":/images/SymRelationViz.png"); }
	
public:
    QString name() { return "SymRelationViz"; }
    QString description() { return "SymRelationViz"; }

	bool isApplicable(){ return true; }
	
	void create();
    void destroy(){}
    void decorate();

	bool keyPressEvent(QKeyEvent*);
	bool mouseMoveEvent(QMouseEvent*);
	bool mousePressEvent(QMouseEvent* event);

    SymRelationViz() : widget(NULL), dockwidget(NULL) {}
    QWidget * widget;
	ModePluginDockWidget * dockwidget;

public slots:

};
