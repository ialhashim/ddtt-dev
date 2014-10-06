#pragma once

#include "SurfaceMeshPlugins.h"
#include "RichParameterSet.h"
#include "interfaces/ModePluginDockWidget.h"

#include "deform-handle.h"

class deform: public SurfaceMeshModePlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "deform.plugin.starlab")
    Q_INTERFACES(ModePlugin)

    QIcon icon(){ return QIcon(":/images/deform.png"); }

public:
    deform() : widget(NULL), dockwidget(NULL) {}

    // Main functionality
    void create();
    void destroy(){}
    void decorate();

    // Selection
    void drawWithNames();
    bool postSelection( const QPoint& );

	bool keyPressEvent(QKeyEvent*);
	bool mouseMoveEvent(QMouseEvent*);
    bool mousePressEvent(QMouseEvent* event);

    QWidget * widget;
	ModePluginDockWidget * dockwidget;

    // Always usable
    bool isApplicable() { return true; }

    // Deformation
    QVector< QSharedPointer<DeformHandle> > handles;

public slots:
	void apply_deformation();
	
signals:
};
