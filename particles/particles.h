#pragma once

#include "SurfaceMeshPlugins.h"
#include "RichParameterSet.h"

class particles: public SurfaceMeshModePlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "particles.plugin.starlab")
    Q_INTERFACES(ModePlugin)

    QIcon icon(){ return QIcon(":/images/particles.png"); }

public:
    particles() { widget = NULL; }

    // Main functionality
    void create();
    void destroy(){}
    void decorate();

    QWidget * widget;

    // Always usable
    bool isApplicable() { return true; }
};
