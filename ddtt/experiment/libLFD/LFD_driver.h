#pragma once
#include "SurfaceMeshPlugins.h"
class lfd_driver : public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "lfd_driver.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "LFD method"; }
    QString description() { return "Computes shape descriptor"; }
    void applyFilter(RichParameterSet*);
};
