#pragma once

#include "SurfaceMeshPlugins.h"
#include "RichParameterSet.h"

class voxelize: public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "voxelize.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Voxelize mesh"; }
    QString description() { return "Voxelize mesh"; }

    void initParameters(RichParameterSet* pars);
    void applyFilter(RichParameterSet* pars);
};
