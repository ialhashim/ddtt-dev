#include "LFD.h"
#include "LFDWidget.h"

void LFD::generate(SurfaceMesh::SurfaceMeshModel * model)
{
    LFDWidget * widget = new LFDWidget(model);
    widget->show();
}
