#include "LFD_driver.h"
#include "LFD.h"

void lfd_driver::applyFilter(RichParameterSet *)
{
    LFD::generate((SurfaceMesh::SurfaceMeshModel*)model());
}
