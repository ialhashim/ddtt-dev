#include "StructureAnalysis.h"
#include "Segmentation.h"

StructureAnalysis::StructureAnalysis(ParticleMesh * pmesh) : s(pmesh)
{
    Segmentation seg( pmesh );
    for(auto d : seg.debug) debug << d;
}
