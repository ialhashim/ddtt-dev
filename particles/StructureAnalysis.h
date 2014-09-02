#pragma once

#include "ParticleMesh.h"

class StructureAnalysis
{
public:
    StructureAnalysis( ParticleMesh * pmesh );
    ParticleMesh * s;

	QVector<RenderObject::Base*> debug;
};
