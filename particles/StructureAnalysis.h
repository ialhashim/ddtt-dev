#pragma once

#include "ParticleMesh.h"

class StructureAnalysis
{
public:
    StructureAnalysis(ParticleMesh * pmesh);

	ParticleMesh * s;

	QMap<size_t, SegmentGraph*> allSegs;
	int pcount;

	QVector<RenderObject::Base*> debug;
};
