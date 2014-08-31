#pragma once

#include "ParticleMesh.h"

class StructureAnalysis
{
public:
	StructureAnalysis(ParticleMesh * pmesh);
	ParticleMesh * s;

	void mergeConvex();
	void mergeSimilar();

	QVector<RenderObject::Base*> debug;
};
