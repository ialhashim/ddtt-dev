#pragma once
#include "ParticleMesh.h"

class ParticleCorresponder
{
public:
    ParticleCorresponder( ParticleMesh * pmeshA, ParticleMesh * pmeshB );

    ParticleMesh *sA, *sB;

	void descriptorCorrespondence();
	void basicCorrespondence();

    QVector<RenderObject::Base*> debug;
};
