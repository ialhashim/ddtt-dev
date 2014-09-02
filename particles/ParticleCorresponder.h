#pragma once
#include "ParticleMesh.h"

class ParticleCorresponder
{
public:
    ParticleCorresponder( ParticleMesh * pmeshA, ParticleMesh * pmeshB );

    ParticleMesh *sA, *sB;

	void basicCorrespondence();

    QVector<RenderObject::Base*> debug;
};
