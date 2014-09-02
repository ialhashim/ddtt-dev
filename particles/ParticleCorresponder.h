#pragma once
#include "ParticleMesh.h"

class ParticleCorresponder
{
public:
    ParticleCorresponder( ParticleMesh * pmeshA, ParticleMesh * pmeshB );

    ParticleMesh *sA, *sB;

    QVector<RenderObject::Base*> debug;
};
