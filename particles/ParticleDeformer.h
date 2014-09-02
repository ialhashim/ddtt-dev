#pragma once
#include "ParticleMesh.h"

class ParticleDeformer
{
public:
    ParticleDeformer( ParticleMesh * pmeshA, ParticleMesh * pmeshB );

    ParticleMesh *sA, *sB;

    QVector<RenderObject::Base*> debug;
};
