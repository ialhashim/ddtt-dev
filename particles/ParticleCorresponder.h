#pragma once
#include "ParticleMesh.h"
#include "PartCorresponder.h"

class ParticleCorresponder
{
public:
    ParticleCorresponder( ParticleMesh * pmeshA, ParticleMesh * pmeshB );
	Particles partToPartCorrespondence( const QVector< std::pair< size_t,size_t> > & partToPartAssignments );
    ParticleMesh *sA, *sB;

    QVector<RenderObject::Base*> debug;
};
