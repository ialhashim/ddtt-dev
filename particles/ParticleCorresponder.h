#pragma once
#include "ParticleMesh.h"
#include "PartCorresponder.h"

class ParticleCorresponder
{
public:
    ParticleCorresponder( ParticleMesh * pmeshA, ParticleMesh * pmeshB );

	ParticleMesh *sA, *sB;

	QVector<RenderObject::Base*> debug;

	std::vector< std::vector<size_t> > computeGroups( ParticleMesh * input );

	Particles segmentToSegmentCorrespondence( const QVector< std::pair< size_t,size_t> > & segToSegAssignments );
};
