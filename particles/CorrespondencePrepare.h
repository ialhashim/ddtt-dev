#pragma once
#include "ParticleMesh.h"

class CorrespondencePrepare
{
public:
    CorrespondencePrepare( ParticleMesh * pmeshA, ParticleMesh * pmeshB );

	ParticleMesh *sA, *sB;

	QVector<RenderObject::Base*> debug;

	std::vector< std::vector<size_t> > computeGroups( ParticleMesh * input );

	Particles segmentToSegmentCorrespondence( const QVector< QPair< size_t,size_t> > & segToSegAssignments );
};
