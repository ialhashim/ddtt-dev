#pragma once
#include "ParticleMesh.h"

class CorrespondencePrepare
{
public:
    CorrespondencePrepare( std::vector<ParticleMesh*> meshes );

	QVector<RenderObject::Base*> debug;

	std::vector< std::vector<size_t> > computeGroups( ParticleMesh * input );
};
