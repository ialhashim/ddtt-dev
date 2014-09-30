#pragma once
#include "ParticleMesh.h"

class CorrespondencePrepare
{
public:
    CorrespondencePrepare( std::vector<ParticleMesh*> meshes );

	QVector<RenderObject::Base*> debug;

	static std::vector< std::vector<size_t> > computeGroups( ParticleMesh * input, const Segments & segments );
};
