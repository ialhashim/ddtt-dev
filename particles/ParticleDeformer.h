#pragma once
#include "ParticleMesh.h"

typedef QVector< QPair<size_t, size_t> >  Pairings;
typedef QVector< QVector<size_t> > Assignments;

class ParticleDeformer
{
public:
    ParticleDeformer( ParticleMesh * pmeshA, ParticleMesh * pmeshB );

    ParticleMesh *sA, *sB;

    QVector<RenderObject::Base*> debug;

	Assignments generateGroupAssignments();
	QVector<Pairings> segmentAssignFromGroupAssign( Assignments groupAssignments );
};
