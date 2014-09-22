#pragma once
#include "ParticleMesh.h"

typedef QVector< QPair<size_t, size_t> >  Pairings;
typedef QVector< QVector<size_t> > Assignments;

class CorrespondenceGenerator
{
public:
    CorrespondenceGenerator( ParticleMesh * pmeshA, ParticleMesh * pmeshB);

    ParticleMesh *sA, *sB;
	QVector<Pairings> computedAssignments;

    QVector<RenderObject::Base*> debug;

	Assignments generateGroupAssignments();
	QVector<Pairings> segmentAssignFromGroupAssign( Assignments groupAssignments );
};
