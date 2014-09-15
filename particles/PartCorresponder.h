#pragma once
#include "ParticleMesh.h"

class PartCorresponder
{
public:
    PartCorresponder( ParticleMesh * pmeshA, SegmentGraph segA,
                      ParticleMesh * pmeshB, SegmentGraph segB );

    ParticleMesh *sA, *sB;
    SegmentGraph segA, segB;

    QVector<RenderObject::Base*> debug;
};
