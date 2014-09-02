#pragma once
#include "ParticleMesh.h"

class Segmentation
{
public:
    Segmentation(ParticleMesh * fromMesh);
    ParticleMesh * s;

    void mergeConvex();
    void mergeSimilar();

    QVector<RenderObject::Base*> debug;
};

