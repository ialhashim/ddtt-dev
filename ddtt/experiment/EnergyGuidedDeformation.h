#pragma once

#include "ShapeGraph.h"
#include "RenderObjectExt.h"

class EnergyGuidedDeformation
{
public:
    EnergyGuidedDeformation(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB,
                               const QVector<QStringList> &a_landmarks, const QVector<QStringList> &b_landmarks,
                               bool debugging);
	Structure::ShapeGraph *a, *b;

	// Utility:

    // Debug:
    bool debugging;
    QVector<RenderObject::Base*> debug;
};
