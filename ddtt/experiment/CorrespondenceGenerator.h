#pragma once

#include "ShapeGraph.h"
#include "myglobals.h"

typedef QPair < QVector<QStringList>, QVector<QStringList> > AssignmentPair;
typedef QVector< AssignmentPair > Paths;

class CorrespondenceGenerator
{
public:
    CorrespondenceGenerator(Structure::ShapeGraph * shapeA, Structure::ShapeGraph * shapeB);
    Structure::ShapeGraph *a, *b;

	Paths generate();
};
