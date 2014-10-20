#pragma once

#include "ShapeGraph.h"
#include "myglobals.h"

typedef QVector<QPair<QVector<QStringList>, QVector<QStringList> > > Paths;

class CorrespondenceGenerator
{
public:
    CorrespondenceGenerator(Structure::ShapeGraph * shapeA, Structure::ShapeGraph * shapeB);
    Structure::ShapeGraph *a, *b;

	Paths generate();
};
