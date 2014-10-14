#pragma once

#include "ShapeGraph.h"
#include "RenderObjectExt.h"

class Deformer
{
public:
	Deformer(Structure::ShapeGraph * a, Structure::ShapeGraph * b, int num_solve_iterations);

	QVector<RenderObject::Base*> debug;
};
