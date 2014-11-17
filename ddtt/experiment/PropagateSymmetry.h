#pragma once

#include "ShapeGraph.h"

class PropagateSymmetry
{
public:
	static void propagate(const QStringList & fixedNodes, Structure::ShapeGraph * graph);
};
