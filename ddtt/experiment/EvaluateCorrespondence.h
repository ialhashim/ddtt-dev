#pragma once
#include "ShapeGraph.h"

struct EvaluateCorrespondence
{
	static void prepare( Structure::ShapeGraph * shape );
    static void evaluate( QStringList fixedNodes, Structure::ShapeGraph * shape );
};
