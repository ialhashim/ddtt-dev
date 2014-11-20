#pragma once
#include "ShapeGraph.h"

struct EvaluateCorrespondence
{
	static void prepare( Structure::ShapeGraph * shape );
	static double evaluate(Structure::ShapeGraph * shape);

	// Utility:
	static Array1D_Vector3 spokesFromLink(Structure::Link * link);
};
