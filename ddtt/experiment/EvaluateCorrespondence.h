#pragma once
#include "ShapeGraph.h"

struct EvaluateCorrespondence
{
	static void prepare( Structure::ShapeGraph * shape );
    static void evaluate( QStringList fixedNodes, Structure::ShapeGraph * shape );

	// Utility:
	static Array1D_Vector3 spokesFromLink(Structure::Link * link);
	static Array1D_Vector3 spokesFromShape(Structure::ShapeGraph * shape);
};
