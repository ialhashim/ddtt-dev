#pragma once
#include "ShapeGraph.h"

struct EvaluateCorrespondence
{
	static void prepare( Structure::ShapeGraph * shape );
	static double evaluate(Structure::ShapeGraph * shape);

	// Utility:
	static Array1D_Vector3 spokesFromLink(Structure::Link * link);
	static void sampleNode(Structure::Node * n, double resolution);
};

Q_DECLARE_METATYPE(Vector3);
Q_DECLARE_METATYPE(Array1D_Vector3);
Q_DECLARE_METATYPE(Array2D_Vector3);
Q_DECLARE_METATYPE(Array2D_Vector4d);
