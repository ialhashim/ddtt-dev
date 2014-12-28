#pragma once
#include "ShapeGraph.h"
#include "NanoKdTree.h"
#include "EnergyGuidedDeformation.h"

struct EvaluateCorrespondence
{
	static void prepare( Structure::ShapeGraph * shape );
	static double evaluate(Energy::SearchNode * searchNode, bool isRecordDetails = false);

	// Utility:
	static Array1D_Vector3 spokesFromLink( Structure::Link * link );
	static Array2D_Vector4d sampleNode(Structure::Node * n, double resolution);
	static QMap<QString, NanoKdTree*> kdTreesNodes( Structure::ShapeGraph * shape );
	static QMap<QString, QMap<QString, double> > hausdroffDistance( Structure::ShapeGraph * shapeA, Structure::ShapeGraph * shapeB );
	static double RMSD(Structure::ShapeGraph * shapeA, Structure::ShapeGraph * shapeB);
};

Q_DECLARE_METATYPE(Vector3);
Q_DECLARE_METATYPE(Array1D_Vector3);
Q_DECLARE_METATYPE(Array2D_Vector3);
Q_DECLARE_METATYPE(Array2D_Vector4d);
