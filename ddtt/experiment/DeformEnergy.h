#pragma once

#include "ShapeGraph.h"
#include "RenderObjectExt.h"

class DeformEnergy
{
public:
    DeformEnergy(Structure::ShapeGraph * shapeA, Structure::ShapeGraph * shapeB,
                 const QVector<QStringList> &a_landmarks,
                 const QVector<QStringList> &b_landmarks,
                 bool debugging = true);

	Structure::ShapeGraph *a, *b;

	double deform(Structure::Node * nodeA, Structure::Node * nodeB, bool isTwistTerm);

    static QString convertCurvesToSheet(Structure::Graph * graph, QStringList & nodeIDs);

	static Array2D_Vector4d DeformEnergy::computeSideCoordinates(int resolution);
	static Array2D_Vector4d sideCoordinates;

	double error;

	bool debugging;
	QVector<RenderObject::Base*> debug;
};
