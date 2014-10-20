#pragma once

#include "ShapeGraph.h"
#include "RenderObjectExt.h"

#include "Constraint.h"
#include "Solver.h"

class Deformer
{
public:
	Deformer(Structure::ShapeGraph * a, Structure::ShapeGraph * b, int num_solve_iterations);
	
	static void shapeNodesToPoints(Structure::ShapeGraph * a,
		QMap< QString, QVector<int> > & cpntsMap,
		QVector< Vector3 > & allPnts);

	static void shapeNodesToConstraints(Structure::ShapeGraph * a,
		const QMap< QString, QVector<int> > & cpntsMap,
		const QVector< Vector3 > & allPnts,
		QVector < QPair<int, int> > & pointRelations);

	static void shapeEdgesToConstraints(Structure::ShapeGraph * a,
		const QMap< QString, QVector<int> > & cpntsMap,
		const QVector< Vector3 > & allPnts,
		QVector < QPair<int, int> > & pointRelations);

	QVector<RenderObject::Base*> debug;
};
