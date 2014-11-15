#include "Propagate.h"
#include "ShapeGraph.h"
#include "StructureAnalysis.h"

void Propagate::propagateSymmetry(const QStringList &fixedNodes, Structure::ShapeGraph *graph)
{
	for (auto & relation : graph->relations)
	{
		if (relation.type == Structure::Relation::REFLECTIONAL)
		{
			auto partA = graph->getNode(relation.parts.front()), partB = graph->getNode(relation.parts.back());
			if (!fixedNodes.contains(partA->id)) std::swap(partA, partB);
			auto plane = StructureAnalysis::getReflectionalPlane(partA, partB);

			// Apply reflection to control points of other part
			Array1D_Vector3 cptsA = partA->controlPoints();
			Array1D_Vector3 cptsB;

			for (auto p : cptsA) cptsB.push_back( StructureAnalysis::pointReflection(p, plane.first, plane.second) );

			partB->setControlPoints(cptsB);
		}
	}
}
