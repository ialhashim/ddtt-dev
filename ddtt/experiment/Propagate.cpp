#include "Propagate.h"
#include "StructureGraph.h"
#include "GenericGraph.h"

void Propagate::propagateProximity(const QStringList &fixedNodes, Structure::Graph *graph)
{
	GenericGraphs::Graph<size_t, double> g;
	for (auto n : graph->nodes) g.AddVertex(n->property["index"].toULongLong());

	auto fixedSource = g.vertices.size();
	g.AddVertex(fixedSource);
	g.initAttributes();

	for (auto edge : graph->edges){
		size_t nid1 = edge->n1->property["index"].toULongLong();
		size_t nid2 = edge->n2->property["index"].toULongLong();
		g.AddEdge(nid1, nid2, 1, edge->property["uid"].toInt());

		g.vertices_prop[nid1]["id"] = edge->n1->id;
		g.vertices_prop[nid2]["id"] = edge->n2->id;

		if (fixedNodes.contains(edge->n1->id)) g.AddEdge(fixedSource, nid1, 0);
		if (fixedNodes.contains(edge->n2->id)) g.AddEdge(fixedSource, nid1, 0);
	}

	// Spanning tree from all fixed nodes
	auto propagationTree = g.BFS(fixedSource);
	propagationTree.size();

	// Setup constraints based on case (1,2,3,..) and weather or not a part is fixed

	// Fixed point stretching or translation
}

void Propagate::propagateSymmetry(const QStringList &fixedNodes, Structure::Graph *graph)
{

}
