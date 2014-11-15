#include "EnergyGuidedDeformation.h"
#include <QStack>

#include "DeformToFit.h"
#include "Propagate.h"
#include "StructureAnalysis.h"
#include "disjointset.h"

namespace Energy{
	struct SearchPath{
		QStringList parts, targets;
		QVector < QStringList > remaining, remainingTarget;
		double cost;
	};
}

EnergyGuidedDeformation::EnergyGuidedDeformation(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB,
                                                 const QVector<QStringList> &a_landmarks, const QVector<QStringList> &b_landmarks,
                                                 bool debugging) : a(shapeA), b(shapeB)
{
	// Analyze groups
	StructureAnalysis::analyzeGroups(shapeA, false);

	QStack < Energy::SearchPath > searchNodes;

	// Combine nodes that belong to same structural relation group
	QVector<QStringList> combined_landmarks_a, combined_landmarks_b;
	DisjointSet disjoint(a_landmarks.size());
	for (size_t i = 0; i < a_landmarks.size(); i++){
		for (size_t j = 0; j < a_landmarks.size(); j++){
			if (shapeA->groupsOf(a_landmarks[i].front()) == shapeA->groupsOf(a_landmarks[j].front()))
				disjoint.Union(i, j);
		}
	}
	for (auto group : disjoint.Groups()){
		QStringList al, bl;
		auto matchInNumber = [&](QStringList & i, QStringList & j){
			if (i.size() == j.size()) return;
			assert(i.size() == 1 || j.size() == 1);
			if (i.size() == 1)	for (int x = 0; x < j.size() - i.size(); x++) i << i.front();
			else				for (int x = 0; x < i.size() - j.size(); x++) j << j.front();
		};
		for (auto idx : group){
			auto p = a_landmarks[idx], q = b_landmarks[idx];
			matchInNumber(p,q);
			al += p;
			bl += q;
		}
		combined_landmarks_a << al;
		combined_landmarks_b << bl;
	}

	// Prepare a search path where there are fixed nodes and uncorresponded ones
	for (size_t i = 0; i < combined_landmarks_a.size(); i++)
	{
		Energy::SearchPath searchNode;
		searchNode.parts = combined_landmarks_a[i];
		searchNode.targets = combined_landmarks_b[i];

		auto remaining = combined_landmarks_a;
		auto remainingTarget = combined_landmarks_b;
		remaining.removeAt(i);
		remainingTarget.removeAt(i);
		searchNode.remaining = remaining;
		searchNode.remainingTarget = remainingTarget;

		searchNodes << searchNode;
	}

	// Prepare shape
    PropagateProximity::prepareForProximity(shapeA);

	shapeA->animation.push_back(shapeA->getAllControlPoints()); // DEBUG

	while (!searchNodes.isEmpty())
	{
		auto & searchNode = searchNodes.pop();

		// Deform part to its target
		for (size_t i = 0; i < searchNode.parts.size(); i++)
		{
			auto partID = searchNode.parts[i];
			auto tpartID = searchNode.targets[i];

			Structure::ShapeGraph::correspondTwoNodes(partID, a, tpartID, b); // does this help? why?
			DeformToFit::registerAndDeformNodes(a->getNode(partID), b->getNode(tpartID));
		}

		// Propagate edit by applying structural constraints
		shapeA->saveKeyframe();
        PropagateProximity::propagateProximity(searchNode.parts, shapeA);
		shapeA->saveKeyframe();
		Propagate::propagateSymmetry(searchNode.parts, shapeA);
		shapeA->saveKeyframe();
	}
}
