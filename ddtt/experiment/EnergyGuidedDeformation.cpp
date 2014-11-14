#include "EnergyGuidedDeformation.h"
#include <QStack>

#include "DeformToFit.h"
#include "Propagate.h"

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
	QStack < Energy::SearchPath > searchNodes;

	for (size_t i = 0; i < a_landmarks.size(); i++)
	{
		Energy::SearchPath searchNode;
		searchNode.parts = a_landmarks[i];
		searchNode.targets = b_landmarks[i];

		auto remaining = a_landmarks;
		auto remainingTarget = b_landmarks;
		remaining.removeAt(i);
		remainingTarget.removeAt(i);
		searchNode.remaining = remaining;
		searchNode.remainingTarget = remainingTarget;

		searchNodes << searchNode;
	}

	// Prepare graph
    PropagateProximity::prepareForProximity(shapeA);

	while (!searchNodes.isEmpty())
	{
		auto & searchNode = searchNodes.pop();

		// Deform part to its target
		for (auto partID : searchNode.parts){
			for (auto tpartID : searchNode.targets){
                DeformToFit::registerAndDeformNodes(a->getNode(partID), b->getNode(tpartID));
			}
		}

		// Propagate edit by applying structural constraints
        PropagateProximity::propagateProximity(searchNode.parts, shapeA);

		Propagate::propagateSymmetry(searchNode.parts, shapeA);
	}
}
