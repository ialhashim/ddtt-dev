#include "EnergyGuidedDeformation.h"
#include "DeformToFit.h"
#include <QStack>

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

	while (!searchNodes.isEmpty())
	{
		auto & searchNode = searchNodes.pop();

		// Deform part to its target
		for (auto partID : searchNode.parts)
		{
			for (auto tpartID : searchNode.targets)
			{
                DeformToFit::registerAndDeformNodes(a->getNode(partID), b->getNode(tpartID));
			}
		}
	}
}
