#include "EnergyGuidedDeformation.h"
#include <QStack>

#include "StructureAnalysis.h"
#include "DeformToFit.h"
#include "PropagateSymmetry.h"
#include "PropagateProximity.h"
#include "EvaluateCorrespondence.h"
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
                                                 bool debugging)
{
	shapeA = new Structure::ShapeGraph(*shapeA);
	shapeB = new Structure::ShapeGraph(*shapeB);
	this->a = shapeA;
	this->b = shapeB;

	// Analyze groups
	StructureAnalysis::analyzeGroups(shapeA, false);

	// Combine nodes that belong to same structural relation group
	QVector<QStringList> combined_landmarks_a, combined_landmarks_b;
	
	for (size_t i = 0; i < a_landmarks.size(); i++)
	{
		auto la = a_landmarks[i];
		auto lb = b_landmarks[i];

		// Case: curve to sheet (unrolling)
		if (la.size() == lb.size() && la.size() == 1 
			&& shapeA->getNode(la.front())->type() == Structure::CURVE 
			&& shapeB->getNode(lb.front())->type() == Structure::SHEET)
		{
			auto snode = shapeA->getNode(la.front());
			Array2D_Vector3 surface_cpts(4, snode->controlPoints());
			auto snode_sheet = new Structure::Sheet(NURBS::NURBSRectangled::createSheetFromPoints(surface_cpts), snode->id);

			// Replace curve with a squashed surface
			shapeA->nodes.replace(shapeA->nodes.indexOf(snode), snode_sheet);

			// Fix coordinates (not robust..)
			for (auto l : shapeA->getEdges(snode->id)) {
				Eigen::Vector4d bestUV(0.5, 0.5, 0, 0), minRange(0, 0, 0, 0), maxRange(1, 1, 0, 0);
				double currentDist = snode->area(), threshold = currentDist * 0.01;
				auto coord = snode_sheet->surface.timeAt(l->position(snode->id), bestUV, minRange, maxRange, currentDist, threshold);
				Array1D_Vector4d sheet_coords(1, coord);
				l->replaceForced(snode->id, snode_sheet, sheet_coords);
			}

			// Remove from all relations
			StructureAnalysis::removeFromGroups(shapeA, snode);

			delete snode;
		}

		// Case: one sheet - many curves
		if (la.size() == 1 && lb.size() > 1
			&& shapeA->getNode(la.front())->type() == Structure::SHEET
			&& shapeB->getNode(lb.front())->type() == Structure::CURVE)
		{
			QString newnode = Structure::ShapeGraph::convertCurvesToSheet(shapeB, lb, Structure::ShapeGraph::computeSideCoordinates());
			shapeB->getNode(newnode);

			lb.clear();
			lb << newnode;
		}

		// Case: many curves - one sheet
		if (la.size() > 1 && lb.size() == 1
			&& shapeA->getNode(la.front())->type() == Structure::CURVE
			&& shapeB->getNode(lb.front())->type() == Structure::SHEET)
		{
			auto tnode_sheet = (Structure::Sheet*)shapeB->getNode(lb.front());
			lb.clear();

			QString sheetid = Structure::ShapeGraph::convertCurvesToSheet(shapeA, la, Structure::ShapeGraph::computeSideCoordinates());
			Structure::Sheet * snode_sheet = (Structure::Sheet *)shapeA->getNode(sheetid);

			QVector<QPair<Eigen::Vector4d, Eigen::Vector4d>> coords;
			for (size_t i = 0; i < la.size(); i++)
			{
				auto snode = shapeA->getNode(la[i]);
				Eigen::Vector4d start_c(0, 0, 0, 0), end_c(1, 0, 0, 0);
				coords << qMakePair(snode_sheet->approxCoordinates(snode->position(start_c)), 
									snode_sheet->approxCoordinates(snode->position(end_c)));
			}

			// Should give a nice enough alignment
			Structure::ShapeGraph::correspondTwoNodes(snode_sheet->id, a, tnode_sheet->id, b);
			DeformToFit::registerAndDeformNodes(snode_sheet, tnode_sheet);

			// Create equivalent curves on target sheet
			for (size_t i = 0; i < la.size(); i++)
			{
				auto coord = coords[i];

				Vector3 start_point = snode_sheet->position(coord.first);
				Vector3 end_point = snode_sheet->position(coord.second);
				Vector3 direction = (end_point - start_point).normalized();

				auto tcurve = tnode_sheet->convertToNURBSCurve(start_point, direction);
				auto tcurve_id = tnode_sheet->id + "," + la[i];

				lb << shapeB->addNode(new Structure::Curve(tcurve, tcurve_id))->id;
			}

			shapeA->removeNode(sheetid);
		}

		// Case: many curves - one curve
		if (la.size() > 1 && lb.size() == 1
			&& shapeA->getNode(la.front())->type() == Structure::CURVE
			&& shapeB->getNode(lb.front())->type() == Structure::CURVE)
		{
			auto tnodeID = lb.front();
			lb.clear();

			for (auto partID : la)
			{
				auto snode = shapeA->getNode(partID);
				StructureAnalysis::removeFromGroups(shapeA, snode);

				snode->property["isMerged"].setValue(true);

				lb << tnodeID;
			}
		}

		combined_landmarks_a << la;
		combined_landmarks_b << lb;
	}

	// Prepare a search path where there are fixed nodes and uncorresponded ones
	QStack < Energy::SearchPath > searchNodes;
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

	EvaluateCorrespondence::prepare(shapeA);

	// Save initial configuration
	shapeA->saveKeyframe();

	while (!searchNodes.isEmpty())
	{
		auto & searchNode = searchNodes.pop();

		// Deform part to its target
		for (size_t i = 0; i < searchNode.parts.size(); i++)
		{
			auto partID = searchNode.parts[i];
			auto tpartID = searchNode.targets[i];

			// does this help? why?
			if (a->getNode(partID)->type() == Structure::SHEET && a->getNode(partID)->type() == b->getNode(tpartID)->type())
				Structure::ShapeGraph::correspondTwoNodes(partID, a, tpartID, b);

			DeformToFit::registerAndDeformNodes(a->getNode(partID), b->getNode(tpartID));
			shapeA->saveKeyframe();
			PropagateSymmetry::propagate(searchNode.parts, shapeA);
			shapeA->saveKeyframe();
		}

		// Propagate edit by applying structural constraints
		PropagateProximity::propagate(searchNode.parts, shapeA);
		shapeA->saveKeyframe();
		PropagateSymmetry::propagate(searchNode.parts, shapeA);
		shapeA->saveKeyframe();
		PropagateProximity::propagate(searchNode.parts, shapeA);
		shapeA->saveKeyframe();
	}

	shapeA->saveKeyframe();
}
