#include "EnergyGuidedDeformation.h"

#include "StructureAnalysis.h"
#include "DeformToFit.h"
#include "PropagateSymmetry.h"
#include "PropagateProximity.h"
#include "EvaluateCorrespondence.h"

void Energy::GuidedDeformation::searchAll()
{
	for (auto & root : search_paths)
	{
		// Make copies
		root.shapeA = new Structure::ShapeGraph(*root.shapeA);
		root.shapeB = new Structure::ShapeGraph(*root.shapeB);

		// Analyze symmetry groups
		StructureAnalysis::analyzeGroups(root.shapeA, false);

		// Prepare for proximity propagation
		PropagateProximity::prepareForProximity(root.shapeA);

		// Prepare for structure distortion evaluation
		EvaluateCorrespondence::prepare(root.shapeA);

		explore(root);
	}
}

void Energy::GuidedDeformation::explore( SearchPath & path )
{
	if (path.assignments.empty() && path.unassigned.empty()) return;

	// Go over previously suggested assignments:
	while (!path.assignments.isEmpty())
	{
		// Get assignment pair <source, target>
		auto ap = path.assignments.pop();
		auto la = ap.first, lb = ap.second;

		// Apply any needed topological operations
		GuidedDeformation::topologicalOpeartions(path.shapeA, path.shapeB, la, lb);

		// Freeze these assigned parts
		for (auto partID : ap.first) path.fixed << partID;

		// Deform the assigned
		GuidedDeformation::applyDeformation(path.shapeA, path.shapeB, la, lb, path.fixed);
	}

	// Evaluate distortion of shape
	path.cost = EvaluateCorrespondence::evaluate(path.shapeA);

	// Suggest for current unassigned:
	for (auto partID : path.unassigned)
	{

	}

	// Explore each suggestion:
	for (auto & child : path.children)
		explore(child);
}

void Energy::GuidedDeformation::topologicalOpeartions(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB,
	QStringList & la, QStringList & lb)
{
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

		// Prepare for evaluation
		EvaluateCorrespondence::sampleNode(snode_sheet, shapeA->property["sampling_resolution"].toDouble());

		// Fix coordinates (not robust..)
		for (auto l : shapeA->getEdges(snode->id)) {
			Eigen::Vector4d bestUV(0.5, 0.5, 0, 0), minRange(0, 0, 0, 0), maxRange(1, 1, 0, 0);
			double currentDist = snode->area(), threshold = currentDist * 0.01;
			auto coord = snode_sheet->surface.timeAt(l->position(snode->id), bestUV, minRange, maxRange, currentDist, threshold);
			Array1D_Vector4d sheet_coords(1, coord);
			l->replaceForced(snode->id, snode_sheet, sheet_coords);

			// Prepare for evaluation
			l->property["orig_spokes"].setValue(EvaluateCorrespondence::spokesFromLink(l));
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
		Structure::ShapeGraph::correspondTwoNodes(snode_sheet->id, shapeA, tnode_sheet->id, shapeB);
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
}

void Energy::GuidedDeformation::applyDeformation(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB, 
	const QStringList & la, const QStringList & lb, const QStringList & fixed)
{
	// Save initial configuration
	shapeA->saveKeyframe();

	// Deform part to its target
	for (size_t i = 0; i < la.size(); i++)
	{
		auto partID = la[i];
		auto tpartID = lb[i];

		// does this help? why?
		if (shapeA->getNode(partID)->type() == Structure::SHEET && shapeB->getNode(tpartID)->type() == Structure::SHEET)
			Structure::ShapeGraph::correspondTwoNodes(partID, shapeA, tpartID, shapeB);

		DeformToFit::registerAndDeformNodes(shapeA->getNode(partID), shapeB->getNode(tpartID));
		shapeA->saveKeyframe();
		PropagateSymmetry::propagate(fixed, shapeA);
		shapeA->saveKeyframe();
	}

	// Propagate edit by applying structural constraints
	PropagateProximity::propagate(fixed, shapeA);
	shapeA->saveKeyframe();
	PropagateSymmetry::propagate(fixed, shapeA);
	shapeA->saveKeyframe();
	PropagateProximity::propagate(fixed, shapeA);
	shapeA->saveKeyframe();
}
