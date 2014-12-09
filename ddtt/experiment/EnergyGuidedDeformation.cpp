#include <deque>
#include <stack>

#include "EnergyGuidedDeformation.h"

#include "StructureAnalysis.h"
#include "DeformToFit.h"
#include "PropagateSymmetry.h"
#include "PropagateProximity.h"
#include "EvaluateCorrespondence.h"

#include "hausdorff.h"

#include "myglobals.h"

void Energy::GuidedDeformation::preprocess(Structure::ShapeGraph * shapeA, Structure::ShapeGraph * shapeB)
{
	// Analyze symmetry groups
	StructureAnalysis::analyzeGroups(shapeA, false);
	StructureAnalysis::analyzeGroups(shapeB, false);

	// Prepare for proximity propagation
	PropagateProximity::prepareForProximity(shapeA);
	PropagateProximity::prepareForProximity(shapeB);

	// Prepare for structure distortion evaluation
	EvaluateCorrespondence::prepare(shapeA);
	EvaluateCorrespondence::prepare(shapeB);
}

void Energy::GuidedDeformation::searchAll(Structure::ShapeGraph * shapeA, Structure::ShapeGraph * shapeB, QVector<Energy::SearchNode> & roots)
{
	if (roots.empty()) return;

	// Insert given roots of search trees
	for (auto & root : roots)
		searchTrees.push_back(SearchTree(root));

	// Pre-processing
	{
		origShapeA = new Structure::ShapeGraph(*shapeA);
		origShapeB = new Structure::ShapeGraph(*shapeB);

		preprocess(origShapeA, origShapeB);
	}

	// Start search from roots
	for (auto & searchTree : searchTrees)
	{
		auto & root = *searchTree.begin();

		// Make copies
		root.shapeA = new Structure::ShapeGraph(*origShapeA);
		root.shapeB = new Structure::ShapeGraph(*origShapeB);

		// Set as unassigned anything not set for this root
		root.unassigned = root.unassignedList();

		// Explore
		std::stack < Energy::SearchTree::iterator_base > path_stack;
		path_stack.push(searchTree.begin());

		while (!path_stack.empty())
		{
			auto pathItr = path_stack.top();
			path_stack.pop();

			auto & path = *pathItr;

			// Apply and evaluate deformation given current assignment
			applyAssignment(path, false);

			// Collect valid suggestions
			auto suggested_children = suggestChildren(path);
			for (auto & child : suggested_children)
				searchTree.append_child(pathItr, child);
			path.num_children = suggested_children.size();

			// Explore each suggestion:
			SearchTree::sibling_iterator child = searchTree.begin(pathItr);
			while (child != searchTree.end(pathItr)) {
				path_stack.push(child);
				++child;
			}

			// Clean up:
			if (pathItr.node->parent)
			{
				delete path.shapeA;
				delete path.shapeB;
				path.shapeA = path.shapeB = NULL;
			}
		}
	}
}

void Energy::GuidedDeformation::applyAssignment(Energy::SearchNode & path, bool isSaveKeyframes)
{
	// Go over and apply the suggested assignments:
	for (auto ap : path.assignments)
	{
		// Get assignment pair <source, target>
		auto la = ap.first, lb = ap.second;
		assert(la.size() && lb.size());

		// Apply any needed topological operations
		topologicalOpeartions(path.shapeA, path.shapeB, la, lb);

		// Assigned parts will be fixed
		for (auto partID : ap.first) path.current << partID;

		// Deform the assigned
		applyDeformation(path.shapeA, path.shapeB, la, lb, path.fixed + path.current, isSaveKeyframes);

		// Track established correspondence
		for (size_t i = 0; i < la.size(); i++) path.mapping[la[i]] = lb[i].split(",").front();
	}

	// Evaluate distortion of shape
	path.cost = EvaluateCorrespondence::evaluate(path.shapeA);
}

QVector<Energy::SearchNode> Energy::GuidedDeformation::suggestChildren(const Energy::SearchNode & path)
{
	// Hard coded thresholding to limit search space
	double candidate_threshold = 0.5;
	double cost_threshold = 0.5;

	// Perform full search if we can afford it
	int shape_relations_threshold = 6;
	auto max_num_relations = std::max(path.shapeA->relations.size(), path.shapeB->relations.size());
	int expand_by = 1;

	/// Suggest for next unassigned:
	QVector<Structure::Relation> candidatesA;

	// Start process from remaining unassigned parts if needed
	if (!path.unassigned.isEmpty())
	{
		// Suggest for all unassigned
		for (auto partID : path.unassigned){
			auto r = path.shapeA->relationOf(partID);
			if (!candidatesA.contains(r)) candidatesA << r;
		}

		if (max_num_relations < shape_relations_threshold)
			candidate_threshold = cost_threshold = 2.0;
		else
			candidatesA.resize(std::min(expand_by, candidatesA.size()));
	}

	// Suggest for each candidate
	QVector < QPair<Structure::Relation, Structure::Relation> > pairings;
	for (auto relationA : candidatesA){
		// Find candidate target groups
		for (auto relationB : path.shapeB->relations)
			pairings.push_back(qMakePair(relationA, relationB));
	}

	// Current shapes bounding boxes
	auto boxA = path.shapeA->bbox(), boxB = path.shapeB->bbox();

	// Output:
	QVector<Energy::SearchNode> suggested_children(pairings.size());

#ifndef QT_DEBUG
#pragma omp parallel for
#endif
	for (int r = 0; r < pairings.size(); r++)
	{
		auto & pairing = pairings[r];

		auto & relationA = pairing.first;
		auto & relationB = pairing.second;

		// Relative position of centroid
		auto rboxA = path.shapeA->relationBBox(relationA);
		Vector3 rboxCenterA = (rboxA.center() - boxA.min()).array() / boxA.sizes().array();

		// Special case: many-to-null:
		if (relationB.type == Structure::Relation::NULLRELATION)
		{
			// Make copy and match up
			Structure::ShapeGraph shapeA(*path.shapeA);
			QStringList la, lb;
			for (auto p : relationA.parts) lb << Structure::null_part;

			// Collapse part geometry to single point
			for (auto p : relationA.parts){
				auto n = shapeA.getNode(p);
				auto cpts = n->controlPoints();
				auto centroid = n->position(Eigen::Vector4d(1, 0, 0, 0));
				for (auto & p : cpts) p = centroid;
				n->setControlPoints(cpts);

				n->property["isAssignedNull"].setValue(true);
			}

			// Evaluate and consider
			double cost = EvaluateCorrespondence::evaluate(&shapeA);

			auto modifiedShapeA = new Structure::ShapeGraph(*path.shapeA);
			auto modifiedShapeB = new Structure::ShapeGraph(*path.shapeB);

			la = relationA.parts;
			QStringList canadidate_unassigned = path.unassigned;
			for (auto p : la)
			{
				canadidate_unassigned.removeAll(p);
				modifiedShapeA->getNode(p)->property["isAssignedNull"].setValue(true);
			}

			Assignments assignment;
			assignment << qMakePair(la, lb);

			assert(la.size() && lb.size());

			suggested_children[r] = SearchNode(modifiedShapeA, modifiedShapeB, path.fixed + path.current, assignment, canadidate_unassigned, path.mapping, cost);
		}
		else
		{
			auto rboxB = path.shapeB->relationBBox(relationB);
			Vector3 rboxCenterB = (rboxB.center() - boxB.min()).array() / boxB.sizes().array();

			// Thresholding: Skip assignment if spatially too far
			auto dist = (rboxCenterA - rboxCenterB).norm();
			if (dist > candidate_threshold) continue;

			// Try and evaluate suggestion
			double cost = 0;
			QStringList la = relationA.parts, lb = relationB.parts;

			// Case: many-to-many, find best matching between two sets
			if (la.size() != 1 && lb.size() != 1)
			{
				Eigen::Vector4d centroid_coordinate(0.5, 0.5, 0, 0);

				lb.clear();

				for (size_t i = 0; i < relationA.parts.size(); i++)
				{
					auto partID = relationA.parts[i];
					auto partA = path.shapeA->getNode(partID);
					Vector3 partCenterA = (partA->position(centroid_coordinate) - rboxA.center()).array() / rboxA.sizes().array();

					QMap<double, QString> dists;
					for (auto tpartID : relationB.parts)
					{
						Vector3 partCenterB = (path.shapeB->getNode(tpartID)->position(centroid_coordinate) - rboxB.center()).array() / rboxB.sizes().array();
						dists[(partCenterA - partCenterB).norm()] = tpartID;
					}

					lb << dists.values().front();
				}
			}
			else
			{
				la = relationA.parts;
				lb = relationB.parts;
			}

			// Make copies
			Structure::ShapeGraph shapeA(*path.shapeA), shapeB(*path.shapeB);

			// Evaluate cost of applying deformation and propagation
			auto copy_la = la, copy_lb = lb;
			topologicalOpeartions(&shapeA, &shapeB, copy_la, copy_lb);
			applyDeformation(&shapeA, &shapeB, copy_la, copy_lb, path.fixed + path.current + copy_la);
			cost = EvaluateCorrespondence::evaluate(&shapeA);

			// Thresholding: Skip really bad assignments
			double diff_cost = abs(cost - path.cost);
			if (diff_cost < cost_threshold)
			{
				auto modifiedShapeA = new Structure::ShapeGraph(*path.shapeA);
				auto modifiedShapeB = new Structure::ShapeGraph(*path.shapeB);

				QStringList unassigned = path.unassigned;
				for (auto p : la) unassigned.removeAll(p);

				Assignments assignment;
				assignment << qMakePair(la, lb);
				assert(la.size() && lb.size());

				suggested_children[r] = SearchNode(modifiedShapeA, modifiedShapeB, path.fixed + path.current, assignment, unassigned, path.mapping, cost);
			}
		}
	}

	QVector<Energy::SearchNode> accepted_children;

	for (auto & child : suggested_children)
		if (child.shapeA)
			accepted_children.push_back(child);

	return accepted_children;
}

QVector<Energy::SearchNode*> Energy::GuidedDeformation::solutions()
{
	auto & t = searchTrees.front();

	QVector<Energy::SearchNode*> result;
	for (auto leaf = t.begin_leaf(); leaf != t.end_leaf(); leaf++)
		if (leaf->unassigned.isEmpty())
			result.push_back(&(*leaf));

	/* Case: All paths where prematurely terminated */
	if (result.empty())
		for (auto leaf = t.begin_leaf(); leaf != t.end_leaf(); leaf++)
			result.push_back(&(*leaf));

	return result;
}

void Energy::GuidedDeformation::topologicalOpeartions(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB,
	QStringList & la, QStringList & lb)
{
	// Utility:
	auto aproxProjection = [](const Vector3 p, Structure::Sheet * sheet){
		Eigen::Vector4d bestUV(0.5, 0.5, 0, 0), minRange(0, 0, 0, 0), maxRange(1, 1, 0, 0);
		double avgEdge = sheet->avgEdgeLength(), threshold = avgEdge * 0.5;
		return sheet->surface.timeAt(p, bestUV, minRange, maxRange, avgEdge, threshold);
	};

	// NULL Special case: any-to-null
	if (lb.contains(Structure::null_part))
	{
		return;
	}

	// Detect and handle many-many cases:
	if (la.size() > 1 && lb.size() > 1)
	{
		// Unique parts
		auto la_set = shapeA->relationOf(la.front()).parts;
		auto lb_set = shapeB->relationOf(lb.front()).parts;

		if (la.size() != la_set.size() || lb.size() != lb_set.size())
		{
			// Experimental: allow sliding
			for (auto partID : la_set)
			{
				auto partA = shapeA->getNode(partID);
				QString situation = la_set.size() > lb_set.size() ? "isMerged" : "isSplit";
				partA->property[situation].setValue(true);
			}
		}
	}

	// CONVERT Case: curve to sheet (unrolling)
	if (la.size() == lb.size() && la.size() == 1
		&& shapeA->getNode(la.front())->type() == Structure::CURVE
		&& shapeB->getNode(lb.front())->type() == Structure::SHEET)
	{
		auto snode = shapeA->getNode(la.front());
		Array2D_Vector3 surface_cpts(4, snode->controlPoints());
		auto snode_sheet = new Structure::Sheet(NURBS::NURBSRectangled::createSheetFromPoints(surface_cpts), snode->id);

		// Replace curve with a squashed sheet
		shapeA->nodes.replace(shapeA->nodes.indexOf(snode), snode_sheet);

		// Prepare for evaluation
		EvaluateCorrespondence::sampleNode(snode_sheet, shapeA->property["sampling_resolution"].toDouble());

		// Fix coordinates (not robust..)
		for (auto l : shapeA->getEdges(snode->id))
		{
			Array1D_Vector4d sheet_coords(1, aproxProjection(l->position(snode->id), snode_sheet));
			l->replaceForced(snode->id, snode_sheet, sheet_coords);

			// Prepare for evaluation
			l->property["orig_spokes"].setValue(EvaluateCorrespondence::spokesFromLink(l));
		}

		// Remove from all relations
		StructureAnalysis::removeFromGroups(shapeA, snode);

		delete snode;
	}

	// MERGE Case: many curves - one sheet
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

			auto c1 = aproxProjection(snode->position(start_c), snode_sheet);
			auto c2 = aproxProjection(snode->position(end_c), snode_sheet);

			coords << qMakePair(c1, c2);
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

	// MERGE Case: many curves - one curve
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

	// MERGE Case: many sheets - one sheet
	if (la.size() > 1 && lb.size() == 1
		&& shapeA->getNode(la.front())->type() == Structure::SHEET
		&& shapeB->getNode(lb.front())->type() == Structure::SHEET)
	{
		auto tnodeID = lb.front();
		lb.clear();

		for (auto partID : la)
		{
			auto snode = shapeA->getNode(partID);
			StructureAnalysis::removeFromGroups(shapeA, snode);
			if (partID != la.front()) snode->property["isMerged"].setValue(true);
			lb << tnodeID;
		}
	}

	// MERGE Case: many sheets - one curve
	if (la.size() > 1 && lb.size() == 1
		&& shapeA->getNode(la.front())->type() == Structure::SHEET
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

	// SPLIT Case: one sheet - many curves
	if (la.size() == 1 && lb.size() > 1
		&& shapeA->getNode(la.front())->type() == Structure::SHEET
		&& shapeB->getNode(lb.front())->type() == Structure::CURVE)
	{
		QString newnode = Structure::ShapeGraph::convertCurvesToSheet(shapeB, lb, Structure::ShapeGraph::computeSideCoordinates());

		lb.clear();
		lb << newnode;
	}

	// SPLIT Case: one curve - many curves
	bool one_many_curves = la.size() == 1 && lb.size() > 1
		&& shapeA->getNode(la.front())->type() == Structure::CURVE
		&& shapeB->getNode(lb.front())->type() == Structure::CURVE;

	// SPLIT Case: one sheet - many sheets
	bool one_many_sheets = la.size() == 1 && lb.size() > 1
		&& shapeA->getNode(la.front())->type() == Structure::SHEET
		&& shapeB->getNode(lb.front())->type() == Structure::SHEET;

	// Perform split:
	if (one_many_curves || one_many_sheets)
	{
		auto snode = shapeA->getNode(la.front());

		QVector<Structure::Node *> new_nodes;
		new_nodes << snode;

		// Match the one to its closest counter-part
		auto boxA = shapeA->bbox(), boxB = shapeB->bbox();
		Vector3 snode_center = (snode->position(Eigen::Vector4d(0.5,0.5,0,0)) - boxA.min()).array() / boxA.sizes().array();
		QMap<double, QString> dists;
		for (auto partID : lb){
			Vector3 tnode_center = (shapeB->getNode(partID)->position(Eigen::Vector4d(0.5, 0.5, 0, 0)) - boxB.min()).array() / boxB.sizes().array();
			dists[(tnode_center - snode_center).norm()] = partID;
		}
		QString matched_tnode = dists.values().front();
		std::swap(lb[0], lb[lb.indexOf(matched_tnode)]);
		auto copy_lb = lb;
		copy_lb.removeAll(matched_tnode);

		// Remaining copies
		for (auto partID : copy_lb)
		{
			// Make node copy
			auto new_snode = snode->clone();
			new_snode->id += QString("@%1").arg(new_nodes.size());

			// Copy underlying geometry
			{
				auto orig_mesh = new_snode->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >().data();
				QSharedPointer<SurfaceMeshModel> new_mesh_ptr(orig_mesh->clone());
				new_mesh_ptr->updateBoundingBox();
				new_snode->property["mesh"].setValue(new_mesh_ptr);
			}

			shapeA->addNode(new_snode);
			new_nodes << new_snode;
		}

		la.clear();

		for (auto n : new_nodes)
		{
			n->property["isSplit"].setValue(true);
			la << n->id;
		}
	}
}

void Energy::GuidedDeformation::applyDeformation(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB,
	const QStringList & la, const QStringList & lb, const QStringList & fixed, bool isSaveKeyframes)
{
	// Special case: many-to-null
	if (lb.contains(Structure::null_part)) return;

	// Save initial configuration
	if (isSaveKeyframes) shapeA->pushKeyframeDebug(new RenderObject::Text(30, 30, "now deform", 15));
	if (isSaveKeyframes) shapeA->saveKeyframe();

	// Deform part to its target
	for (size_t i = 0; i < la.size(); i++)
	{
		auto partID = la[i];
		auto tpartID = lb[i];

		// does this help? why?
		if (shapeA->getNode(partID)->type() == Structure::SHEET && shapeB->getNode(tpartID)->type() == Structure::SHEET)
			Structure::ShapeGraph::correspondTwoNodes(partID, shapeA, tpartID, shapeB);

		DeformToFit::registerAndDeformNodes(shapeA->getNode(partID), shapeB->getNode(tpartID)); if (isSaveKeyframes) shapeA->saveKeyframe();
	}

	if (isSaveKeyframes) shapeA->pushKeyframeDebug(new RenderObject::Text(30, 30, "now symmetry", 15));
	PropagateSymmetry::propagate(fixed, shapeA); if (isSaveKeyframes) shapeA->saveKeyframe();

	// Propagate edit by applying structural constraints
	if (isSaveKeyframes) shapeA->pushKeyframeDebug(new RenderObject::Text(30, 30, "now proximity", 15));
	PropagateProximity::propagate(fixed, shapeA); if (isSaveKeyframes) { shapeA->saveKeyframe(); }
	//PropagateSymmetry::propagate(fixed, shapeA); if (isSaveKeyframes) { shapeA->saveKeyframe(); }
	//PropagateProximity::propagate(fixed, shapeA); if (isSaveKeyframes) { shapeA->saveKeyframe(); }

	postDeformation(shapeA, fixed);
}

void Energy::GuidedDeformation::postDeformation(Structure::ShapeGraph * shape, const QStringList & fixed)
{
	auto fixedSet = fixed.toSet();

	for (auto & r : shape->relations)
	{
		if (r.parts.empty() || !r.parts.toSet().intersect(fixedSet).empty()) continue;

		StructureAnalysis::updateRelation(shape, r);
	}
}

QVector<Energy::SearchNode*> Energy::GuidedDeformation::childrenOf(Energy::SearchNode * path)
{
	QVector<Energy::SearchNode*> children;

	auto & searchTree = searchTrees.front();
	auto pathItr = searchTree.begin();
	for (; pathItr != searchTree.end(); pathItr++) if (&(*pathItr) == path) break;

	SearchTree::sibling_iterator child = searchTree.begin(pathItr);
	while (child != searchTree.end(pathItr)) {
		children << &(*child);
		++child;
	}

	return children;
}

QVector<Energy::SearchNode*> Energy::GuidedDeformation::getEntirePath(Energy::SearchNode * path)
{
	QVector<Energy::SearchNode*> entirePath;

	auto & t = searchTrees.front();
	auto itr = t.begin();
	for (; itr != t.end(); itr++) if (&(*itr) == path) break;

	// Find ancestors
	auto current = itr;
	while (current != 0){
		auto & p = *current;
		entirePath.push_front(&p);
		current = t.parent(current);
	}

	return entirePath;
}

void Energy::GuidedDeformation::applySearchPath(const QVector<Energy::SearchNode*> & path)
{
	auto shapeA = new Structure::ShapeGraph(*origShapeA);
	auto shapeB = new Structure::ShapeGraph(*origShapeB);

	for (auto & p : path)
	{
		for (auto ap : p->assignments)
		{
			// Get assignment pair <source, target>
			auto la = ap.first, lb = ap.second;
			topologicalOpeartions(shapeA, shapeB, la, lb);
			applyDeformation(shapeA, shapeB, la, lb, p->fixed + p->current, true);

			p->shapeA = shapeA;
			p->shapeB = shapeB;
		}
	}
}
