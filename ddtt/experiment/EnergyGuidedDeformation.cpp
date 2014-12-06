#include <deque>
#include <stack>

#include "EnergyGuidedDeformation.h"

#include "StructureAnalysis.h"
#include "DeformToFit.h"
#include "PropagateSymmetry.h"
#include "PropagateProximity.h"
#include "EvaluateCorrespondence.h"

#include "hausdorff.h"

void Energy::GuidedDeformation::searchAll()
{
	if (search_paths.empty() || !search_paths.front().shapeA) return;

	origShapeA = new Structure::ShapeGraph(*search_paths.front().shapeA);
	origShapeB = new Structure::ShapeGraph(*search_paths.front().shapeB);

	// Analyze symmetry groups
	StructureAnalysis::analyzeGroups(origShapeA, false);
	StructureAnalysis::analyzeGroups(origShapeB, false);

	// Prepare for proximity propagation
	PropagateProximity::prepareForProximity(origShapeA);

	// Prepare for structure distortion evaluation
	EvaluateCorrespondence::prepare(origShapeA);
	EvaluateCorrespondence::prepare(origShapeB);

	// Compression:
	{
		for (auto n : origShapeA->nodes){
			auto i = (byte)idxPartMapA.size();
			idxPartMapA[i] = n->id;
			partIdxMapA[n->id] = i;
		}

		for (auto n : origShapeB->nodes){
			auto i = (byte)idxPartMapB.size();
			idxPartMapB[i] = n->id;
			partIdxMapB[n->id] = i;
		}
	}

	for (auto & root : search_paths)
	{
		// Make copies
		root.shapeA = new Structure::ShapeGraph(*origShapeA);
		root.shapeB = new Structure::ShapeGraph(*origShapeB);

		root.unassigned = root.unassignedList();

		explore(root);
	}
}

void Energy::GuidedDeformation::explore( Energy::SearchPath & root )
{
	if (root.assignments.empty() && root.unassigned.empty()) return;

	std::stack < Energy::SearchPath* > path_stack;
	path_stack.push(&root);
	bool isRoot = true;

	while (!path_stack.empty())
	{
		auto & path = *path_stack.top();
		path_stack.pop();

		// Go over and apply previously suggested assignments:
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
			applyDeformation(path.shapeA, path.shapeB, la, lb, path.fixed + path.current);

			// Track established correspondence
			for (size_t i = 0; i < la.size(); i++) path.mapping[la[i]] = lb[i].split(",").front();
		}

		// Evaluate distortion of shape
		path.cost = EvaluateCorrespondence::evaluate(path.shapeA);

		double candidate_threshold = 0.5;
		double cost_threshold = 0.5;

		// Perform full search if we can afford it
		int shape_relations_threshold = 7;
		auto max_num_relations = std::max(path.shapeA->relations.size(), path.shapeB->relations.size());
		if (max_num_relations < shape_relations_threshold) candidate_threshold = cost_threshold = 2.0;

		/// Suggest for next unassigned:
		QVector<Structure::Relation> candidatesA;

		// Collect set of next candidates to be assigned by proximity to currently assigned
		for (auto partID : path.current)
		{
			for (auto edge : path.shapeA->getEdges(partID)){
				auto other = edge->otherNode(partID);
				if (path.fixed.contains(other->id)) continue;

				auto r = path.shapeA->relationOf(other->id);
				if (!candidatesA.contains(r)) candidatesA << r;
			}
		}

		// Collect next candidate to be assigned by finding part with good equivalent on target
		/*if (!path.unassigned.isEmpty())
		{
			auto dists = EvaluateCorrespondence::hausdroffDistance(path.shapeA, path.shapeB);

			double min_dist = DBL_MAX;
			QString min_part = path.unassigned.front();

			for (auto partA : dists.keys()){
				if (!path.unassigned.contains(partA)) continue;

				for (auto partB : dists[partA].keys()){
					double dist = dists[partA][partB];
					if (dist < min_dist){
						min_dist = dist;
						min_part = partA;
					}
				}
			}

			candidatesA << path.shapeA->relationOf(min_part);
		}*/

		// Start process from remaining unassigned parts if needed
		if (candidatesA.empty() && !path.unassigned.isEmpty())
		{
			// Suggest for all unassigned
			//for (auto partID : path.unassigned){
			//	auto r = path.shapeA->relationOf(partID);
			//	if (!candidatesA.contains(r)) candidatesA << r;
			//}

			// Arbitrary choice
			//candidatesA << path.shapeA->relationOf(path.unassigned.front());

			// Pick next largest unassigned
			//QMap < double, QString > unassigned_size;
			//for (auto partID : path.unassigned) unassigned_size[path.shapeA->robustBBox(partID).diagonal().norm()] = partID;
			//candidatesA << path.shapeA->relationOf( unassigned_size.values().back() );

			// Start from an arbitrary node in fixed set
			QStringList nowFixed = path.fixed + path.current;
			for (auto partID : path.mapping.keys()){
				for (auto edge : path.shapeA->getEdges(partID)){
					auto other = edge->otherNode(partID);
					if (nowFixed.contains(other->id)) continue;

					auto r = path.shapeA->relationOf(other->id);
					if (!candidatesA.contains(r)) candidatesA << r;
					break;
				}
				if (!candidatesA.empty())
					break;
			}
		}

		// Current shapes bounding boxes
		auto boxA = path.shapeA->bbox(), boxB = path.shapeB->bbox();

		// Suggest for each candidate
		QVector < QPair<Structure::Relation, Structure::Relation> > pairings;
		for (auto relationA : candidatesA)
		{
			// Create a null
			Structure::Relation null_relation;
			null_relation.type = Structure::Relation::NULLRELATION;
			null_relation.parts.push_back(Structure::null_part);

			// Find candidate target groups
			for (auto relationB : path.shapeB->relations /*+ (QVector<Structure::Relation>() << null_relation)*/)
				pairings.push_back(qMakePair(relationA, relationB));
		}

		QVector<Energy::SearchPath> suggested_children(pairings.size());

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
				//double cost = EvaluateCorrespondence::evaluate(&shapeA);
				{
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

					suggested_children[r] = SearchPath(modifiedShapeA, modifiedShapeB, path.fixed + path.current, assignment, canadidate_unassigned, path.mapping);
				}
			}
			else
			{
				auto rboxB = path.shapeB->relationBBox(relationB);
				Vector3 rboxCenterB = (rboxB.center() - boxB.min()).array() / boxB.sizes().array();
				auto dist = (rboxCenterA - rboxCenterB).norm();

				// Skip assignment if spatially too far
				if (dist > candidate_threshold) 
					continue;

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

				double diff_cost = abs(cost - path.cost);

				// Consider ones we like
				if (diff_cost < cost_threshold)
				{
					auto modifiedShapeA = new Structure::ShapeGraph(*path.shapeA);
					auto modifiedShapeB = new Structure::ShapeGraph(*path.shapeB);

					QStringList canadidate_unassigned = path.unassigned;
					for (auto p : la) canadidate_unassigned.removeAll(p);

					Assignments assignment;
					assignment << qMakePair(la, lb);

					assert(la.size() && lb.size());

					suggested_children[r] = SearchPath(modifiedShapeA, modifiedShapeB, path.fixed + path.current, assignment, canadidate_unassigned, path.mapping);
				}
			}
		}

		// Collect valid suggestions
		for (auto & child : suggested_children) 
			if (child.shapeA) 
				path.children.push_back(child);

		// Explore each suggestion:
		for (auto & child : path.children)
			path_stack.push(&child);

		// Memory saving:
		if (isRoot) isRoot = false;
		else
		{ 
			delete path.shapeA; 
			delete path.shapeB; 
			path.shapeA = path.shapeB = NULL;
		}
	}
}

QVector<Energy::SearchPath*> Energy::GuidedDeformation::solutions()
{
	auto t = Energy::SearchPath::exploreAsTree(search_paths);

	QVector<Energy::SearchPath*> result;
	for (auto leaf = t.begin_leaf(); leaf != t.end_leaf(); leaf++)
		if((*leaf)->unassigned.isEmpty()) 
			result.push_back(*leaf);

	/* Case: All paths where prematurely terminated */
	if (result.empty())
		for (auto leaf = t.begin_leaf(); leaf != t.end_leaf(); leaf++)
			result.push_back(*leaf);

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

		// Compression:
		if (!partIdxMapA.contains(snode_sheet->id))
		{
			auto i = (byte)idxPartMapA.size();
			idxPartMapA[i] = snode_sheet->id;
			partIdxMapA[snode_sheet->id] = i;
		}

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

		// Compression:
		if (!partIdxMapA.contains(snode_sheet->id))
		{
			auto i = (byte)idxPartMapA.size();
			idxPartMapA[i] = snode_sheet->id;
			partIdxMapA[snode_sheet->id] = i;
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
			if(partID != la.front()) snode->property["isMerged"].setValue(true);
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

		// Compression:
		if (!partIdxMapB.contains(newnode))
		{
			auto i = (byte)idxPartMapB.size();
			idxPartMapB[i] = newnode;
			partIdxMapB[newnode] = i;
		}
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
		int count = 0;

		for (auto partID : lb)
		{
			count++;
			if (count == 1) continue; // skip first

			auto new_snode = snode->clone();
			new_snode->id += QString("@%1").arg(count);

			shapeA->addNode(new_snode);
			new_nodes << new_snode;

			// Compression:
			if (!partIdxMapB.contains(new_snode->id))
			{
				auto i = (byte)idxPartMapB.size();
				idxPartMapB[i] = new_snode->id;
				partIdxMapB[new_snode->id] = i;
			}
		}

		la.clear();

		for (auto n : new_nodes){
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

void Energy::GuidedDeformation::applySearchPath(const QVector<Energy::SearchPath*> & path)
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

void Energy::SearchPath::compress(const QMap<QString, PartIndex> & mapA, const QMap<QString, PartIndex> & mapB)
{
	assert(_assignments.empty());

	// Convert part IDs from strings to small numbers
	for (auto & a : assignments){
		QVector < SearchPath::PartIndex > mappedA, mappedB;
		for (auto p : a.first) mappedA << mapA[p];
		for (auto p : a.second) mappedB << mapB[p];
		_assignments.push_back(qMakePair(mappedA, mappedB));
	}

	for (auto & p : fixed) _fixed.push_back(mapA[p]);
	for (auto & p : current) _current.push_back(mapA[p]);
	for (auto & p : unassigned) _unassigned.push_back(mapA[p]);
	for (auto & k : mapping.keys()) _mapping[mapA[k]] = mapB[mapping[k]];

	// Clean up
	assignments.clear();
	fixed.clear();
	current.clear();
	unassigned.clear();
	mapping.clear();
}

void Energy::SearchPath::decompress(const QMap<PartIndex, QString> & mapA, const QMap<PartIndex, QString> & mapB)
{
	// Reconstruct part IDs from small numbers
	for (auto a : _assignments){
		QStringList mappedA, mappedB;
		for (auto p : a.first) mappedA << mapA[p];
		for (auto p : a.second) mappedB << mapB[p];
		assignments.push_back(qMakePair(mappedA, mappedB));
	}

	for (auto & p : _fixed) fixed.push_back(mapA[p]);
	for (auto & p : _current) current.push_back(mapA[p]);
	for (auto & p : _unassigned) unassigned.push_back(mapA[p]);
	for (auto & k : _mapping.keys()) mapping[mapA[k]] = mapB[_mapping[k]];

	// Clean up
	_assignments.clear();
	_fixed.clear();
	_current.clear();
	_unassigned.clear();
	_mapping.clear();
}
