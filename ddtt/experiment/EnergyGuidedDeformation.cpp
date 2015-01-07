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
		origShapeA = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(*shapeA));
		origShapeB = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(*shapeB));

		preprocess(origShapeA.data(), origShapeB.data());
	}

	// Start search from roots
	for (auto & searchTree : searchTrees)
	{
		auto & root = *searchTree.begin();

		// Make copies
		root.shapeA = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(*origShapeA));
		root.shapeB = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(*origShapeB));

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

			// Apply deformation given current assignment
			applyAssignment(&path, false);

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
		}
	}
}

void Energy::GuidedDeformation::applyAssignment(Energy::SearchNode * path, bool isSaveKeyframes)
{
	double prevEnergy = path->energy;

	// Go over and apply the suggested assignments:
	for (auto ap : path->assignments)
	{
		// Get assignment pair <source, target>
		auto la = ap.first, lb = ap.second;
		assert(la.size() && lb.size());

		// Apply any needed topological operations
		topologicalOpeartions(path->shapeA.data(), path->shapeB.data(), la, lb);

		// Assigned parts will be fixed
		for (auto partID : ap.first) path->current << partID;

		// Deform the assigned
		applyDeformation(path->shapeA.data(), path->shapeB.data(), la, lb, path->fixed + path->current, isSaveKeyframes);

		// Track established correspondence
		for (size_t i = 0; i < la.size(); i++) path->mapping[la[i]] = lb[i].split(",").front();
	}

	// Evaluate distortion of shape
	double curEnergy = EvaluateCorrespondence::evaluate(path);

	path->cost = curEnergy - prevEnergy;
	path->energy = path->energy + path->cost;
}

QVector<Energy::SearchNode> Energy::GuidedDeformation::suggestChildren(Energy::SearchNode & path)
{
	// Hard coded thresholding to limit search space
	double candidate_threshold = 0.5;
	double cost_threshold = 0.3;
	int k_top_candidates = 5;

	/// Suggest for next unassigned:
	QVector<Structure::Relation> candidatesA;

	// Start process from remaining unassigned parts if needed
	if (!path.unassigned.isEmpty())
	{
		// Suggest for all remaining unassigned
		for (auto & partID : path.unassigned){
			if (!path.shapeA->hasRelation(partID)) continue;
			auto r = path.shapeA->relationOf(partID);
			if (!candidatesA.contains(r)) candidatesA << r;
		}
	}

	// Suggest for each candidate
	QVector < QPair<Structure::Relation, Structure::Relation> > pairings;
	for (auto & relationA : candidatesA)
	{
		// Candidate target groups
		/// Thresholding [0]: possibly skip geometrically very different ones
		for (auto & relationB : path.shapeB->relations)
			pairings.push_back(qMakePair(relationA, relationB));
	}

	// Current shapes bounding boxes
	auto boxA = path.shapeA->bbox(), boxB = path.shapeB->bbox();

	// Output:
	QList< QPair<double, Energy::SearchNode> > evaluated_children;

	//#ifndef QT_DEBUG
	//#pragma omp parallel for
	//#endif
	for (int r = 0; r < pairings.size(); r++)
	{
		auto & pairing = pairings[r];

		auto & relationA = pairing.first;
		auto & relationB = pairing.second;

		// Relative position of centroid
		auto rboxA = path.shapeA->relationBBox(relationA);
		Vector3 rboxCenterA = (rboxA.center() - boxA.min()).array() / boxA.sizes().array();

		auto rboxB = path.shapeB->relationBBox(relationB);
		Vector3 rboxCenterB = (rboxB.center() - boxB.min()).array() / boxB.sizes().array();

		/// Thresholding [1]: Skip assignment if spatially too far
		auto dist = (rboxCenterA - rboxCenterB).norm();
		if (dist < candidate_threshold)
		{
			double curEnergy = 0;
			QStringList la = relationA.parts, lb = relationB.parts;

			// Case: many-to-many, find best matching between two sets
			if (la.size() != 1 && lb.size() != 1)
			{
				//Eigen::Vector4d centroid_coordinate(0.5, 0.5, 0, 0);
				auto intrestingCentroid = [&](Structure::ShapeGraph * shape, QString nodeID){
					Eigen::Vector4d c(0, 0, 0, 0);
					auto edges = shape->getEdges(nodeID);
					if (edges.empty()) return Eigen::Vector4d(0.5,0.5,0,0);
					for (auto e : edges) c += e->getCoord(nodeID).front();
					return Eigen::Vector4d(c / edges.size());
				};

				lb.clear();

				for (size_t i = 0; i < relationA.parts.size(); i++)
				{
					auto partID = relationA.parts[i];
					auto partA = path.shapeA->getNode(partID);

					auto centroid_coordinateA = intrestingCentroid(path.shapeA.data(), partID);
					Vector3 partCenterA = (partA->position(centroid_coordinateA) - rboxA.min()).array() / rboxA.sizes().array();

					partCenterA.array() *= rboxA.diagonal().normalized().array();

					QMap<double, QString> dists;
					for (auto tpartID : relationB.parts)
					{
						auto centroid_coordinateB = intrestingCentroid(path.shapeB.data(), tpartID);
						Vector3 partCenterB = (path.shapeB->getNode(tpartID)->position(centroid_coordinateB) - rboxB.min()).array() / rboxB.sizes().array();

						partCenterB.array() *= rboxA.diagonal().normalized().array();

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
			auto shapeA = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(*path.shapeA));
			auto shapeB = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(*path.shapeB));

			// Apply then evaluate cost of current suggestion
			auto copy_la = la, copy_lb = lb;
			topologicalOpeartions(shapeA.data(), shapeB.data(), copy_la, copy_lb);
			applyDeformation(shapeA.data(), shapeB.data(), copy_la, copy_lb, path.fixed + path.current + copy_la.toSet());

			// Evaluate:
			{
				SearchNode tempNode;
				tempNode.shapeA = shapeA;
				tempNode.shapeB = shapeB;
				tempNode.fixed = path.fixed + path.current + copy_la.toSet();
				tempNode.mapping = path.mapping;
				for (size_t i = 0; i < la.size(); i++) tempNode.mapping[copy_la[i]] = copy_lb[i].split(",").front();
				tempNode.unassigned = path.unassigned;
				for (auto p : la) tempNode.unassigned.remove(p);
				curEnergy = EvaluateCorrespondence::evaluate(&tempNode);
			}

			/// Thresholding [2]: Skip really bad assignments
			double cost = curEnergy - path.energy;
			if (cost < cost_threshold)
			{
				auto modifiedShapeA = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(*path.shapeA));
				auto modifiedShapeB = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(*path.shapeB));

				auto unassigned = path.unassigned;
				for (auto p : la) unassigned.remove(p);

				Assignments assignment;
				assignment << qMakePair(la, lb);
				assert(la.size() && lb.size());

				//#pragma omp critical
				evaluated_children.push_back(qMakePair(cost, SearchNode(modifiedShapeA, modifiedShapeB, 
					path.fixed + path.current, assignment, unassigned, path.mapping, cost, path.energy)));
			}
		}
	}

	/// Thresholding [3] : only accept the 'k' top suggestions
	qSort(evaluated_children);
	auto sorted_children = evaluated_children.toVector();
	sorted_children.resize(std::min(sorted_children.size(), k_top_candidates));

	QVector < Energy::SearchNode > accepted_children;
	for (auto & child : sorted_children) accepted_children.push_back(child.second);

	// Clean up of no longer needed data:
	path.shapeA.clear();
	path.shapeB.clear();

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

	// Many-many case:
	if (la.size() > 1 && lb.size() > 1)
	{
		// Unique parts
		QStringList la_set = la, lb_set = lb;
		
		if (shapeA->hasRelation(la.front())) la_set = shapeA->relationOf(la.front()).parts;
		if (shapeB->hasRelation(lb.front())) lb_set = shapeB->relationOf(lb.front()).parts;

		if (la.size() != la_set.size() || lb.size() != lb_set.size())
		{
			// Experimental: will allow sliding
			for (auto partID : la_set)
			{
				auto partA = shapeA->getNode(partID);
				QString situation = la_set.size() > lb_set.size() ? "isMerged" : "isSplit";
				partA->property[situation].setValue(true);
				partA->property["isManyMany"].setValue(true);
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
		snode_sheet->property["mesh"].setValue(snode->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >());

		// Replace curve with a squashed sheet
		shapeA->nodes.replace(shapeA->nodes.indexOf(snode), snode_sheet);

		// Prepare for evaluation
		EvaluateCorrespondence::sampleNode(shapeA, snode_sheet, shapeA->property["sampling_resolution"].toDouble());

		// Fix coordinates (not robust..)
		for (auto l : shapeA->getEdges(snode->id))
		{
			Array1D_Vector4d sheet_coords(1, aproxProjection(l->position(snode->id), snode_sheet));
			l->replaceForced(snode->id, snode_sheet, sheet_coords);

			// Prepare for evaluation
			l->property["orig_spokes"].setValue(EvaluateCorrespondence::spokesFromLink(shapeA, l));
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

			snode->property["isMerged"].setValue(true);
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

		// Clean up
		shapeA->removeNode(sheetid);
	}

	// MERGE Case: many curves - one curve
	if (la.size() > 1 && lb.size() == 1
		&& shapeA->getNode(la.front())->type() == Structure::CURVE
		&& shapeB->getNode(lb.front())->type() == Structure::CURVE)
	{
		auto tnodeID = lb.front();
		lb.clear();

		// Experimental: will allow sliding
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

		// Experimental: will allow sliding
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

		// Experimental: will allow sliding
		for (auto partID : la)
		{
			auto snode = shapeA->getNode(partID);
			StructureAnalysis::removeFromGroups(shapeA, snode);
			snode->property["isMerged"].setValue(true);
			lb << tnodeID;
		}
	}

	// pseudo-SPLIT Case: one sheet - many curves
	if (la.size() == 1 && lb.size() > 1
		&& shapeA->getNode(la.front())->type() == Structure::SHEET
		&& shapeB->getNode(lb.front())->type() == Structure::CURVE)
	{
		auto snode = shapeA->getNode(la.front());
		snode->property["isSplit"].setValue(true);

		QString newnode = Structure::ShapeGraph::convertCurvesToSheet(shapeB, lb, Structure::ShapeGraph::computeSideCoordinates());

		// Effectively perform a split by moving samples used for correspondence evaluation:
		{
			// Make a copy of the sheet
			auto snode_sheet = shapeA->getNode(la.front());
			auto snode_sheet_copy = snode_sheet->clone();
			snode_sheet_copy->id += "clone";
			shapeA->addNode(snode_sheet_copy);

			// Turn curves on target to a new sheet
			auto tnode_sheet = shapeB->getNode(newnode);

			// Deform the copy to target curves
			Structure::ShapeGraph::correspondTwoNodes(snode_sheet_copy->id, shapeA, tnode_sheet->id, shapeB);
			DeformToFit::registerAndDeformNodes(snode_sheet_copy, tnode_sheet);

			// Snap source sheet samples to closest target curve projection
			auto samples = snode_sheet->property["samples_coords"].value<Array2D_Vector4d>();
			assert(!samples.empty());

			auto tcurve = shapeB->getNode(lb.front());
			auto start_c = aproxProjection(tcurve->position(Eigen::Vector4d(0, 0, 0, 0)), (Structure::Sheet*)snode_sheet_copy);
			auto end_c = aproxProjection(tcurve->position(Eigen::Vector4d(1, 0, 0, 0)), (Structure::Sheet*)snode_sheet_copy);
			Eigen::Vector4d range = (end_c - start_c).cwiseAbs();
			int snapIDX = range[0] > range[1] ? 1 : 0;

			auto coordSnap = [&](Eigen::Vector4d coord, int dim, int samplesCount, int curvesCount){
				if (curvesCount >= samplesCount) return coord;
				double interval = 1.0 / curvesCount;
				double v = coord[dim] < 0.5 ? floor(coord[dim] / interval) : ceil(coord[dim] / interval);
				coord[dim] = v * interval;
				return coord;
			};

			int samplesCount = samples.size();
			int curvesCount = lb.size();

			for (auto & row : samples){
				for (auto & c : row){
					c = coordSnap(c, snapIDX, samplesCount, curvesCount);
				}
			}

			// Replace samples
			snode_sheet->property["samples_coords"].setValue(samples);
				
			// Clean up
			shapeA->removeNode(snode_sheet_copy->id);
		}

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
		auto copy_lb = lb;

		auto snode = shapeA->getNode(la.front());

		QVector<Structure::Node *> new_nodes;
		new_nodes << snode;

		// Match the first one to its closest counter-part
		{
			auto boxA = shapeA->bbox(), boxB = shapeB->bbox();
			Vector3 snode_center = (snode->position(Eigen::Vector4d(0.5, 0.5, 0, 0)) - boxA.min()).array() / boxA.sizes().array();
			QMap<double, QString> dists;
			for (auto partID : lb){
				Vector3 tnode_center = (shapeB->getNode(partID)->position(Eigen::Vector4d(0.5, 0.5, 0, 0)) - boxB.min()).array() / boxB.sizes().array();
				dists[(tnode_center - snode_center).norm()] = partID;
			}
			QString matched_tnode = dists.values().front();
			std::swap(lb[0], lb[lb.indexOf(matched_tnode)]);
			copy_lb.removeAll(matched_tnode);
		}

		// Make remaining copies
		for (auto partID : copy_lb)
		{
			// Clone first node
			auto new_snode = snode->clone();
			new_snode->id += QString("@%1").arg(new_nodes.size());

			// Clone underlying geometry
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
	const QStringList & la, const QStringList & lb, const QSetString & fixed, bool isSaveKeyframes)
{
	// Special case: many-to-null
	if (lb.contains(Structure::null_part)) return;

	// Save initial configuration
	//if (isSaveKeyframes) shapeA->pushKeyframeDebug(new RenderObject::Text(30, 30, "now deform", 15));
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

	//if (isSaveKeyframes) shapeA->pushKeyframeDebug(new RenderObject::Text(30, 30, "now symmetry", 15));
	PropagateSymmetry::propagate(fixed, shapeA); if (isSaveKeyframes) shapeA->saveKeyframe();

	// Propagate edit by applying structural constraints
	//if (isSaveKeyframes) shapeA->pushKeyframeDebug(new RenderObject::Text(30, 30, "now proximity", 15));
	PropagateProximity::propagate(fixed, shapeA); if (isSaveKeyframes) { shapeA->saveKeyframe(); }
	//PropagateSymmetry::propagate(fixed, shapeA); if (isSaveKeyframes) { shapeA->saveKeyframe(); }
	//PropagateProximity::propagate(fixed, shapeA); if (isSaveKeyframes) { shapeA->saveKeyframe(); }

	postDeformation(shapeA, fixed);
}

void Energy::GuidedDeformation::postDeformation(Structure::ShapeGraph * shape, const QSet<QString> & fixedSet)
{
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

	std::sort(children.begin(), children.end(), [&](Energy::SearchNode * a, Energy::SearchNode * b){ return a->cost < b->cost; });

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

void Energy::GuidedDeformation::applySearchPath(QVector<Energy::SearchNode*> path)
{
	path.front()->shapeA = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(*origShapeA));
	path.front()->shapeB = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(*origShapeB));

	for (size_t i = 0; i < path.size(); i++)
	{
		auto p = path[i];

		// Prepare using previous step
		if (i > 0)
		{
			p->shapeA = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(*path[i - 1]->shapeA.data()));
			p->shapeB = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(*path[i - 1]->shapeB.data()));
			p->energy = path[i - 1]->energy;
		}

		applyAssignment(p, true);
	}
}
