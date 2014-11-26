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
		StructureAnalysis::analyzeGroups(root.shapeB, false);

		// Prepare for proximity propagation
		PropagateProximity::prepareForProximity(root.shapeA);

		// Prepare for structure distortion evaluation
		EvaluateCorrespondence::prepare(root.shapeA);

		root.unassigned = root.unassignedList();

		explore(root);
	}
}

void Energy::GuidedDeformation::explore( Energy::SearchPath & path )
{
	if (path.assignments.empty() && path.unassigned.empty()) return;

	// Go over and apply previously suggested assignments:
	for (auto ap : path.assignments)
	{
		// Get assignment pair <source, target>
		auto la = ap.first, lb = ap.second;

		assert(la.size() && lb.size());

		// Apply any needed topological operations
		GuidedDeformation::topologicalOpeartions(path.shapeA, path.shapeB, la, lb);

		// Assigned parts will be fixed
		for (auto partID : ap.first) path.current << partID;

		// Deform the assigned
		GuidedDeformation::applyDeformation(path.shapeA, path.shapeB, la, lb, path.fixed + path.current);

		// Track established correspondence
		for (size_t i = 0; i < la.size(); i++) path.mapping[la[i]] = lb[i].split(",").front();
	}

	// Evaluate distortion of shape
	path.cost = EvaluateCorrespondence::evaluate(path.shapeA);

	double candidate_threshold = 0.3;
	double cost_threshold = 0.3;

	// Suggest for current unassigned:
	{
		// Collect set of next candidates to be assigned
		QVector<Structure::Relation> candidatesA;
		for (auto partID : path.current)
		{
			for (auto edge : path.shapeA->getEdges(partID)){
				auto other = edge->otherNode(partID);
				if (path.fixed.contains(other->id)) continue;

				auto r = path.shapeA->relationOf(other->id);
				if (!candidatesA.contains(r)) candidatesA << r;
			}
		}

		// Start process from remaining unassigned parts when needed
		if (candidatesA.empty() && !path.unassigned.isEmpty())
			candidatesA << path.shapeA->relationOf(path.unassigned.front());

		// Suggest for each candidate
		for (auto relationA : candidatesA)
		{
			// Relative position of centroid
			auto boxA = path.shapeA->bbox(), boxB = path.shapeB->bbox();
			auto rboxA = path.shapeA->relationBBox(relationA);
			Vector3 rboxCenterA = (rboxA.center() - boxA.min()).array() / boxA.sizes().array();

			Structure::Relation null_relation;
			null_relation.type = Structure::Relation::NULLRELATION;
			null_relation.parts.push_back( Structure::null_part );

			// Find candidate target groups
			for (auto & relationB : path.shapeB->relations + (QVector<Structure::Relation>() << null_relation))
			{
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

					// Evaluate
					double cost = EvaluateCorrespondence::evaluate(&shapeA);

					// Consider
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

						SearchPath child(modifiedShapeA, modifiedShapeB, path.fixed + path.current, assignment, canadidate_unassigned, path.mapping, cost);

						path.children.push_back(child);
					}

					continue;
				}

				auto rboxB = path.shapeB->relationBBox(relationB);
				Vector3 rboxCenterB = (rboxB.center() - boxB.min()).array() / boxB.sizes().array();
				auto dist = (rboxCenterA - rboxCenterB).norm();

				// Suggest or skip assignment
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
						Vector3 partCenterA = (path.shapeA->getNode(partID)->position(centroid_coordinate) - rboxA.center()).array() / rboxA.sizes().array();

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

				// Evaluate
				auto copy_la = la, copy_lb = lb;
				GuidedDeformation::topologicalOpeartions(&shapeA, &shapeB, copy_la, copy_lb);
				GuidedDeformation::applyDeformation(&shapeA, &shapeB, copy_la, copy_lb, path.fixed + la);
				cost = EvaluateCorrespondence::evaluate(&shapeA);
				
				// Consider ones we like
				if (cost < cost_threshold)
				{
					auto modifiedShapeA = new Structure::ShapeGraph(*path.shapeA);
					auto modifiedShapeB = new Structure::ShapeGraph(*path.shapeB);

					QStringList canadidate_unassigned = path.unassigned;
					for (auto p : la) canadidate_unassigned.removeAll(p);

					Assignments assignment;
					assignment << qMakePair(la,lb);

					assert(la.size() && lb.size());

					SearchPath child(modifiedShapeA, modifiedShapeB, path.fixed + path.current, assignment, canadidate_unassigned, path.mapping);

					path.children.push_back(child);
				}
			}
		}
	}

	// Explore each suggestion:
	for (auto & child : path.children)
		explore(child);
}

QVector<Energy::SearchPath*> Energy::GuidedDeformation::solutions()
{
	auto t = Energy::SearchPath::exploreAsTree(search_paths);

	QVector<Energy::SearchPath*> result;
	for (auto leaf = t.begin_leaf(); leaf != t.end_leaf(); leaf++)
		if((*leaf)->unassigned.isEmpty()) 
			result.push_back(*leaf);

	return result;
}

void Energy::GuidedDeformation::topologicalOpeartions(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB,
	QStringList & la, QStringList & lb)
{
	// Special case: many-to-null
	if (lb.contains(Structure::null_part))
	{
		return;
	}

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
	// Special case: many-to-null
	if (lb.contains(Structure::null_part)) return;

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
