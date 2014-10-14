#include "Deformer.h"
#include "myglobals.h"

#include "Constraint.h"
#include "Solver.h"

Deformer::Deformer(Structure::ShapeGraph * a, Structure::ShapeGraph * b, int num_solve_iterations)
{
	// Collect control points of shape
	QMap< QString, QVector<int> > cpntsMap;
	QVector<Vector3> allPnts;
	for (auto n : a->nodes)
	{
		auto cpnts = n->controlPoints();
		
		for (int i = 0; i < (int)cpnts.size(); i++){
			cpntsMap[n->id].push_back( allPnts.size() );
			allPnts << cpnts[i];
		}
	}

	// Collect relationships between points
	QVector < QPair<int,int> > edgeRelations;
	for (auto n : a->nodes)
	{
		if (n->type() == Structure::SHEET)
		{
			auto s = (Structure::Sheet*) n;
			auto surface = s->surface;

			for (size_t u = 0; u < surface.mNumUCtrlPoints; u++){
				for (size_t v = 0; v < surface.mNumVCtrlPoints; v++){
					int idx0 = cpntsMap[n->id][v + (surface.mNumVCtrlPoints * u)];

					if (u + 1 < surface.mNumUCtrlPoints){
						int idx1 = cpntsMap[n->id][v + (surface.mNumVCtrlPoints * (u + 1))];
						edgeRelations.push_back(qMakePair(idx0, idx1));
					}

					if (v + 1 < surface.mNumVCtrlPoints){
						int idx2 = cpntsMap[n->id][(v + 1) + (surface.mNumVCtrlPoints * u)];
						edgeRelations.push_back(qMakePair(idx0, idx2));
					}
				}
			}
		}
		else
		{
			auto c = (Structure::Curve*) n;
			int count_u = c->numCtrlPnts();

			for (int ui = 0; ui + 1 < count_u; ui++)
			{
				int uj = ui + 1;

				int idx_ui = cpntsMap[n->id][ui];
				int idx_uj = cpntsMap[n->id][uj];

				edgeRelations.push_back( qMakePair(idx_ui, idx_uj) );
			}
		}
	}

	// Collect relationships between nodes
	for (auto edge : a->edges)
	{
		auto nid1 = edge->n1->id, nid2 = edge->n2->id;
		auto epos1 = edge->position(nid1), epos2 = edge->position(nid2);
		QMap <int, double> dists_ids1, dists_ids2;

		// Compute distance between encoded edge and control points
		for (size_t i = 0; i < cpntsMap[nid1].size(); i++){
			int idx = cpntsMap[nid1][i];
			dists_ids1[idx] = (allPnts[idx] - epos1).norm();
		}
		for (size_t i = 0; i < cpntsMap[nid2].size(); i++){
			int idx = cpntsMap[nid2][i];
			dists_ids2[idx] = (allPnts[idx] - epos2).norm();
		}
		
		// Find closest control points to encoded edge
		auto dists1 = dists_ids1.values(), dists2 = dists_ids2.values();
		auto smallest1 = std::min_element(dists1.begin(), dists1.end()) - dists1.begin();
		auto smallest2 = std::min_element(dists2.begin(), dists2.end()) - dists2.begin();
		auto smallest_id1 = dists_ids1.keys().at(smallest1);
		auto smallest_id2 = dists_ids2.keys().at(smallest2);

		// Add edge constraint across two nodes
		edgeRelations.push_back(qMakePair(smallest_id1, smallest_id2));
	}

	// Pack points into a matrix
	ShapeOp::Matrix3X points(3, allPnts.size());
	for (size_t i = 0; i < allPnts.size(); i++) points.col(i) = allPnts[i];

	// Create solver
	ShapeOp::Solver s;
	s.setPoints(points);

	// Setup edge constarints
	for (auto edge : edgeRelations)
	{
		double edge_weight = 1.0;
		std::vector<int> id_vector;

		id_vector.push_back(edge.first);
		id_vector.push_back(edge.second);

		auto c = std::make_shared<ShapeOp::EdgeStrainConstraint>(id_vector, edge_weight, points);
		s.addConstraint(c);
	}

	/// Apply handles (closeness constraints):
	/* By some arbitrary first control point */
	/*for (auto n : a->nodes){
		int firstElement = cpntsMap[n->id][0];
		std::vector<int> id_vector;
		id_vector.push_back( firstElement );
		auto c = std::make_shared<ShapeOp::ClosenessConstraint>(id_vector, close_weight, points);
		s.addConstraint(c);
	}*/
	/* By given sparse correspondence */
	for (auto & landmark : a->landmarks)
	{
		// Compute control point index closest to landmark
		QMap <int, double> dists_ids;
		for (size_t i = 0; i < cpntsMap[landmark.partid].size(); i++){
			int idx = cpntsMap[landmark.partid][i];
			dists_ids[idx] = (allPnts[idx] - landmark).norm();
		}
		auto dists = dists_ids.values();
		auto smallest_id = dists_ids.keys().at(std::min_element(dists.begin(), dists.end()) - dists.begin());

		// Add constraint
		double close_weight = 1.0;
		std::vector<int> id_vector;
		id_vector.push_back(smallest_id);

		auto c = std::make_shared<ShapeOp::ClosenessConstraint>(id_vector, close_weight, points);
		landmark.constraint_id = s.addConstraint(c);
	}

	s.initialize();

	/* Apply handle deformation */
	for (size_t l = 0; l < a->landmarks.size(); l++)
	{
		auto c = std::dynamic_pointer_cast <ShapeOp::ClosenessConstraint>(s.getConstraint(a->landmarks[l].constraint_id));
		c->setPosition( b->landmarks[l] );

		for(int i = 0; i < 3; i++) a->landmarks[l][i] = b->landmarks[l][i];
	}

	s.solve( num_solve_iterations );

	/// Apply new positions back to shape:
	ShapeOp::Matrix3X new_points = s.getPoints();

	for (auto nid : cpntsMap.keys())
	{
		auto n = a->getNode(nid);
		QVector<int> cpnt_ids = cpntsMap[nid];
		
		std::vector<Vector3> newCtrlPoints;
		for (int i = 0; i < cpnt_ids.size(); i++) newCtrlPoints.push_back(new_points.col(cpnt_ids[i]));
		n->setControlPoints(newCtrlPoints);

		// Re-bulid a surface
		if (n->type() == Structure::SHEET) ((Structure::Sheet*)n)->surface.quads.clear();
	}
}
