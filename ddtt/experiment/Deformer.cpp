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
			int count_u = s->numUCtrlPnts();
			int count_v = s->numVCtrlPnts();


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

	// Pack points into a matrix
	ShapeOp::Matrix3X points(3, allPnts.size());
	for (size_t i = 0; i < allPnts.size(); i++) points.col(i) = allPnts[i];

	// Create solver
	ShapeOp::Solver s;
	s.setPoints(points);

	// Setup constarints
	for (auto edge : edgeRelations)
	{
		double edge_weight = 1.0;
		std::vector<int> id_vector;

		id_vector.push_back(edge.first);
		id_vector.push_back(edge.second);

		auto c = std::make_shared<ShapeOp::EdgeStrainConstraint>(id_vector, edge_weight, points);
		s.addConstraint(c);
	}

	// Apply handles (closeness constraints)


	s.initialize();
	s.solve( num_solve_iterations );

	// Apply new positions back to shape
	ShapeOp::Matrix3X new_points = s.getPoints();

	for (auto nid : cpntsMap.keys())
	{
		auto n = a->getNode(nid);
		QVector<int> cpnt_ids = cpntsMap[nid];
		
		std::vector<Vector3> newCtrlPoints;
		for (int i = 0; i < cpnt_ids.size(); i++) newCtrlPoints.push_back(new_points.col(cpnt_ids[i]));
		n->setControlPoints(newCtrlPoints);
	}
}
