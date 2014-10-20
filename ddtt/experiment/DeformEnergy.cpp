#include "DeformEnergy.h"
#include <Eigen/Geometry>
#include "myglobals.h"

Array2D_Vector4d DeformEnergy::sideCoordinates = DeformEnergy::computeSideCoordinates(20);

DeformEnergy::DeformEnergy(Structure::ShapeGraph * shapeA, Structure::ShapeGraph * shapeB,
                           const QVector<QStringList> & a_landmarks,
                           const QVector<QStringList> & b_landmarks,
                           bool debugging) : a(shapeA), b(shapeB), error(0), debugging(debugging)
{
	// Deform geometry
    for (size_t i = 0; i < a_landmarks.size(); i++)
	{
        auto la = a_landmarks[i];
        auto lb = b_landmarks[i];
		bool isOneA = la.size() == 1;
		bool isOneB = lb.size() == 1;
        bool isSheetA = a->getNode(la.front())->type() == Structure::SHEET;
        bool isSheetB = b->getNode(lb.front())->type() == Structure::SHEET;
		bool isOneToMany = (isOneA != isOneB);
		bool isOneIsSheet = (isSheetA || isSheetB);
		bool isOtherNotSheet = !isSheetA || !isSheetB;

		bool isNodeAdded = false;
		Structure::ShapeGraph * graph;
        QStringList * landmarks;

		if (isOneToMany && isOneIsSheet && isOtherNotSheet)
		{
			graph = (isOneA) ? shapeB : shapeA;
			landmarks = (isOneA) ? &lb : &la;

            QString newnode = DeformEnergy::convertCurvesToSheet(graph, *landmarks);

			landmarks->clear();
            landmarks->push_back(newnode);

			isNodeAdded = true;
		}

		for (size_t u = 0; u < la.size(); u++)
		{
			for (size_t v = 0; v < lb.size(); v++)
			{
                auto nodeA = a->getNode(la[u]);
                auto nodeB = b->getNode(lb[v]);

				double area = deform(nodeA, nodeB);

				error += area;
			}
		}

		if (isNodeAdded)
            graph->removeNode(landmarks->front());
	}
}

QString DeformEnergy::convertCurvesToSheet(Structure::Graph * graph, QStringList & nodeIDs)
{
	// Fit centers into 3D line
	Vector3 point, direction;
    MatrixXd curveCenters(nodeIDs.size(), 3);
    for (size_t r = 0; r < nodeIDs.size(); r++)
        curveCenters.row(r) = graph->getNode(nodeIDs[r])->position(Eigen::Vector4d(0.5, 0.5, 0, 0));
	point = Vector3(curveCenters.colwise().mean());
	curveCenters = curveCenters.rowwise() - point.transpose();
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(curveCenters, Eigen::ComputeThinU | Eigen::ComputeThinV);
	direction = Vector3(svd.matrixV().col(0)).normalized();

	// Sort curves
	std::vector <size_t> sorted;
	QMap<size_t, double> dists;
    for (size_t r = 0; r < nodeIDs.size(); r++) dists[r] = curveCenters.row(r).dot(direction);
	for (auto p : sortQMapByValue(dists)) sorted.push_back(p.second);

	// Build sheet control points
	Array2D_Vector3 cpnts;
	for (size_t i = 0; i < sorted.size(); i++){
		size_t idx = sorted[i];
		cpnts.push_back(graph->getNode(nodeIDs[idx])->getPoints(std::vector<Array1D_Vector4d>(1, sideCoordinates[0])).front());
	}

	// Requirment for NURBS is minimum 4 rows
	if (cpnts.size() < 4){
		if (cpnts.size() == 2)
		{
			Array1D_Vector3 m1, m2;
			for (size_t i = 0; i < cpnts.front().size(); i++) m1.push_back(AlphaBlend(1.0 / 3.0, cpnts.front()[i], cpnts.back()[i]));
			for (size_t i = 0; i < cpnts.front().size(); i++) m2.push_back(AlphaBlend(2.0 / 3.0, cpnts.front()[i], cpnts.back()[i]));
			cpnts.insert(cpnts.begin() + 1, m1);
			cpnts.insert(cpnts.begin() + 2, m2);
		}
		else
		{
			Array1D_Vector3 m1, m2;
			for (size_t i = 0; i < cpnts.front().size(); i++) m1.push_back(AlphaBlend(1.0 / 2.0, cpnts[0][i], cpnts[1][i]));
			for (size_t i = 0; i < cpnts.front().size(); i++) m2.push_back(AlphaBlend(1.0 / 2.0, cpnts[1][i], cpnts[2][i]));
			cpnts.insert(cpnts.begin() + 1, m1);
			cpnts.insert(cpnts.begin() + 3, m2);
		}
	}
	NURBS::NURBSRectangled sheet = NURBS::NURBSRectangled::createSheetFromPoints(cpnts);
    Structure::Sheet * newSheet = new Structure::Sheet(sheet, nodeIDs.front() + nodeIDs.back());

	return graph->addNode(newSheet)->id;
}

Array2D_Vector4d DeformEnergy::computeSideCoordinates( int resolution )
{
	Array2D_Vector4d coords(4);
	for (int i = 0; i < resolution; i++) coords[0].push_back(Eigen::Vector4d(double(i) / (resolution - 1), 0, 0, 0));
	for (int i = 0; i < resolution; i++) coords[1].push_back(Eigen::Vector4d(1, double(i) / (resolution - 1), 0, 0));
	for (int i = 0; i < resolution; i++) coords[2].push_back(Eigen::Vector4d(1 - (double(i) / (resolution - 1)), 1, 0, 0));
	for (int i = 0; i < resolution; i++) coords[3].push_back(Eigen::Vector4d(0, 1 - (double(i) / (resolution - 1)), 0, 0));
	return coords;
}

double DeformEnergy::deform(Structure::Node * inputNodeA, Structure::Node * inputNodeB)
{
	Structure::Node * nodeA = inputNodeA->clone();
	Structure::Node * nodeB = inputNodeB->clone();

	// Point-to-point correspondence
	if (nodeA->type() == Structure::SHEET && nodeB->type() == Structure::SHEET)
		Structure::ShapeGraph::correspondTwoSheets((Structure::Sheet*)nodeA, (Structure::Sheet*)nodeB);

	if (nodeA->type() == Structure::CURVE && nodeB->type() == Structure::CURVE)
		Structure::ShapeGraph::correspondTwoCurves((Structure::Curve*)nodeA, (Structure::Curve*)nodeB);

	// Number of sides to compute
	int numSides = (nodeA->type() == Structure::SHEET || nodeB->type() == Structure::SHEET) ? 4 : 1;

	QVector< QVector<Vector3> > quads;

	// Transformation
	{
		int num_steps = 10;

		auto cpntsA = nodeA->controlPoints();
		auto cpntsB = nodeB->controlPoints();

		std::vector<Array1D_Vector3> cur(numSides), prev(numSides);

		for (int si = 0; si < numSides; si++)
		{
			prev[si] = nodeA->getPoints(std::vector<Array1D_Vector4d>(1, sideCoordinates[si])).front();
			cur[si] = prev[si];

			for (int i = 1; i < num_steps; i++)
			{
				double t = double(i) / (num_steps - 1);

				auto pI = nodeA->getPoints(std::vector<Array1D_Vector4d>(1, sideCoordinates[si])).front();
				auto pJ = nodeB->getPoints(std::vector<Array1D_Vector4d>(1, sideCoordinates[si])).front();

				for (int u = 0; u < prev[si].size(); u++)
				{
					cur[si][u] = AlphaBlend(t, pI[u], pJ[u]);
				}

				// Build faces
				for (int i = 0; i < cur[si].size() - 1; i++)
				{
					QVector<Vector3> quad;

					quad << prev[si][i];
					quad << cur[si][i];
					quad << cur[si][i + 1];
					quad << prev[si][i + 1];

					quads << quad;
				}

				prev[si] = cur[si];
			}
		}
	}

	double area = 0.0;

	starlab::PolygonSoup * ps;
	if (debugging) ps = new starlab::PolygonSoup();

	for(auto quad : quads)
	{
		double triArea1 = (quad[1] - quad[0]).cross(quad[2] - quad[0]).norm() / 2.0;
		double triArea2 = (quad[2] - quad[0]).cross(quad[3] - quad[0]).norm() / 2.0;
		double quadArea = triArea1 + triArea2;

		area += quadArea;

		if (debugging)
		{
			QVector<starlab::QVector3> pnts;
			for (auto p : quad) pnts << p;
			ps->addPoly(pnts);
		}
	}
	if (debugging) debug << ps;

	delete nodeA;
	delete nodeB;

	return area;
}
