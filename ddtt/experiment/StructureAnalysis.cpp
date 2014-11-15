#include "StructureAnalysis.h"
#include "ShapeGraph.h"

void StructureAnalysis::analyzeGroups(Structure::ShapeGraph * shape, bool isDebug)
{
	shape->relations.clear();

	for (auto g : shape->groups)
	{
		Structure::Relation r;

		r.parts = QStringList::fromVector(g);

		if (g.size() == 1)
		{
			r.type = Structure::Relation::SELF;
			// TODO: figure out the plane, or select the most similar to global reflectional
		}

		if (g.size() == 2)
		{
			r.type = Structure::Relation::REFLECTIONAL;

			auto partA = shape->getNode(g.front()), partB = shape->getNode(g.back());

			// Since fitted input might not always be perfect, we equalize parts geometry (resample) here
			{
				if (partA->numCtrlPnts() < partB->numCtrlPnts())
					partA->equalizeControlPoints(partB);
				else
					partB->equalizeControlPoints(partA);
			}

			auto plane = getReflectionalPlane(partA, partB);
			r.point = plane.first;
			r.axis = plane.second;
		}

		std::vector<Vector3> centers;			
		Vector3 centroid(0, 0, 0);

		if (g.size() > 2)
		{
			// Compute distances from group centroid to parts
			for (auto part : r.parts) centers.push_back( shape->getNode(part)->position(Eigen::Vector4d(0.5,0.5,0,0)) );

			for (auto c : centers) centroid += c;
			centroid /= centers.size();
			Array1D_Real dists;
			for (auto c : centers) dists.push_back((c-centroid).norm());

			// Rotational relations have consistent distance
			double threshold = 0.02;
			double sigma = stdev(dists);
			if (sigma < threshold) r.type = Structure::Relation::ROTATIONAL;
			else r.type = Structure::Relation::TRANSLATION;
		}

		if (r.type == Structure::Relation::TRANSLATION)
		{
			auto line = best_line_from_points(centers);
			r.point = line.first;
			r.axis = line.second;
		}

		if (r.type == Structure::Relation::ROTATIONAL)
		{
			auto plane = best_plane_from_points(centers);
			r.point = plane.first;
			r.axis = plane.second;
		}

		// Add relation to shape:
		shape->relations.push_back(r);

		// Visualize:
		if (isDebug)
		{
			shape->debug << starlab::PointSoup::drawPoint(r.point, 12, Qt::red);

			if (r.type != Structure::Relation::TRANSLATION)
			{
				auto plane = new starlab::PlaneSoup(0.1, true, r.type == Structure::Relation::REFLECTIONAL ? Qt::red : Qt::green);
				plane->addPlane(r.point, r.axis);
				shape->debug << plane;
			}
			else
			{
				auto line = new starlab::LineSegments(3);
				line->addLine(r.point, Vector3(r.point + r.axis));
				shape->debug << line;
			}
		}
	}
}

std::pair<Vector3, Vector3> StructureAnalysis::getReflectionalPlane(Structure::Node * n1, Structure::Node * n2)
{
	Vector3 point(0, 0, 0), axis(0, 0, 0);

	int num_samples = 10;

	for (int i = 0; i < num_samples; i++)
	{
		double t = double(i) / (num_samples - 1);
		auto p = n1->position(Eigen::Vector4d(t, t, 0, 0));
		auto q = n2->position(Eigen::Vector4d(t, t, 0, 0));

		point += (p + q) * 0.5;
		axis += (p - q).normalized();
	}

	point /= num_samples;
	axis /= num_samples;
	axis.normalize();

	return std::make_pair(point, axis);
}

SurfaceMesh::Vector3 StructureAnalysis::pointReflection(const Vector3 & p, const Vector3 & planePoint, const Vector3 & planeNormal)
{
	double d = (p - planePoint).dot(planeNormal);
	Vector3 v = planeNormal * d;
	return p - (v * 2);
}
