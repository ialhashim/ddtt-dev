#include "StructureAnalysis.h"
#include "ShapeGraph.h"
#include "GenericGraph.h"

void StructureAnalysis::analyzeGroups(Structure::ShapeGraph * shape, bool isDebug)
{
	shape->relations.clear();

	QStringList nodesIDs;
	Structure::NodeGroups tmpA;
	for (auto n : shape->nodes) nodesIDs << n->id;
	for (auto g : shape->groups) for (auto nid : g) nodesIDs.removeAll(nid);
	for (auto nid : nodesIDs) tmpA.push_back(QVector<QString>() << nid);
	for (auto g : shape->groups) tmpA.push_back(g);
	shape->groups = tmpA;

	for (auto g : shape->groups)
	{
		Structure::Relation r;

		r.parts = QStringList::fromVector(g);

		if (g.size() == 1)
		{
			r.type = Structure::Relation::SELF;
			// TODO: figure out the plane, or select the most similar to global reflectional
		}

		// Since fitted input might not always be perfect, we equalize parts geometry (resample) here
		{
			auto firstPart = shape->getNode(r.parts.front());
			for (auto partID : r.parts)
				if (shape->getNode(partID)->numCtrlPnts() > firstPart->numCtrlPnts())
					firstPart = shape->getNode(partID);

			for (auto partID : r.parts)
			{
				auto smaller = firstPart;
				auto larger = shape->getNode(partID);
				if (smaller->numCtrlPnts() > larger->numCtrlPnts()) std::swap(smaller, larger);
				smaller->equalizeControlPoints(larger);
			}
		}

		if (g.size() == 2)
		{
			r.type = Structure::Relation::REFLECTIONAL;

			auto partA = shape->getNode(g.front()), partB = shape->getNode(g.back());

			auto plane = getReflectionalPlane(partA, partB);
			r.point = plane.first;
			r.axis = plane.second;
		}

		std::vector<Vector3> centers;			
		Vector3 centroid(0, 0, 0);

		if (g.size() > 2)
		{
			// Compute distances from group centroid to parts
			for (auto part : r.parts) centers.push_back( shape->getNode(part)->position(Eigen::Vector4d(0,0,0,0)) );

			for (auto c : centers) centroid += c;
			centroid /= centers.size();
			Array1D_Real dists;
			for (auto c : centers) dists.push_back((c-centroid).norm());

			// Median absolute deviation (MAD) 
			double avg_dist = std::accumulate(dists.begin(), dists.end(), 0.0) / dists.size();
			QVector<double> absolute_deviations;
			for (auto dist : dists) absolute_deviations << abs(dist - avg_dist);
			std::sort(absolute_deviations.begin(), absolute_deviations.end());
			double mad = 0;
			if (absolute_deviations.size() % 2 == 1)
				mad = absolute_deviations[(absolute_deviations.size() - 1) / 2];
			else
			{
				auto a = absolute_deviations[(absolute_deviations.size() / 2) - 1];
				auto b = absolute_deviations[(absolute_deviations.size() / 2)];
				mad = (a + b) / 2;
			}

			// Rotational relations have consistent distance
			double threshold = avg_dist * 0.4;
			if (mad < threshold) 
				r.type = Structure::Relation::ROTATIONAL;
			else 
				r.type = Structure::Relation::TRANSLATIONAL;

			if (r.type == Structure::Relation::TRANSLATIONAL)
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

				// Sort parts by angle around axis
				QStringList sorted;
				QMap<QString, double> angles;
				for (size_t i = 0; i < r.parts.size(); i++)
				{
					double angle = signedAngle(centers[0], centers[i], r.axis);
					if (angle < 0) angle = (M_PI * 2) + angle;
					angles[r.parts[i]] = angle;
				}
				for (auto pair : sortQMapByValue(angles)) sorted << pair.second;
				r.parts = sorted;

				// Vector from part's head to centroid
				for (size_t i = 0; i < r.parts.size(); i++)
					r.deltas.push_back(centroid - shape->getNode(r.parts[i])->position(Eigen::Vector4d(0, 0, 0, 0)));
			}
		}

		// Add relation to shape:
		shape->relations.push_back(r);

		// Visualize:
		if (isDebug)
		{
			shape->debug << starlab::PointSoup::drawPoint(r.point, 12, Qt::red);

			if (r.type != Structure::Relation::TRANSLATIONAL)
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

void StructureAnalysis::removeFromGroups(Structure::ShapeGraph * shape, Structure::Node * node)
{
	for (auto & r : shape->relations) r.parts.removeAll(node->id);
}