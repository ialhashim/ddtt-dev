#include "PropagateProximity.h"
#include "ShapeGraph.h"

#include "convexhull2d.h"
#include "mean_value_coordinates.h"

Q_DECLARE_METATYPE(Vector3);

struct ProximityConstraint{
    Vector3 d;
    Structure::Node *from, *to;
    Structure::Link *link;
    ProximityConstraint(Structure::Node *from = nullptr, Structure::Link *link = nullptr) :from(from), link(link){
        if (!from || !link) return;
        to = link->otherNode(from->id);
        d = link->property[from->id].value<Vector3>();
    }
    inline Vector3 start(){ return link->position(from->id); }
	inline Vector3 start2(){ return link->getNode(from->id)->position(coords().back()); }
    inline Vector3 delta(){ return d; }
    inline Eigen::Vector4d coord(){ return link->getCoord(to->id).front(); }
	inline Array1D_Vector4d coords(){ return link->getCoord(to->id); }
};

void PropagateProximity::prepareForProximity(Structure::Graph * graph)
{
    for (auto & edge : graph->edges){
        auto p1 = edge->position(edge->n1->id), p2 = edge->position(edge->n2->id);
        edge->property[edge->n1->id].setValue(Vector3(p2 - p1));
        edge->property[edge->n2->id].setValue(Vector3(p1 - p2));
    }
}

void PropagateProximity::propagateProximity(const QStringList &fixedNodes, Structure::ShapeGraph *graph)
{
    // Constraints per part
    QMap < QString, QVector< ProximityConstraint > > constraints;

    // Initialize propagation state
    for (auto n : graph->nodes){
        n->property["propagated"].setValue( false );
        n->property["fixed"].setValue( fixedNodes.contains(n->id) );
    }

    /// Find propagation levels:
    // First level:
    QVector < QVector<Structure::Node*> > propagationLevel(1);
    for (auto nid : fixedNodes) {
        auto n = graph->getNode(nid);
        n->property["propagated"].setValue(true);
        propagationLevel.front().push_back(n);
    }
    // Remaining levels:
    forever{
        QVector<Structure::Node*> curLevel;
        for (auto & nid : propagationLevel.back())
        {
            for (auto & edge : graph->getEdges(nid->id))
            {
                if (edge->n1->id == edge->n2->id) continue;
                auto otherNode = edge->otherNode(nid->id);
                if (otherNode->property["propagated"].toBool()) continue;

                if (!curLevel.contains(otherNode)) curLevel.push_back(otherNode);

                constraints[otherNode->id].push_front( ProximityConstraint(nid,edge) );
            }
        }

        for (auto & nid : curLevel)	nid->property["propagated"].setValue(true);
        if (curLevel.isEmpty())	break;
        propagationLevel.push_back(curLevel);
    };

    // Apply constraints
    propagationLevel.removeFirst(); // Fixed parts do not change
    for (int i = 0; i < propagationLevel.size(); i++)
    {
        for (auto & n : propagationLevel[i])
        {
            assert(constraints[n->id].size());

            if (constraints[n->id].size() == 1)
			{
				auto & c = constraints[n->id].front();

				if (c.link->type == Structure::POINT_EDGE)
					n->deformTo(c.coord(), c.start() + c.delta(), true);

				if (c.link->type == Structure::LINE_EDGE)
					n->deformTwoHandles(c.coords().front(), c.start() + c.delta(), c.coords().back(), c.start2() + c.delta());
            }
            else
            {
                auto & c_list = constraints[n->id];

                if (c_list.size() == 2)
                {
                    auto & ca = c_list.front();
                    auto & cb = c_list.back();

                    n->deformTwoHandles(ca.coord(), ca.start() + ca.delta(), cb.coord(), cb.start() + cb.delta());
                }
                else
				{
                    // Only consider outer most constraints
                    QVector<ProximityConstraint> filtered_constraints;
                    std::vector<Eigen::Vector2d> coords;
                    for (auto & c : c_list) coords.push_back(Eigen::Vector2d(c.coord()[0], c.coord()[1]));
					for (auto idx : chull2d::convex_hull_2d(coords)) filtered_constraints << c_list[idx];

					bool isFixedAround = true;
					for (auto & c : c_list) isFixedAround &= fixedNodes.contains(c.from->id);
					if (isFixedAround)
					{
						filtered_constraints = c_list;                    
						
						std::vector < Eigen::Vector3d > cage;
						for (auto & c : filtered_constraints) cage.push_back(c.link->position(n->id));

						// Rotate cage towards z-axis around center of constraints
						auto best_plane = best_plane_from_points(cage);
						double cage_plane_dot = best_plane.second.dot(Vector3::UnitZ());
						if (cage_plane_dot < 0) best_plane.second *= 1;
						Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(best_plane.second, Vector3::UnitZ());
						for (auto & p : cage) p = (q * (p - best_plane.first)) + best_plane.first;  // now a rotated cage

						QVector<ProximityConstraint> ordered_constraints;
						auto ordered_indices = chull2d::convex_hull_2d(cage);
						for (auto idx : ordered_indices) ordered_constraints << filtered_constraints[idx];
						filtered_constraints = ordered_constraints;

						// DEBUG: cage
						if (false){
							auto ps = new starlab::PolygonSoup;
							QVector<starlab::QVector3> pts;
							for (auto p : ordered_indices) pts << cage[p];
							ps->addPoly(pts);
							graph->debug << ps;
						}
					}

					if (filtered_constraints.size() == 2)
					{
						auto & ca = c_list.front();
						auto & cb = c_list.back();
						n->deformTwoHandles(ca.coord(), ca.start() + ca.delta(), cb.coord(), cb.start() + cb.delta());
						continue;
					}

                    // Build cage from convex hull of constraints
                    std::vector < Eigen::Vector3d > cage;
                    for (auto & c : filtered_constraints) cage.push_back(c.link->position(n->id));

                    // Translation
                    Vector3 oldCenter(0, 0, 0), newCenter(0, 0, 0);
                    for (auto & p : cage) oldCenter += p;
                    oldCenter /= cage.size();

                    // Rotate cage towards z-axis around center of constraints
                    auto best_plane = best_plane_from_points(cage);
					double cage_plane_dot = best_plane.second.dot(Vector3::UnitZ());
					if (cage_plane_dot < 0) best_plane.second *= 1;
					Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(best_plane.second, Vector3::UnitZ());
					for (auto & p : cage) p = (q * (p - best_plane.first)) + best_plane.first;  // now a rotated cage

                    // Compute weights and 'heights' for control points
                    auto cpts = n->controlPoints();
                    std::vector< std::vector<double> > weights;
                    std::vector<double> heights;

                    for (auto & p : cpts)
                    {
                        // rotate control points towards z-axis using same above rotation
                        p = (q * (p - best_plane.first)) + best_plane.first;

                        // Record "height"
                        heights.push_back(p.z());

                        // Compute weights
                        weights.push_back(MeanValueCoordinates::computeWeights(p[0], p[1], cage));
                    }

                    // Modify rotated cage to expected positions
                    for (size_t i = 0; i < cage.size(); i++)
                    {
                        Vector3 pj = filtered_constraints[i].start() + filtered_constraints[i].delta();
                        Vector3 p = (q * (pj - best_plane.first)) + best_plane.first;
                        cage[i][0] = p[0];
                        cage[i][1] = p[1];

                        newCenter += pj;
                    }

                    newCenter /= cage.size();
                    Vector3 translation = newCenter - oldCenter;

                    // Compute new locations along z-plane:
                    for (size_t i = 0; i < cpts.size(); i++){
                        auto p = MeanValueCoordinates::interpolate2d(weights[i], cage);
                        cpts[i][0] = p.first;
                        cpts[i][1] = p.second;
                        cpts[i][2] = heights[i];
                    }

                    // Rotate back
                    auto qinv = q.inverse();
                    for (auto & p : cpts) p = translation + ((qinv * (p - best_plane.first)) + best_plane.first);

                    n->setControlPoints(cpts);
                }
            }
        }

		graph->saveKeyframe();
    }
}
