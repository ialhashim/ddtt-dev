#include "Propagate.h"
#include "StructureGraph.h"
#include "GenericGraph.h"
Q_DECLARE_METATYPE(Vector3);

#include "convexhull2d.h"

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
	inline Vector3 delta(){ return d; }
	inline Eigen::Vector4d coord(){ return link->getCoord(to->id).front(); }
};

void Propagate::prepareForProximity(Structure::Graph * graph)
{
	for (auto & edge : graph->edges){
		auto p1 = edge->position(edge->n1->id), p2 = edge->position(edge->n2->id);
		edge->property[edge->n1->id].setValue(Vector3(p2 - p1));
		edge->property[edge->n2->id].setValue(Vector3(p1 - p2));
	}
}

void Propagate::propagateProximity(const QStringList &fixedNodes, Structure::Graph *graph)
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
				n->deformTo(c.coord(), c.start() + c.delta(), true);
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
					QMap < double, QPair<ProximityConstraint, ProximityConstraint> > furthestConstraints;

					for (auto & ci : c_list)
						for (auto & cj : c_list)
							furthestConstraints[(ci.coord() - cj.coord()).norm()] = qMakePair(ci,cj);

					auto selectedTwo = furthestConstraints.values().back();

					auto & ca = selectedTwo.first;
					auto & cb = selectedTwo.second;
					n->deformTwoHandles(ca.coord(), ca.start() + ca.delta(), cb.coord(), cb.start() + cb.delta());
					
					/*
					// Only consider outer most constraints
					QVector<ProximityConstraint> filtered_constraints;
					QVector<Eigen::Vector2d> coords;
					for (auto & c : c_list) coords.push_back(Eigen::Vector2d(c.coord()[0], c.coord()[1]));
					for (auto idx : convexhull2d_indices(coords)) filtered_constraints << c_list[idx];

					filtered_constraints.size();*/
				}
			}
		}
	}
}

void Propagate::propagateSymmetry(const QStringList &fixedNodes, Structure::Graph *graph)
{

}
