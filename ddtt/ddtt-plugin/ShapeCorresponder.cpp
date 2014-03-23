#pragma warning(disable:4267)

#include <QDebug>
#include "next_combination.h"
#include "ShapeCorresponder.h"

#include "GenericGraph.h"
#include "StructureGraph.h"
#include "GraphDistance.h"
using namespace Structure;

#include "GraphCorresponder.h"
#include "Task.h"
#include "Scheduler.h"
#include "TopoBlender.h"

struct DeformationPath{
	Structure::Graph *source, *target;
	QVector<Pairing> pairs;
	GraphCorresponder * gcorr;
	QSharedPointer<Scheduler> scheduler;
	QSharedPointer<TopoBlender> blender;
	double weight;
	DeformationPath(){ source = target = NULL; gcorr = NULL; weight = 0.0; }
};
bool DeformationPathCompare (const DeformationPath & i, const DeformationPath & j) { return (i.weight < j.weight); }

void measure_position(Structure::Graph * graph)
{
	Eigen::AlignedBox3d bbox = graph->bbox();
	Vector3 bbox_min = bbox.min();
	Vector3 bbox_max = bbox.max();

	foreach(Node * n, graph->nodes)
	{
		if(n->property.contains("taskType") && n->property["taskType"].toInt() == Task::GROW) continue;

		Vector3 p = n->bbox().center();
		std::vector<Vector3> curve = n->discretizedAsCurve( DIST_RESOLUTION );

		n->meta["geo_x"] = (p[0] - bbox_min[0]) / (bbox_max[0] - bbox_min[0]);
		n->meta["geo_y"] = (p[1] - bbox_min[1]) / (bbox_max[1] - bbox_min[1]);
		n->meta["geo_z"] = (p[2] - bbox_min[2]) / (bbox_max[2] - bbox_min[2]);
	}
}

void measure_ground_parts(Structure::Graph * graph)
{
	/// Find set of ground nodes
	std::vector<Structure::Node*> floorNodes;
	{
		Eigen::AlignedBox3d box = graph->bbox( true );
		double bottom = box.min().z();
		double height = box.max().z() - bottom;

		double threshold = height * 0.1;

		foreach(Node * n, graph->nodes)
		{
			double minz = n->bbox().min().z();
			double dist = abs(minz - bottom);

			if( dist > threshold )
				continue;
			else
				floorNodes.push_back(n);
		}
	}

	/// Build a graph
	GenericGraphs::Graph<size_t, double> g;
	int floorNode = graph->nodes.size() * 10;
	{
		// Add regular edges
		for(Link * e : graph->edges){
			g.AddEdge(e->n1->property["index"].toUInt(), e->n2->property["index"].toUInt(), 1, e->property["uid"].toInt());
		}

		for(Node * n : floorNodes){
			g.AddEdge(floorNode, n->property["index"].toUInt(), 0);
		}
	}

	/// Assign based on distance to ground
	{
		QMap<Node*, double> dists;

		foreach(Node * n, graph->nodes)	
			dists[n] = g.NodeDistance(floorNode, n->property["index"].toUInt());

		QList<double> dists_sorted = dists.values();
		qSort(dists_sorted);

		foreach(Node * n, graph->nodes)
		{
			if(n->property.contains("taskType") && n->property["taskType"].toInt() == Task::GROW) continue;

			std::vector<Vector3> curve = n->discretizedAsCurve( DIST_RESOLUTION );
			double val = (dists[n] - dists_sorted.front()) / (dists_sorted.back() - dists_sorted.front());
			n->meta["structure_ground"] = QString::number(val);
		}
	}
}

void compute_part_measures(Structure::Graph * graph)
{
	measure_ground_parts( graph );
	measure_position( graph );
}

double partDifference(QString sid, QString tid, Structure::Graph * source, Structure::Graph * target )
{
	std::vector<double> sim( target->nodes.size(), DBL_MAX );

	Structure::Node * sn = source->getNode(sid);
	Structure::Node * tn = target->getNode(tid);

	int k = sn->meta.size();
	QList<QString> keys = sn->meta.keys();
	QList<QString> tkeys = tn->meta.keys();
	keys = (keys.toSet().intersect(tkeys.toSet())).toList(); // make sure only test shared attributes

	// Source vector
	Eigen::VectorXd Vs( k ), Vt ( k );

	for(int i = 0; i < k; i++){
		Vs(i) = sn->meta[ keys.at(i) ].toDouble();
		Vt(i) = tn->meta[ keys.at(i) ].toDouble();
	}

	return (Vs - Vt).norm();
}

mat buildDifferenceMatrix( Structure::Graph * source, Structure::Graph * target )
{
	int N = source->nodes.size();
	int M = target->nodes.size();

	int extra_N = (M > N) ? M-N : 0;
	int extra_M = (N > M) ? N-M : 0;

	mat m( N + extra_N, mat_row(M + extra_M, AssignmentLib::Edge::WORST_WEIGHT) );

	for(int i = 0; i < N; i++)
		for(int j = 0; j < M; j++)
			m[i][j] = -partDifference( source->nodes[i]->id, target->nodes[j]->id, source, target );

	return m;
}

QVector<Pairing> ShapeCorresponder::findPairing( mat m, QVector<Pairing> fixedPairs )
{
	QVector<Pairing> result;

	/// Build weights matrix
	AssignmentLib::Matrix mm( m.size(), std::vector<AssignmentLib::Edge>( m.front().size() ));
	
	int N = source->nodes.size();
	int M = target->nodes.size();

	/// Fixed weights
	QMap<QString, int> sindices, tindices;
	for(auto n : source->nodes) sindices[n->id] = n->property["index"].toUInt();
	for(auto n : target->nodes) tindices[n->id] = n->property["index"].toUInt();

	for(auto pairs : fixedPairs){
		for(auto sid : pairs.first)
		{
			int i = sindices[sid];

			// Block assignment to all
			for(int j = 0; j < M; j++)
				m[i][j] = AssignmentLib::Edge::WORST_WEIGHT;
			
			// Only allow assignment to fixed set
			for(auto tid : pairs.second)
				m[i][ tindices[tid] ] = 0.0;
		}
	}

	/// Setup bipartite graph matrix
	for(int i = 0; i < N; i++)
		for(int j = 0; j < M; j++)
			mm[i][j] = AssignmentLib::Edge( pair<size_t, size_t>(i,j), m[i][j] );

	AssignmentLib::BipartiteGraph bg( mm );
	AssignmentLib::Hungarian h(bg);
	h.HungarianAlgo();

	QVector< QPair< QPair<QString,QString>, double > > assignments;
	for( auto p : h.M ){
		if(p.first >= N || p.second >= M) continue;
		assignments.push_back( qMakePair( qMakePair(source->nodes[p.first]->id, target->nodes[p.second]->id), m[p.first][p.second] ) );
	}

	// Record found pairings
	for(auto a : assignments)
	{
		QPair<QString,QString> pair = a.first;
		//double cost = a.second;

		result.push_back( qMakePair( QStringList() << pair.first, QStringList() << pair.second ) );
	}

	return result;
}

ShapeCorresponder::ShapeCorresponder(Structure::Graph * g1, Structure::Graph * g2) : source(g1), target(g2)
{
	QVector<Structure::Graph*> graphs;
	graphs << source << target;

	/// Graphs preprocessing
	for(auto g : graphs) compute_part_measures( g );

	mat m = buildDifferenceMatrix(source, target);

	QMap< QString, QVector< QPair<double, QString> > > candidates;
	for(size_t i = 0; i < source->nodes.size(); i++)
	{
		Node * sn = source->nodes[i];

		// Collect candidates
		QVector< QPair<double, QString> > diffs;
		for(size_t j = 0; j < target->nodes.size(); j++)
			diffs.push_back( qMakePair( m[i][j], target->nodes[j]->id ) );
		
		std::sort(diffs.begin(), diffs.end());
		std::reverse(diffs.begin(), diffs.end());

		candidates[sn->id] = diffs;
	}

	/// Build set of candidate correspondences
	std::vector<DeformationPath> paths;
	{
		for(auto sourceNode : source->nodes)
		{
			int k = 3;

			for(int i = 0; i < k; i++)
			{
				QVector<Pairing> fixedSet;

				fixedSet.push_back( Pairing(QStringList() << sourceNode->id, QStringList() << candidates[sourceNode->id][i].second ) );

				paths.push_back( DeformationPath() );
				DeformationPath & path = paths.back();

				// Set correspondence
				{
					path.pairs = findPairing( m, fixedSet );

					path.gcorr = new GraphCorresponder(source, target);

					for(auto p : path.pairs){
						QStringList sourceNodes = p.first;
						QStringList targetNodes = p.second;

						path.gcorr->addCorrespondences( sourceNodes.toVector() , targetNodes.toVector(), -1 );
					}

					path.gcorr->isReady = true;
					path.gcorr->correspondAllNodes();
				}

				// Prepare blending
				{
					path.scheduler = QSharedPointer<Scheduler>( new Scheduler );
					path.blender = QSharedPointer<TopoBlender>( new TopoBlender( path.gcorr, path.scheduler.data() ) );
				}
			}
		}
	}

	/// Find best correspondence
	{
		// Evaluate deformations
		for(auto & path : paths)
		{
			// Deform
			path.scheduler->timeStep = 0.1;
			path.scheduler->property["isDisableGrow"] = true;
			path.scheduler->executeAll();

			// Collect error
			double error = 0;

			for(auto g : path.scheduler->allGraphs)
			{
				compute_part_measures( g );

				for(auto n_orig : source->nodes)
				{
					auto n = g->nodesWithProperty("original_ID", n_orig->id).front();
					if(n == NULL) continue;

					error += partDifference(n_orig->id, n->id, source, g);
				}
			}

			// Record error
			path.weight = error;
		}

		std::sort(paths.begin(), paths.end(), DeformationPathCompare);

		DeformationPath & bestPath = paths.front();
	}
}
