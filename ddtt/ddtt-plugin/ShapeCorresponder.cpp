#pragma warning(disable:4267)

#include <QDebug>
#include "next_combination.h"
#include "ShapeCorresponder.h"

#include "GenericGraph.h"
#include "StructureGraph.h"
#include "GraphDistance.h"
using namespace Structure;

#include <QStack>
QStack<double> nurbsQuality;
static void inline beginFastNURBS(){
	nurbsQuality.clear();
	nurbsQuality.push(TIME_ITERATIONS);
	nurbsQuality.push(CURVE_TOLERANCE);
	nurbsQuality.push(RombergIntegralOrder);

	TIME_ITERATIONS			= 6;
	CURVE_TOLERANCE			= 1e-05;
	RombergIntegralOrder	= 5;
}
static void inline endFastNURBS(){
	if(nurbsQuality.size() < 3) return;
	RombergIntegralOrder = nurbsQuality.pop();
	CURVE_TOLERANCE = nurbsQuality.pop();
	TIME_ITERATIONS = nurbsQuality.pop();
}

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

			double range = (dists_sorted.back() - dists_sorted.front());
			double val = (dists[n] - dists_sorted.front()) / range;

			if(range == 0.0) val = 0;

			n->meta["structure_ground"] = val;
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

	QSet<QString> toNothingPairs;

	for(auto pairs : fixedPairs){
		for(auto sid : pairs.first)
		{
			if(!sindices.contains(sid)) continue;

			int i = sindices[sid];

			// Block assignment to all
			for(int j = 0; j < M; j++)
				m[i][j] = AssignmentLib::Edge::WORST_WEIGHT;
			
			// Only allow assignment to fixed set
			for(auto tid : pairs.second)
			{
				if(!tindices.contains(tid)) 
				{
					toNothingPairs.insert(sid);
					continue;
				}

				m[i][ tindices[tid] ] = 0.0;
			}
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
	for(size_t i = 0; i < h.M.size(); i++)
	{
		auto p = h.M[i];

		if(p.first >= N || p.second >= M) continue;

		QString sid = source->nodes[p.first]->id;
		QString tid = target->nodes[p.second]->id;

		if( toNothingPairs.contains(sid) )
			tid = "NOTHING";
		
		assignments.push_back( qMakePair( qMakePair(sid, tid), m[p.first][p.second] ) );
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

VectorPairStrings pairsDebugging( QVector<Pairing> pairs )
{
	VectorPairStrings result;
	for(auto p : pairs) result.push_back( qMakePair(p.first.join(","), p.second.join(",")) );
	return result;
}

ShapeCorresponder::ShapeCorresponder(Structure::Graph * g1, Structure::Graph * g2) : source(g1), target(g2)
{
	// Set of random colors
	QVector<QColor> colors;
	for(int i = 0; i < (int)source->nodes.size() * 2; i++) colors.push_back(starlab::qRandomColor2());

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
	{
		for(auto sourceNode : source->nodes)
		{
			QVector<Pairing> fixedSets;

			// K-nearest neighbors
			{
				int k = 3;

				for(int i = 0; i < k; i++)
					fixedSets.push_back( Pairing(QStringList() << sourceNode->id, QStringList() << candidates[sourceNode->id][i].second ) );
			}

			// To nothing
			{
				fixedSets.push_back( Pairing(QStringList() << sourceNode->id, QStringList() << "NOTHING" ) );
			}

			// Prepare deformation paths
			for( auto fixedSet : fixedSets )
			{
				paths.push_back( DeformationPath() );
				DeformationPath & path = paths.back();

				QVector<Pairing> curSet;
				curSet.push_back( fixedSet );

				path.pairs = findPairing( m, curSet );
				path.pairsDebug = pairsDebugging( path.pairs );

				path.gcorr = new GraphCorresponder(source, target);

				for(auto p : path.pairs){
					QStringList sourceNodes = p.first;
					QStringList targetNodes = p.second;

					path.gcorr->addCorrespondences( sourceNodes.toVector() , targetNodes.toVector(), -1 );
				}

				path.gcorr->isReady = true;
				path.gcorr->correspondAllNodes();
			}
		}
	}

	/// Find best correspondence
	DeformationPath bestPath;
	{
		// Evaluate deformations
		#pragma omp parallel for
		for(int pi = 0; pi < (int)paths.size(); pi++)
		{
			auto & path = paths[pi];

			// Prepare blending
			path.scheduler = QSharedPointer<Scheduler>( new Scheduler );
			path.blender = QSharedPointer<TopoBlender>( new TopoBlender( path.gcorr, path.scheduler.data() ) );

			// Deform
			path.scheduler->timeStep = 0.2;
			path.scheduler->property["isDisableGrow"] = true;

			beginFastNURBS();

			path.scheduler->executeAll();

			endFastNURBS();

			// Collect error
			double error = 0;

			for(auto g : path.scheduler->allGraphs)
			{
				compute_part_measures( g );

				for(auto n_orig : source->nodes)
				{
					QVector<Node*> nodesWithID = g->nodesWithProperty("original_ID", n_orig->id);
					if(nodesWithID.size() < 1) continue;

					auto n = nodesWithID.front();
					
					double partDiff = partDifference(n_orig->id, n->id, source, g);
					if( !std::isfinite(partDiff) ) partDiff = 1e20;

					path.errors.push_back( partDiff );

					error += partDiff;
				}
			}

			// Record error
			path.weight = error;

			path.colors = colors;
		}

		// Best = lowest error
		std::sort(paths.begin(), paths.end(), DeformationPathCompare);
		bestPath = paths.front();

		// set indices
		int j = 0;
		for( auto & p : paths )	p.idx = j++;

		source->setColorAll(Qt::black);
		target->setColorAll(Qt::black);

		for( auto p : bestPath.pairs )
		{
			QColor c = starlab::qRandomColor2();
			for( auto sid : p.first ) source->setColorFor(sid, c);
			for( auto tid : p.second ) target->setColorFor(tid, c);
		}
	}
}
