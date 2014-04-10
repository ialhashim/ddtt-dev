#pragma warning(disable:4267)

#include <QDebug>
#include "next_combination.h"
#include "ShapeCorresponder.h"

#include "GenericGraph.h"
#include "StructureGraph.h"
#include "GraphDistance.h"
using namespace Structure;

QStringList toQStringList( const QVector<QString> & v ){
	QStringList l;
	for(auto & s : v) l << s;
	return l;
}

#include "cartesian.h"

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
		Vs(i) = sn->meta.value(keys.at(i)).toDouble();
		Vt(i) = tn->meta.value(keys.at(i)).toDouble();
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

		result.push_back( qMakePair( QVector<QString>() << pair.first, QVector<QString>() << pair.second ) );
	}

	return result;
}

VectorPairStrings pairsDebugging( QVector<Pairing> pairs )
{
	VectorPairStrings result;
	for(auto p : pairs) result.push_back( qMakePair(toQStringList(p.first).join(","), toQStringList(p.second).join(",")) );
	return result;
}

QVector<QString> nothingSet(){
	return QVector<QString>() << "NOTHING";
}

VectorPairings splitMany( QVector<QString> A, QVector<QString> B )
{
	VectorPairings result;

	// Pick first and last from smaller set
	QPair< QString, QString > smaller( A.front(), A.back() );

	// Divide larger into two
	QPair< QVector<QString>, QVector<QString> > larger;
	std::vector<QString> b1(B.begin(), B.begin() + B.size()/2);
	std::vector<QString> b2(B.begin() + B.size()/2, B.end());
	for(auto n : b1) larger.first.push_back(n);
	for(auto n : b2) larger.second.push_back(n);

	result.push_back( qMakePair(QVector<QString>() << smaller.first, larger.first) );
	result.push_back( qMakePair(QVector<QString>() << smaller.second, larger.second) );

	return result;
}

QVector<VectorPairings> manyToMany(Structure::Graph * sg, Structure::Graph * tg, QVector<QString> snodes, QVector<QString> tnodes)
{
	QVector<VectorPairings> allResolved;

	// Sort based on 'x' axis
	QVector< QPair<double, Structure::Node*> > sortedS, sortedT;
	for(auto nid : snodes) {
		Structure::Node * n = sg->getNode(nid); 
		sortedS.push_back( qMakePair(n->meta.contains("geo_x") ? n->meta["geo_x"].toDouble() : sortedS.size(), n) );
	}
	for(auto nid : tnodes) {
		Structure::Node * n = tg->getNode(nid); 
		sortedT.push_back( qMakePair(n->meta.contains("geo_x") ? n->meta["geo_x"].toDouble() : sortedT.size(), n) );
	}
	qSort(sortedS);
	qSort(sortedT);

	bool isReversed = false;

	QVector<QString> sortedA, sortedB;
	for(auto nv : sortedS) sortedA.push_back(nv.second->id);
	for(auto nv : sortedT) sortedB.push_back(nv.second->id);

	if(sortedS.size() > sortedT.size()){
		isReversed = true;
		std::swap(sortedA, sortedB);
	}

	// Match using 'splitMany'
	{
		VectorPairings resolved;
		VectorPairings split = splitMany(sortedA, sortedB);
		if( isReversed ) for(auto & s : split) std::swap(s.first, s.second);
		for(auto p : split) resolved.push_back(p);
		allResolved.push_back(resolved);
	}

	// Try a one-to-one
	if(sortedA.size() == sortedB.size())
	{
		VectorPairings resolved;

		for(auto sid: sortedA){
			//Array1D_Vector3 sPoints = sg->getNode(sid)->controlPoints();
			
			Vector3 ps = sg->getNode(sid)->controlPoints().front();

			double minDist = DBL_MAX;
			QString bestID;

			for(auto tid: sortedB){
				//Array1D_Vector3 tPoints = tg->getNode(tid)->controlPoints();
				//double dist = HausdorffDistance(sPoints, tPoints);
				Vector3 pt = tg->getNode(tid)->controlPoints().front();
				double dist = (ps - pt).norm();

				// Closest
				if(dist < minDist){
					minDist = dist;
					bestID = tid;
				}
			}

			resolved.push_back( Pairing( QVector<QString>()<< sid, QVector<QString>()<< bestID ) );
		}

		allResolved.push_back( resolved );
	}

	return allResolved;
}

Assignments allAssignments( QVector< QVector<QString> > sgroups, QVector< QVector<QString> > tgroups, 
						   mat m, Structure::Graph * sg, Structure::Graph * tg )
{
	Assignments assignments;

	// Map to index
	QMap< int, QVector<QString> > sourceItem, targetItem;
	for(auto items : sgroups) sourceItem[sourceItem.size()] = items;
	for(auto items : tgroups) targetItem[targetItem.size()] = items;
	sourceItem[sourceItem.size()] = nothingSet();
	targetItem[targetItem.size()] = nothingSet();

	// Record all initial costs
	QMap< int, QVector< QPair<double, int> > > candidates;
	for(size_t i = 0; i < sgroups.size(); i++)
	{
		QVector< QPair<double, int> > diffs;
		for(size_t j = 0; j < tgroups.size(); j++) diffs.push_back( qMakePair( m[i][j], j ) );
		std::sort(diffs.begin(), diffs.end());
		std::reverse(diffs.begin(), diffs.end());
		candidates[i] = diffs;
	}

	// Collect only 'k' nearest candidates + nothing
	int K = 3;
	
	std::map< int, std::set<int> > allowed;
	for(size_t i = 0; i < sgroups.size(); i++){
		allowed[i] = std::set<int>();
		for(int j = 0; j < K; j++) allowed[i].insert(candidates[i][j].second); // Good candidates
		allowed[i].insert( tgroups.size() ); // Nothing
	}

	// Build full set of good assignments
	std::vector< std::vector<int> > goodCombs;
	{
		// Flatten sets
		std::vector< std::vector<int> > setAllowed;

		for(auto a : allowed){
			std::vector<int> s;
			for(auto e : a.second) s.push_back(e);
			setAllowed.push_back(s);
		}

		int NOTHING_ITEM = tgroups.size();

		blitz::CartesianProduct< std::vector<int>, std::vector<int> > product( setAllowed );

		// Go over all possible combinations, only keep good ones (i.e. single matching to targets)
		int totalCombs = 0;
		for(auto comb : product){
			totalCombs++;
			bool isGood = true;
			std::set<int> titems;
			for(auto i : comb) isGood = isGood && ((i == NOTHING_ITEM) || titems.insert(i).second);
			if( isGood ) goodCombs.push_back( comb );
		}
	}

	// Generate final assignments
	for(auto assignment : goodCombs)
	{
		VectorPairings paring;
		QVector<Pairing> manyManyCases;

		for(size_t i = 0; i < assignment.size(); i++)
		{
			QVector<QString> snodes = sourceItem[i];
			QVector<QString> tnodes = targetItem[assignment[i]];

			// Add all but many-to-many
			if(snodes.size() > 1 && tnodes.size() > 1)
				manyManyCases.push_back( qMakePair(snodes, tnodes) );
			else
				paring.push_back( qMakePair(snodes,tnodes) );
		}

		// Now resolve many-to-many using different possibilities
		if( manyManyCases.size() )
		{
			QVector< QVector<VectorPairings> > manyManyResolution;

			for(auto candidate : manyManyCases)
				manyManyResolution.push_back( manyToMany(sg, tg, candidate.first, candidate.second) );
			
			blitz::CartesianProduct< QVector<VectorPairings>, QVector<VectorPairings> > product( manyManyResolution );

			for(auto comb : product){
				VectorPairings curParing = paring;
				for(auto pairs : comb) 
					for(auto pair : pairs) 
						curParing.push_back(pair);
				
				assignments.push_back( curParing );
			}
		}
		else
			assignments.push_back( paring );
	}

	return assignments;
}

ShapeCorresponder::ShapeCorresponder(Structure::Graph * g1, Structure::Graph * g2, bool) : source(g1), target(g2)
{
	// Set of random colors
	QVector<QColor> colors;
	for(int i = 0; i < (int)source->nodes.size() * 2; i++) colors.push_back(starlab::qRandomColor2());

	QElapsedTimer prepareTimer, computeTimer;
	prepareTimer.start();

	QVector<Structure::Graph*> graphs;
	graphs << source << target;

	/// Graphs preprocessing
	for(auto g : graphs) compute_part_measures( g );

	mat m = buildDifferenceMatrix(source, target);

	Assignments assignments = allAssignments(source->nodesAsGroups(), target->nodesAsGroups(), m, source, target);

	// Ground-truth
	{
		// Check for existing correspondence file
		QString path = QFileInfo(source->property["name"].toString()).absolutePath() + "/", ext = ".txt";
		QString g1n = source->name(), g2n = target->name();
		QStringList files; files << (path+g1n+"_"+g2n+ext) << (path+g2n+"_"+g1n+ext);
		int fileidx = -1;
		for(int i = 0; i < 2; i++){
			QString f = files[i];
			if(!QFileInfo (f).exists()) continue;
			fileidx = i;
		}
		if( fileidx != -1 ){

			bool corrReversed = (fileidx == 0) ? false : true;
			GraphCorresponder * bestCorrespond = new GraphCorresponder(source, target);
			bestCorrespond->loadCorrespondences(files[fileidx], corrReversed);
			
			VectorPairings v1;
			QVector< PART_LANDMARK > vp = bestCorrespond->correspondences;
			for(auto p : vp) v1.push_back( Pairing(p.first, p.second) );

			bool isFound = false;

			for( VectorPairings v2 : assignments ){
				if(v1.size() == v2.size() && std::is_permutation(v1.begin(), v1.end(), v2.begin())){
					isFound = true;
					break;
				}
			}
		}
	}
	
	// Prepare deformation paths
	for( auto a : assignments )
	{
		paths.push_back( DeformationPath() );
		DeformationPath & path = paths.back();

		path.pairs = a;
		path.pairsDebug = pairsDebugging( path.pairs );

		path.gcorr = new GraphCorresponder(source, target);

		int i = 0;

		for(auto p : path.pairs)
		{
			path.gcorr->addCorrespondences( p.first, p.second, -1 );

			// Nodes coloring
			for(auto sid : p.first) path.scolors[sid] = colors[i];
			for(auto tid : p.second) path.tcolors[tid] = colors[i];
			i++;
		}

		path.gcorr->isReady = true;
		path.gcorr->correspondAllNodes();
	}

	// Timing
	property["prepareTime"].setValue( (int)prepareTimer.elapsed() );
	computeTimer.start();

	/// Find best correspondence
	DeformationPath bestPath;
	{
		beginFastNURBS();

		// Evaluate deformations
		#pragma omp parallel for
		for(int pi = 0; pi < (int)paths.size(); pi++)
		{
			auto & path = paths[pi];

			// Prepare blending
			path.scheduler = QSharedPointer<Scheduler>( new Scheduler );
			path.blender = QSharedPointer<TopoBlender>( new TopoBlender( path.gcorr, path.scheduler.data() ) );

			// Deform
			path.scheduler->timeStep = 0.1;
			path.scheduler->executeAll();

			// Collect error
			double error = 0;
			
			for(auto g : path.scheduler->allGraphs)
			{
				g->moveBottomCenterToOrigin();

				compute_part_measures( g );

				for(auto n_orig : source->nodes)
				{
					Structure::Node * n = NULL;
					for(auto ni : g->nodes)
						if(ni->property.contains("original_ID") && ni->property.value("original_ID") == n_orig->id)
							n = ni;
					if(!n) continue;

					double partDiff = partDifference(n_orig->id, n->id, source, g);
					if( !std::isfinite(partDiff) ) partDiff = 1e20;

					path.errors.push_back( partDiff );

					error += partDiff;
				}
			}

			// Record error
			path.weight = error;

			// Clean up
			{
				path.scheduler.clear();
				path.blender.clear();
				path.errors.clear();
			}
		}
	}

	// Timing
	property["computeTime"].setValue( (int)computeTimer.elapsed() );

	if( paths.size() )
	{
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

	// Stat
	property["pathsCount"].setValue( (int)paths.size() );
}
