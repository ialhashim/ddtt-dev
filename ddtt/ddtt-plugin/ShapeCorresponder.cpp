#pragma warning(disable:4267)

#include <random>
#include <algorithm>
#include <iterator>

#include <QApplication>
#include <QTableWidget>
#include <QHeaderView>
#include <QDebug>
#include "next_combination.h"
#include "ShapeCorresponder.h"
#include "SynthesisManager.h"

#include "GenericGraph.h"
#include "StructureGraph.h"
#include "GraphDistance.h"
using namespace Structure;

static inline QString shortName(QString name){
	if(name.length() < 3) return name;
	return QString("%1%2%3").arg(name.at(0)).arg(name.at(1)).arg(name.at(name.length()-1));
}

QStringList toQStringList( const QVector<QString> & v ){
	QStringList l;
	for(auto & s : v) l << s;
	return l;
}

#include "cartesian.h"

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
	//measure_ground_parts( graph );
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

double partDifference2(QString sid, QString tid, Structure::Graph * source, Structure::Graph * target )
{
	Array1D_Vector3 sPoints = source->getNode(sid)->controlPoints();
	Array1D_Vector3 tPoints = target->getNode(tid)->controlPoints();
	return HausdorffDistance(sPoints, tPoints);
}

mat buildDifferenceMatrix( Structure::Graph * source, Structure::Graph * target )
{
	int N = source->nodes.size();
	int M = target->nodes.size();

	int extra_N = (M > N) ? M-N : 0;
	int extra_M = (N > M) ? N-M : 0;

	mat m( N + extra_N, mat_row(M + extra_M, DBL_MAX) );

	double minVal = DBL_MAX, maxVal = -DBL_MAX;

	for(int i = 0; i < N; i++){
		for(int j = 0; j < M; j++){
			double val = partDifference2( source->nodes[i]->id, target->nodes[j]->id, source, target );
			m[i][j] = val;

			// Track limits
			minVal = std::min(minVal, val);
			maxVal = std::max(maxVal, val);
		}
	}

	// Normalize
	for(int i = 0; i < (int)m.size(); i++){
		for(int j = 0; j < (int)m.front().size(); j++){
			if(m[i][j] == DBL_MAX) m[i][j] = maxVal;
			m[i][j] = (m[i][j] - minVal) / (maxVal - minVal);
		}
	}

	return m;
}

void visualizeDifferenceMatrix(mat m, Structure::Graph * sg, Structure::Graph * tg)
{
	int rows = sg->nodes.size(), cols = tg->nodes.size();
	QTableWidget * tw = new QTableWidget(rows, cols);

	// Cells
	int fixedSize = 40;
	tw->horizontalHeader()->setDefaultSectionSize(fixedSize);
	tw->verticalHeader()->setDefaultSectionSize(fixedSize);
	tw->horizontalHeader()->setSectionResizeMode(QHeaderView::Fixed);
	tw->verticalHeader()->setSectionResizeMode(QHeaderView::Fixed);

	// Labels
	QStringList hlabels, vlabels;
	for(auto n : sg->nodes) vlabels << shortName(n->id); tw->setVerticalHeaderLabels(vlabels);
	for(auto n : tg->nodes) hlabels << shortName(n->id); tw->setHorizontalHeaderLabels(hlabels);
	
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){
			double val = abs( m[i][j] );
			tw->setItem(i, j, new QTableWidgetItem( QString::number(val, 'g', 1) ));
			QColor c = starlab::qtJetColor( val );
			tw->item(i,j)->setBackground( starlab::qtJetColor( val ) );
			tw->item(i,j)->setTextColor( c.lighter() );
		}
	}

	tw->show();
	qApp->processEvents();
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
						   mat m, Structure::Graph * sg, Structure::Graph * tg, int K)
{
	Assignments assignments;

	// Map to index
	std::map< int, QVector<QString> > sourceItem, targetItem;
	for(auto items : sgroups) sourceItem[sourceItem.size()] = items;
	for(auto items : tgroups) targetItem[targetItem.size()] = items;
	sourceItem[sourceItem.size()] = nothingSet();
	targetItem[targetItem.size()] = nothingSet();

	// Bounds
	K = qMin(K, tgroups.size());

	// Record all initial costs
	std::map< int, QVector< QPair<double, int> > > candidates;
	std::map< QString, QString > candidates_debug;
	for(size_t i = 0; i < sgroups.size(); i++)
	{
		QVector< QPair<double, int> > diffs;

		// Score based on minimum matching of the source and target group
		for(size_t jj = 0; jj < tgroups.size(); jj++)
		{
			double minVal = DBL_MAX;
			int si = 0;
			int tj = 0;

			for(auto sid : sgroups[i])
			{
				int i = sg->indexOfNode( sg->getNode(sid) );

				for(auto tid : tgroups[jj])
				{
					int j = tg->indexOfNode( tg->getNode(tid) );
					double val = abs(m[i][j]);

					if(val < minVal){
						si = i;
						tj = j;
						minVal = val;
					}
				}
			}

			diffs.push_back( qMakePair( m[si][tj], jj ) );
		}

		std::sort(diffs.begin(), diffs.end());
		//std::reverse(diffs.begin(), diffs.end());
		candidates[i] = diffs;

		// DEBUG:
		QString curGroup = toQStringList(sgroups[i]).join("");
		for(size_t j = 0; j < tgroups.size(); j++){
			QString targetGroup = tgroups[diffs[j].second].front();
			candidates_debug[ curGroup ] += targetGroup + "-";
		}
	}

	// Collect only 'k' nearest candidates + nothing
	std::map< int, std::vector<int> > allowed;
	for(size_t i = 0; i < sgroups.size(); i++){
		allowed[i] = std::vector<int>();
		for(int j = 0; j < K; j++) allowed[i].push_back(candidates[i][j].second); // Good candidates
		allowed[i].push_back( tgroups.size() ); // Nothing
	}

	// DEBUG:
	std::map<QString, QString> allowed_debug;
	for(size_t i = 0; i < sgroups.size(); i++){
		for(int j = 0; j < K; j++)
			allowed_debug[ toQStringList(sgroups[i]).join("") ] += tgroups[ allowed[i][j] ].front() + "-";
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
	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(colors.begin(), colors.end(), g);

	QElapsedTimer prepareTimer, computeTimer;
	prepareTimer.start();

	QVector<Structure::Graph*> graphs;
	graphs << source << target;

	/// Graphs preprocessing
	for(auto g : graphs) compute_part_measures( g );

	mat m = buildDifferenceMatrix(source, target);
	if(true) visualizeDifferenceMatrix(m, source, target);

	// Considered neighbors
	int K = 2;
	bool isSubsample = false;

	Assignments assignments = allAssignments(source->nodesAsGroups(), target->nodesAsGroups(), m, source, target, K);

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
		int pathIndex = paths.size();

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
			if( !p.first.contains("NOTHING") && !p.second.contains("NOTHING") )
			{
				for(auto sid : p.first) path.scolors[sid] = colors[i];
				for(auto tid : p.second) path.tcolors[tid] = colors[i];

				i++;
			}
		}

		path.gcorr->isReady = true;
		path.gcorr->correspondAllNodes();
		
		path.i = pathIndex;
	}

	// Timing
	property["prepareTime"].setValue( (int)prepareTimer.elapsed() );
	computeTimer.start();

	// Stats
	property["pathsCount"].setValue( (int)paths.size() );

	// Subsample paths
	if( isSubsample )
	{
		int MaxNumPaths = 2;
		std::vector<bool> mask = subsampleMask(MaxNumPaths, paths.size());
		std::vector<DeformationPath> subsampled;
		for(auto & path : paths){
			if(mask[path.i]) 
				subsampled.push_back(path);
		}
		paths = subsampled;
	}

	/// Find best correspondence
	DeformationPath bestPath;
	{
		beginFastNURBS();

		// Evaluate deformations
		#pragma omp parallel for
		for(int pi = 0; pi < (int)paths.size(); pi++)
		{
			auto & path = paths[pi];
			path.i = pi;

			// Prepare blending
			path.scheduler = QSharedPointer<Scheduler>( new Scheduler );
			path.blender = QSharedPointer<TopoBlender>( new TopoBlender( path.gcorr, path.scheduler.data() ) );

			// Deform
			path.scheduler->timeStep = 0.1;
			path.scheduler->executeAll();

			// Collect error
			double error = 0;

			Structure::Graph sourceCopy( *source );
			
			for(auto g : path.scheduler->allGraphs)
			{
				g->moveBottomCenterToOrigin( true );

				compute_part_measures( g );

				for(auto n_orig : sourceCopy.nodes)
				{
					Structure::Node * n = NULL;
					for(auto ni : g->nodes){
						if(ni->property.contains("original_ID") && ni->property.value("original_ID") == n_orig->id){
							n = ni; break;
						}
					}
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

		source->setColorAll(Qt::lightGray);
		target->setColorAll(Qt::lightGray);

		for( auto p : bestPath.pairs )
		{
			QColor c = starlab::qRandomColor2();
			for( auto sid : p.first ) source->setColorFor(sid, c);
			for( auto tid : p.second ) target->setColorFor(tid, c);
		}
	}
}
