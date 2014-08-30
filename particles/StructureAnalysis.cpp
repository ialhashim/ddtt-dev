#include <QMessageBox>
#include "StructureAnalysis.h"
#include "myglobals.h"
#include "Planes.h"
#include "Bounds.h"
#include "graph_helper.h"

#include "disjointset.h"

#include "SplitOperation.h"

// Declare property types
Q_DECLARE_METATYPE(Eigen::Vector3d)
Q_DECLARE_METATYPE(Eigen::Vector3f)
Q_DECLARE_METATYPE(Eigen::AlignedBox3d)
Q_DECLARE_METATYPE(Boundsd)

StructureAnalysis::StructureAnalysis(ParticleMesh * pmesh) : s(pmesh)
{
	std::random_shuffle(ParticleMesh::rndcolors.begin(), ParticleMesh::rndcolors.end());

	rc = 0; // reset random color
	globalDebug.clear();

	SplitOperation op( pmesh, s->toGraph() );
	op.split();
	
	std::vector<SplitOperation*> clusters;
	op.collectClusters( clusters );

	std::map<size_t,size_t> mappedClusters;
	for(auto c : clusters)
		mappedClusters[c->seg.uid] = mappedClusters.size();

	// Assign new clusters
	for(auto c : clusters)
		for(auto v : c->seg.vertices)
			s->particles[v].segment = mappedClusters[c->seg.uid];

	bool isDone = false;

	while( !isDone )
	{
		isDone = true;

		// Get candidate good segments:
		SegmentGraph neiGraph;
		auto candidates = s->segmentToComponents( s->toGraph(), neiGraph );

		// Pre-compute convex hulls:
		QMap< size_t, ConvexHull<Vector3> > hulls;
		for(auto & seg : candidates)
			hulls[seg.uid] = ConvexHull<Vector3>( s->particlesCorners(seg.vertices) );

		// Merge smaller clusters
		bool isMerge = true;
		while( isMerge )
		{
			isMerge = false;

			// Sort candidates by their size
			std::vector<SegmentGraph*> sorted;
			for(auto c : candidates.keys()) if(candidates[c].vertices.size()) sorted.push_back( &candidates[c] );
			std::sort( sorted.begin(), sorted.end(), [&](const SegmentGraph* a, const SegmentGraph* b){ 
				return a->vertices.size() < b->vertices.size(); 
			});

			// Test if small cluster merges to a better solidity score
			for(auto seg : sorted)
			{
				auto & hull = hulls[seg->uid];
				std::map<int, ConvexHull<Vector3> > newHulls;
				int bestJ = seg->uid;
				double bestScore = solidity_threshold;

				double mysol = hull.solidity(s->grid.unitlength);
				if(mysol > 0.9) bestScore = mysol;

				for(auto j : neiGraph.getEdges(seg->uid))
				{
					if(candidates[j].vertices.empty()) continue;

					auto newHull = hull.merged(hulls[j]);
					double solidity = newHull.solidity(s->grid.unitlength);
				
					if(solidity > bestScore){
						bestScore = solidity;
						bestJ = j;
						newHulls[j] = newHull;
					}
				}

				if(bestJ == seg->uid) continue; // no merge is good

				auto big = &candidates[seg->uid];
				auto smaller = &candidates[bestJ];
				if(big->vertices.size() < smaller->vertices.size()) std::swap(big,smaller);

				// Migrate from small to big
				for(auto v : smaller->vertices) 
				{
					big->AddVertex(v);
					s->particles[v].segment = big->sid;
				}
			
				smaller->vertices.clear();

				hulls[big->uid] = newHulls[bestJ];

				isMerge = true;
				isDone = false;

				break;
			}
		}

		// Assign final segment IDs
		int sid = 0;
		for(auto & seg : candidates)
		{
			for(auto v : seg.vertices)
				s->particles[v].segment = sid;
			sid++;
		}
	}

	// Merge similar segments
	if( true )
	{
		SegmentGraph neiGraph;
		auto candidates = s->segmentToComponents( s->toGraph(), neiGraph );
		QMap<size_t,size_t> segMap, mapSeg;
		QMap<size_t,Vector3> segDirection;

		for(auto & seg : candidates)
		{
			mapSeg[segMap.size()] = seg.uid;
			segMap[seg.uid] = segMap.size();

			Eigen::AlignedBox3d bbox;
			for(auto v : seg.vertices) bbox.extend( s->particles[v].pos );

			QMap<double,size_t> dists;
			for(auto v : seg.vertices) dists[(s->particles[v].pos-bbox.center()).norm()] = v;
			auto rep = dists.values().first();
			//auto dir = s->mainDirection(rep);

			std::vector<Vector3> pnts;
			for(auto v : seg.vertices) pnts.push_back(s->particles[v].pos);
			auto mat = toEigenMatrix<double>(pnts);

			// PCA
			Eigen::MatrixXd centered = mat.rowwise() - mat.colwise().mean();
			Eigen::MatrixXd cov = centered.adjoint() * centered;
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);

			// Largest
			auto dir = eig.eigenvectors().col(2);

			// DEBUG:
			//auto vs = new starlab::VectorSoup;
			//vs->addVector(s->particles[rep].pos, Vector3(dir * 0.2));
			//debug << vs;

			segDirection[seg.uid] = dir;
		}

		DisjointSet disjoint(candidates.size());

		double similarity_threshold = 0.9;

		for(auto segID : neiGraph.vertices){
			for(auto neiID : neiGraph.GetNeighbours(segID)){
				double similarity = abs( segDirection[segID].dot(segDirection[neiID]) );
				if(similarity > similarity_threshold)
					disjoint.Union(segMap[segID], segMap[neiID]);
			}
		}

		for(size_t i = 0; i < candidates.size(); i++)
		{
			auto segID = mapSeg[i];
			for(auto v: candidates[segID].vertices)
				s->particles[v].segment = disjoint.Parent[i];
		}
	}

	// Draw hull around segments
	if( s->property["showHulls"].toBool() )
	{		
		SegmentGraph neiGraph;
		auto candidates = s->segmentToComponents( s->toGraph(), neiGraph );

		int sid = 0;
		
		for(auto & seg : candidates)
		{
			starlab::PolygonSoup * ps = new starlab::PolygonSoup;
			for( auto f : ConvexHull<Vector3>( s->particlesCorners(seg.vertices) ).faces )
			{
				QVector<starlab::QVector3> face;
				for(auto v : f) face << v;
				ps->addPoly(face, ParticleMesh::rndcolors[sid]);
			}

			sid++;
			debug << ps;
		}
	}

	// Visualize distance
	/*for(auto & seg : candidates)
	{
		std::vector<unsigned int> verts;
		for(auto v : seg.vertices) verts.push_back(v);

		std::vector<double> measures;
		for(auto v : verts) measures.push_back( s->particles[v].measure );

		auto bound = std::minmax_element( measures.begin(), measures.end() );

		int min_i = bound.first - measures.begin();
		int max_i = bound.second - measures.begin();
		double min_val =  measures[min_i];
		double max_val = measures[max_i];
		double range = max_val - min_val;

		for(auto v : seg.vertices)
		{
			double normalized = ((s->particles[v].measure-min_val) / range);
			s->particles[v].weight = pow(abs(normalized - 0.5) * 2,2);
			s->particles[v].flag = VIZ_WEIGHT;
		}
	}*/

	/*
	// Cluster separated segments:
	DisjointSet disjoint( candidates.size() );

	auto candidatesIDs = candidates.keys();
	QMap<size_t,size_t> centroids;

	// Closest particle to centroid
	for(size_t i = 0; i < candidatesIDs.size(); i++){
		uint centroid_id;
		double minDist = DBL_MAX;
		for(auto v : candidates[candidatesIDs[i]].vertices){
			double dist = (s->particles[v].pos - hulls[candidatesIDs[i]].center).norm();
			if(dist < minDist){
				minDist = dist;
				centroid_id = v;
			}
		}

		centroids[i] = centroid_id;
	}

	for(size_t i = 0; i < candidatesIDs.size(); i++)
	{
		for(size_t j = i + 1; j < candidatesIDs.size(); j++)
		{
			double threshold = 0.01;

			auto pi = centroids[i];
			auto pj = centroids[j];
			//auto di = Eigen::Map<Eigen::VectorXf>(&s->desc[pi][0], s->desc[pi].size());
			//auto dj = Eigen::Map<Eigen::VectorXf>(&s->desc[pj][0], s->desc[pj].size());
			//auto dist = (di.normalized() - dj.normalized()).norm();
			//auto dist = abs(s->particles[pi].measure - s->particles[pj].measure);
			
			auto dist = abs(s->particles[pi].flat - s->particles[pj].flat);

			if(dist < threshold) disjoint.Union(i,j);
		}
	}

	for(size_t i = 0; i < candidatesIDs.size(); i++)
	{
		for( auto v : candidates[candidatesIDs[i]].vertices )
		{
			s->particles[v].segment = disjoint.Parent[i];
		}
	}*/

	/*std::set<size_t> seeds;
	for(auto c : clusters)
	{
		// Closest particle to this centroid
		uint centroid_id;
		double minDist = DBL_MAX;
		for(auto v : c->seg.vertices){
			double dist = (s->particles[v].pos - c->hull.center).norm();
			if(dist < minDist){
				minDist = dist;
				centroid_id = v;
			}
		}

		if(seeds.size() < 3) seeds.insert( centroid_id );
	}
	s->cluster(seeds.size(), seeds, false, false);*/

	/*
	op.debugAllChildren(debug);
	QStringList rep; op.report(rep);
	saveToTextFile( "_split_report.txt", rep );
	*/

	for(auto d : globalDebug) debug << d;
}
