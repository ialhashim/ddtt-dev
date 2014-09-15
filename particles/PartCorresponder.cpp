#include "PartCorresponder.h"

#pragma warning( disable : 4005 4100 4267 4616 4700 4291 4189 4267 4996 4267 ) // its not my code..
#include "flann/flann.hpp"
using namespace flann;
typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> Matrixf;
typedef Index< L2<float> > FlannIndex;

PartCorresponder::PartCorresponder(ParticleMesh *pmeshA, SegmentGraph segA,
                                   ParticleMesh *pmeshB, SegmentGraph segB)
    : sA(pmeshA), sB(pmeshB), segA(segA), segB(segB)
{
	// Enable debugging mode
	//sA->property["debug"].setValue( true );

	QVector<ParticleMesh*> input;
	input << sA << sB;

	QVector<SegmentGraph*> segment;
	segment << &segA << &segB;

	QVector<Eigen::AlignedBox3d> bbox( input.size() );
	std::vector< std::vector<SegmentGraph> > selectedEdges( input.size() );

	// Searching
	flann::SearchParams params;
	params.cores = omp_get_num_procs();
	QVector<Matrixf*> matrices;
	QVector<FlannIndex*> trees;

	/// For each input shape:
	for(size_t s = 0; s < input.size(); s++)
	{
		for(auto v : segment[s]->vertices)
		{
			auto & p = input[s]->particles[v];
			p.flag = VIZ_WEIGHT;
			p.weight = 0.0;

			bbox[s].extend( p.pos );
		}

		auto edges = input[s]->getEdgeParticlesOfSegment( *segment[s] );

		/// Compute relative edge heights
		QMap<int, double> edgeHeight;
		QMap<int, Vector3> edgeCenters;

		for(size_t i = 0; i < edges.size(); i++)
		{
			auto & edge = edges[i];

			Vector3 center(0,0,0);
			for(auto v : edge.vertices) center += input[s]->particles[v].pos;
			center /= edge.vertices.size();

			edgeHeight[i] = (center.z() - bbox[s].min().z()) / bbox[s].sizes().z();
			edgeCenters[i] = center;
		}

		auto edgeHeights = edgeHeight.values();

		double bottomEdgeHeight = *std::min_element(edgeHeights.begin(), edgeHeights.end());

		/// Collect bottom most edges
		std::vector<SegmentGraph> curSelected;

		for(size_t i = 0; i < edges.size(); i++)
		{
			if( abs(edgeHeight[i] - bottomEdgeHeight) > 0.25 ) continue;
			if( edges[i].vertices.empty() ) continue;
			
			curSelected.push_back( edges[i] );
			//debug << starlab::PointSoup::drawPoint(edgeCenters[i], 20, Qt::magenta);
		}

		selectedEdges.push_back( curSelected );

		// 1) Compute distances from selected edge
		if( false )
		{
			SegmentGraph curGraph;
			std::map<size_t,size_t> vmap;
			for(auto v : segment[s]->vertices) vmap[v] = vmap.size();
			for(auto e : segment[s]->GetEdgesSet()) curGraph.AddEdge( vmap[e.index], vmap[e.target], 1 );

			auto start_vertex = curGraph.AddVertex( curGraph.vertices.size() );

			// Connect all selected edges to the start vertex
			for(auto & selected : curSelected)
				for(auto v : selected.vertices)
					curGraph.AddEdge( start_vertex, vmap[v], 0 );

			curGraph.DijkstraComputePaths( start_vertex );

			double mindist = DBL_MAX, maxdist = -DBL_MAX;

			for(auto v : segment[s]->vertices)
			{
				auto d = curGraph.min_distance[ vmap[v] ];

				input[s]->particles[v].weight = d;

				mindist = std::min(mindist, d);
				maxdist = std::max(maxdist, d);
			}

			// Add normalized distance to bottom edge
			for(auto v : segment[s]->vertices)
				input[s]->particles[v].data << ( (input[s]->particles[v].weight-mindist) / (maxdist-mindist) );
		}

		// 2) Compute distance to average geodesic center
		if( false )
		{
			auto agd = input[s]->agd( 100, *segment[s] );

			std::vector<size_t> vpool;
			for(auto v : segment[s]->vertices) vpool.push_back( v );

			// Add normalized distance to center of segment
			for(size_t i = 0; i < agd.size(); i++)
				input[s]->particles[vpool[i]].data << agd[i];
		}

		// 3) Compute relative spatial positions
		if( true )
		{
			for(auto v : segment[s]->vertices)
			{
				Vector3 relative = (input[s]->particles[v].pos - bbox[s].min()).array() / bbox[s].sizes().array();

				// Add relative position
				input[s]->particles[v].data << relative.x() << relative.y() << relative.z();

				// Custom
				//input[s]->particles[v].data << input[s]->particles[v].pos.x() << input[s]->particles[v].pos.y();
			}
		}

		// Setup searching
		{
			auto v = *segment[s]->vertices.begin();
			int descSize = input[s]->particles[v].data.size();

			// Build matrix of part descriptors
			auto m = new Matrixf( segment[s]->vertices.size(), descSize );
			int r = 0;
			for(auto v : segment[s]->vertices)
				m->row(r++) = Eigen::Map<Eigen::VectorXf>(&input[s]->particles[v].data[0], descSize);

			matrices << m;

			// Build index tree	
			auto index = new FlannIndex(Matrix<float>( m->data(), m->rows(), m->cols() ), flann::KDTreeIndexParams());
			index->buildIndex();

			trees << index;
		}
	}

	// Find dense correspondence
	int knn = 1;

	for(size_t s = 0; s < input.size(); s++)
	{
		auto si = s, sj = (si+1) % input.size();

		std::vector< std::vector<int> > indices;
		std::vector< std::vector<float> > dists;

		auto m = matrices[ si ];
		trees[ sj ]->knnSearch( Matrix<float>( m->data(), m->rows(), m->cols() ), indices, dists, knn, params );

		std::vector<size_t> vipool, vjpool;
		for(auto v : segment[si]->vertices) vipool.push_back( v );
		for(auto v : segment[sj]->vertices) vjpool.push_back( v );

		for(size_t i = 0; i < vipool.size(); i++)
		{
			auto v = vipool[i];
			auto v_correspond = vjpool[indices[i].front()];

			input[s]->particles[v].correspondence = v_correspond;
		}
	}

	// Clean up
	for(size_t s = 0; s < input.size(); s++)
	{
		delete matrices[s];
		delete trees[s];
	}
}
