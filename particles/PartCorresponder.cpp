#include "PartCorresponder.h"
#include "Bounds.h"

int slice_uid = 0;

PartCorresponder::PartCorresponder(ParticleMesh *pmeshA, SegmentGraph segA,
                                   ParticleMesh *pmeshB, SegmentGraph segB)
    : sA(pmeshA), sB(pmeshB), segA(segA), segB(segB)
{
	// Enable debugging mode
	//sA->property["debug"].setValue( true );

	QVector<ParticleMesh*> input; input << sA << sB;
	QVector<SegmentGraph*> segment; segment << &segA << &segB;
	QVector<Slices> slices( input.size() ) ;

	int MAX_SLICE_COUNT = 15;

	/// Compute slices for each input shape:
	for(size_t s = 0; s < input.size(); s++)
	{
		auto & curSlices = slices[s];

		// Compute range of 'z' grid values
		Bounds<int> zcoords;
		for(auto v : segment[s]->vertices)
		{
			Vector3 pg = (input[s]->particles[v].pos - input[s]->grid.translation.cast<double>()) / input[s]->grid.unitlength;
			zcoords.extend( (int)pg.z() );
		}

		int total_range = 1 + zcoords.range();
		int fixedNumLayers = std::min(total_range, MAX_SLICE_COUNT);
		int perLayer = std::floor((double)total_range / fixedNumLayers);
		int numLayers = std::ceil((double)total_range / perLayer);

		curSlices.resize( numLayers );

		// Divide segment into layers
		for(int i = 0; i < numLayers; i++)
		{
			int bottom = zcoords.minimum + (i * perLayer);
			int top = bottom + perLayer;

			SegmentGraph layerGraph;

			for(auto edge : segment[s]->GetEdgesSet())
			{
				Vector3 pgi = (input[s]->particles[edge.index].pos - input[s]->grid.translation.cast<double>()) / input[s]->grid.unitlength;
				Vector3 pgj = (input[s]->particles[edge.target].pos - input[s]->grid.translation.cast<double>()) / input[s]->grid.unitlength;

				int zi = (int)pgi.z(), zj = (int)pgj.z();

				if(zi >= bottom && zi < top) layerGraph.AddVertex(edge.index);
				if(zj >= bottom && zj < top) layerGraph.AddVertex(edge.target);

				if((zi < bottom || zj < bottom) || (zi >= top || zj >= top)) continue;

				layerGraph.AddEdge( edge.index, edge.target, 1 );
			}

			auto & slice = curSlices[i];
			slice.chunksFromGraphs( layerGraph.toConnectedParts() );

			// Compute chunks parameters
			for(size_t c = 0; c < slice.chunks.size(); c++)
			{
				auto & chunk = slice.chunks[c];

				// Compute bounding box enclosing chunk
				for(auto p : input[s]->particlesCorners(chunk.g.vertices)) 
					chunk.box.extend(p);

				// Compute index of relative particle positions inside box
				chunk.tree = new NanoKdTree;
				for(auto v : chunk.g.vertices)
				{
					input[s]->particles[v].relativePos = (input[s]->particles[v].pos - chunk.box.min()).array() / chunk.box.sizes().array();
					chunk.tree->addPoint( input[s]->particles[v].relativePos );
					chunk.vmap.push_back( v );
				}
				chunk.tree->build();

				//auto bs = new starlab::BoxSoup;bs->addBox( chunk.box );debug<<bs;
			}
		}
	}

	/// Correspond chunks:	
	std::vector< std::pair< int, std::pair<size_t,size_t> > > assignments;
	{
		auto si = 0, sj = 1;

		// i has more slices than j
		if( slices[si].size() < slices[sj].size() ) std::swap(si, sj);

		for(size_t sliceID = 0; sliceID < slices[si].size(); sliceID++)
		{
			auto & slice_i = slices[si][sliceID];
			double a = double(sliceID) / std::max(1, (slices[si].size()-1));
			auto & slice_j = slices[sj][a * (slices[sj].size()-1)];

			/// First project chunks onto diagonal of their slice
			QVector< QVector< size_t > > sortedChunks;

			QVector<Slice*> bothSlices;
			bothSlices << &slice_i << &slice_j;

			for( auto slice : bothSlices )
			{
				Eigen::AlignedBox3d slice_box;
				for(auto & chunk : slice->chunks) slice_box.extend( chunk.box );

				std::vector< std::pair<double,size_t> > projected;

				for(size_t idx = 0; idx < slice->chunks.size(); idx++)
				{
					Vector3 vec1 = slice->chunks[idx].box.center() - slice_box.min();
					Vector3 vec2 = (slice_box.diagonal()).normalized();
					double t = (vec1).dot(vec2);

					projected.push_back( std::make_pair(t, idx) );
				}

				std::sort(projected.begin(), projected.end());

				QVector< size_t > sortedIndices;
				for(auto proj : projected) sortedIndices << proj.second;
				sortedChunks << sortedIndices;
			}

			// Now divide and assign
			auto assignments = distributeVectors( sortedChunks.front().size(), sortedChunks.back().size() );

			/// Match chunks:
			for(auto assignment : assignments)
			{
				auto & chunk_i = slice_i.chunks[assignment.first];
				auto & chunk_j = slice_j.chunks[assignment.second];

				correspondChunks( chunk_i, input[si], chunk_j, input[sj] );
			}
		}
	}
}

void PartCorresponder::correspondChunks( SliceChunk & chunk_i, ParticleMesh* mesh_i, 
										SliceChunk & chunk_j, ParticleMesh* mesh_j )
{
	QVector<ParticleMesh*> input; input << mesh_i << mesh_j;
	QVector<SliceChunk*> chunk; chunk << &chunk_i << &chunk_j;

	for(size_t si = 0; si < input.size(); si++)
	{
		auto sj = (si+1) % input.size();

		std::vector<size_t> closestMap( chunk[si]->vmap.size(), -1 );

		// Find closest particle
		//#pragma omp parallel for
		for(int i = 0; i < (int)chunk[si]->vmap.size(); i++)
		{
			auto vi = chunk[si]->vmap[i];
			auto vj = chunk[sj]->vmap[chunk[sj]->tree->closest( input[si]->particles[vi].relativePos )];

			closestMap[i] = vj;
		}

		// Match particles one-to-one
		for(int i = 0; i < (int)chunk[si]->vmap.size(); i++)
		{
			auto vi = chunk[si]->vmap[i];
			auto vj = closestMap[i];

			auto pi = &input[si]->particles[vi];
			auto pj = &input[sj]->particles[vj];

			if( pi->isMatched )
			{
				input[si]->particles.push_back( *pi );
				pi = &input[si]->particles.back();
				pi->id = input[si]->particles.size()-1;
			}

			if( pj->isMatched )
			{
				input[sj]->particles.push_back( *pj );
				pj = &input[sj]->particles.back();
				pj->id = input[sj]->particles.size()-1;
			}

			pi->correspondence = pj->id;
			pj->correspondence = pi->id;

			pj->isMatched = true;
			pi->isMatched = true;
		}
	}
}

QVector< QPair<int,int> > PartCorresponder::distributeVectors(int x, int y)
{
	QSet< QPair<int,int> > result;

	if(x == y && x == 1){
		result << qMakePair(0,0);
		return result.toList().toVector();
	}

	int start_i = 0, start_j = 0;
	int end_i = x-1, end_j = y-1;

	bool isStopI = false, isStopJ = false;

	while( true ){
		result << qMakePair(start_i,start_j) << qMakePair(end_i, end_j);

		if(end_i - start_i < 2) { isStopI = true; }
		if(end_j - start_j < 2) { isStopJ = true; }
		if( isStopI && isStopJ ) break;

		if( !isStopI ){
			start_i++;
			end_i--;
		}

		if( !isStopJ ){
			start_j++;
			end_j--;
		}
	}

	return result.toList().toVector();
}
