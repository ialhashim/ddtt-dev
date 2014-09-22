#include "PartCorresponder.h"
#include "Bounds.h"

int slice_uid = 0;

Slices PartCorresponder::computeSlices( ParticleMesh * input, const SegmentGraph & seg )
{
	static int MAX_SLICE_COUNT = 20;

	Slices curSlices;
	if(seg.vertices.empty()) return curSlices;

	// Compute range of 'z' grid values
	Bounds<int> zcoords;
	for(auto v : seg.vertices)
	{
		Vector3 pg = (input->particles[v].pos - input->grid.translation.cast<double>()) / input->grid.unitlength;
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

		for(auto & edge : seg.GetEdgesSet())
		{
			Vector3 pgi = (input->particles[edge.index].pos - input->grid.translation.cast<double>()) / input->grid.unitlength;
			Vector3 pgj = (input->particles[edge.target].pos - input->grid.translation.cast<double>()) / input->grid.unitlength;

			int zi = (int)pgi.z(), zj = (int)pgj.z();

			if((zi < bottom || zj < bottom) || (zi >= top || zj >= top)) continue;

			if(zi >= bottom && zi < top) layerGraph.AddVertex(edge.index);
			if(zj >= bottom && zj < top) layerGraph.AddVertex(edge.target);

			layerGraph.AddEdge( edge.index, edge.target, 1 );
		}

		auto & slice = curSlices[i];
		slice.chunksFromGraphs( layerGraph.toConnectedParts() );

		// Compute chunks parameters
		for(size_t c = 0; c < slice.chunks.size(); c++)
		{
			auto & chunk = slice.chunks[c];

			// Compute bounding box enclosing chunk
			for(auto p : input->particlesCorners(chunk.g.vertices)) 
				chunk.box.extend(p);

			// Compute index of relative particle positions inside box
			chunk.tree = QSharedPointer<NanoKdTree>(new NanoKdTree);
			for(auto v : chunk.g.vertices)
			{
				input->particles[v].relativePos = (input->particles[v].pos - chunk.box.min()).array() / chunk.box.sizes().array();
				chunk.tree->addPoint( input->particles[v].relativePos );
				chunk.vmap.push_back( v );
			}
			chunk.tree->build();
		}
	}

	return curSlices;
}

void PartCorresponder::correspondSegments(const QPair<size_t,size_t> & segmentsPair, 
										  const QVector<ParticleMesh *> & input, 
										  QVector<Particles> & particles)
{
	QVector<Slices> slices;
	slices << input[0]->property["segments"].value<Segments>()[segmentsPair.first].property["slices"].value<Slices>();
	slices << input[1]->property["segments"].value<Segments>()[segmentsPair.second].property["slices"].value<Slices>();

	auto ii = 0, ij = 1;

	// i has more slices than j
	bool isSwap = slices[ii].size() < slices[ij].size();

	/// Correspond chunks:
	{
		if( isSwap ) std::swap(ii, ij);

		for(size_t sliceID = 0; sliceID < slices[ii].size(); sliceID++)
		{
			auto & slice_i = slices[ii][sliceID];
			double a = double(sliceID) / std::max(1, (slices[ii].size()-1));
			int j_idx =  std::ceil(a * (slices[ij].size()-1));
			auto & slice_j = slices[ij][j_idx];

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

			if(sortedChunks.front().empty() || sortedChunks.back().empty()) continue;

			// Now divide and assign
			auto assignments = PartCorresponder::distributeVectors( sortedChunks.front().size(), sortedChunks.back().size() );

			/// Match chunks:
			for(auto assignment : assignments)
			{
				auto & chunk_i = slice_i.chunks[assignment.first];
				auto & chunk_j = slice_j.chunks[assignment.second];

				/// Correspond chunks:
				{
					QVector<SliceChunk*> chunk; chunk << &chunk_i << &chunk_j;

					if( isSwap ) std::reverse(chunk.begin(), chunk.end());

					for(size_t si = 0; si < input.size(); si++)
					{
						auto sj = (si+1) % input.size();

						std::vector<size_t> closestMap( chunk[si]->vmap.size(), -1 );

						// Find closest particle
						for(int i = 0; i < (int)chunk[si]->vmap.size(); i++)
						{
							auto vi = chunk[si]->vmap[i];
							auto vj = chunk[sj]->vmap[chunk[sj]->tree->closest( particles[si][vi].relativePos )];

							closestMap[i] = vj;
						}

						// Match particles one-to-one
						for(int i = 0; i < (int)chunk[si]->vmap.size(); i++)
						{
							auto vi = chunk[si]->vmap[i];
							auto vj = closestMap[i];

							auto pi = &particles[si][vi];
							auto pj = &particles[sj][vj];

							if( pi->isMatched )
							{
								particles[si].push_back( *pi );
								pi = &particles[si].back();
								pi->id = particles[si].size()-1;
							}

							if( pj->isMatched )
							{
								particles[sj].push_back( *pj );
								pj = &particles[sj].back();
								pj->id = particles[sj].size()-1;
							}

							pi->correspondence = pj->id;
							pj->correspondence = pi->id;

							pj->isMatched = true;
							pi->isMatched = true;
						}
					}
				}
			}
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
