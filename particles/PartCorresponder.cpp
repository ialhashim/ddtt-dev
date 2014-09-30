#include "PartCorresponder.h"
#include "Bounds.h"
#include "myglobals.h"

int slice_uid = 0;

QVector<Eigen::AlignedBox3d> splitBox( const Eigen::AlignedBox3d & box, Vector3 axis )
{
	Eigen::AlignedBox3d b1 = box, b2 = box;
	Vector3 delta = axis.array() * ((box.sizes() * 0.5).array());
	return QVector<Eigen::AlignedBox3d>() << box.intersection(b1.translate(-delta)) <<
		box.intersection(b2.translate(delta));
}

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

		if( layerGraph.vertices.empty() ) continue;

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

	Slices goodSlices;
	for(auto slice : curSlices) 
		if(!slice.chunks.empty()) 
			goodSlices.push_back(slice);
	return goodSlices;
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
			if( isSwap ) std::reverse(bothSlices.begin(), bothSlices.end());

			// Most basic case, a single chunk matches with another
			if( slice_i.chunks.size() == 1 && slice_j.chunks.size() == 1 )
			{
				QVector<SliceChunk> chunks;
				chunks << slice_i.chunks.front() << slice_j.chunks.front();
				matchChunk( chunks, input, particles );
				continue;
			}

			// Find out if we have chunks on a 2D grid in either slice
			bool isGrid2D = false;
			double grid_spread_threshold = 0.5;

			for(size_t shapei = 0; shapei < input.size(); shapei++)
			{
				auto & slice = bothSlices[shapei];
				if(slice->chunks.size() < 3) continue;

				Eigen::MatrixXd points( slice->chunks.size(), 3 ); size_t r = 0;
				for(auto & chunk : slice->chunks) points.row(r++) = chunk.box.center();
				Vector3 b_center (points.colwise().mean());
				points = points.rowwise() - b_center.transpose();
				Eigen::JacobiSVD<Eigen::MatrixXd> svd(points, Eigen::ComputeThinU | Eigen::ComputeThinV);
				Eigen::Vector3d values = svd.singularValues().normalized();
				double ratio = (values[1] / values[0]);

				if(ratio > grid_spread_threshold)
				{
					isGrid2D = true;
					slice->isGrid2D = true;
				}
			}

			// Split 1D grids along the 'y' axis when matching against 2D grids
			if( isGrid2D )
			{
				for(size_t i = 0; i < bothSlices.size(); i++)
				{
					auto & slice = bothSlices[i];
					if(slice->isGrid2D) continue; // no splits needed

					QVector<SliceChunk> splitChunks;
					
					for(auto & oldChunk : slice->chunks)
					{
						auto halves = splitBox( oldChunk.box, Vector3::UnitY() );

						for(size_t h = 0; h < halves.size(); h++)
						{
							SliceChunk chunk;

							// Migrate vertex to its designated half
							for(auto v : oldChunk.g.vertices)
								if(halves[h].contains( input[i]->particles[v].pos ))
									chunk.g.AddVertex(v);

							// Compute bounding box enclosing chunk
							for(auto p : input[i]->particlesCorners(chunk.g.vertices)) 
								chunk.box.extend(p);

							// Compute index of relative particle positions inside box
							chunk.tree = QSharedPointer<NanoKdTree>(new NanoKdTree);
							for(auto v : chunk.g.vertices)
							{
								particles[i][v].relativePos = (particles[i][v].pos - chunk.box.min()).array() / chunk.box.sizes().array();
								chunk.tree->addPoint( particles[i][v].relativePos );
								chunk.vmap.push_back( v );
							}
							chunk.tree->build();

							splitChunks.push_back(chunk);
						}
					}

					// Replace slice's chunks with newly split ones
					slice->chunks = splitChunks;
				}
			}

			QVector< QVector<SliceChunk> > gridChunks;
			gridChunks << bothSlices.front()->chunks << bothSlices.back()->chunks;

			matchGridChunk(gridChunks, isGrid2D, input, particles);
		}
	}
}

void PartCorresponder::matchGridChunk(const QVector< QVector<SliceChunk> > & chunk, bool isGrid2D,
									  const QVector<ParticleMesh *> & input, QVector<Particles> & particles)
{
	if( !isGrid2D )
	{
		// Sort chunks along 1D grid
		QVector< QVector<SliceChunk> > sorted;
		for(auto gridChunk : chunk)
		{
			// Sort based on coordinates of chunk center
			std::sort(gridChunk.begin(), gridChunk.end(), [&]( const SliceChunk& a, const SliceChunk& b ){ 
				double ma = (a.box.center().x() * 1e6) + a.box.center().y();
				double mb = (b.box.center().x() * 1e6) + b.box.center().y();
				return ma < mb; 
			});

			sorted.push_back( gridChunk );
		}

		// match 1D grid chunks
		match1DGridChunk( sorted, input, particles );
	}
	else
	{		
		// Compute slice bounding box
		QVector< Eigen::AlignedBox3d > sliceChunk(2);
		for(size_t i = 0; i < input.size(); i++){
			auto & inputchunks = chunk[i];
			for(auto & chunk : inputchunks) 
				sliceChunk[i].extend(chunk.box);
		}

		// Split on 'y' into two rows of 1D grids
		QVector< QVector< QVector<SliceChunk> > > splitRows(input.size());

		for(size_t i = 0; i < input.size(); i++)
		{
			splitRows[i].resize(2);

			auto & inputchunks = chunk[i];
			for(size_t ci = 0; ci < inputchunks.size(); ci++)
			{
				if( inputchunks[ci].box.center().y() < sliceChunk[i].center().y() )
					splitRows[i].front().push_back( inputchunks[ci] );
				else
					splitRows[i].back().push_back( inputchunks[ci] );
			}
		}

		// For each row
		for(int r = 0; r < splitRows.front().size(); r++)
		{
			QVector< QVector<SliceChunk> > row;

			for(int i = 0; i < input.size(); i++)
				row.push_back( splitRows[i][r] );

			matchGridChunk(row, false, input, particles);
		}
	}
}

void PartCorresponder::match1DGridChunk( QVector< QVector<SliceChunk> > sortedChunk, const QVector<ParticleMesh *> & input, QVector<Particles> & particles )
{
	QVector< QVector<SliceChunk> > readyChunkPairs;

	auto & sortedChunkFront = sortedChunk.front();
	auto & sortedChunkBack = sortedChunk.back();

	// Different number of chunks
	if( sortedChunkFront.size() != sortedChunkBack.size() )
	{
		int targetSize = std::min( sortedChunkFront.size(), sortedChunkBack.size() );

		if(targetSize == 1)
		{
			QVector<SliceChunk> chunkPairs;
			if( sortedChunkFront.size() == 1 ) chunkPairs.push_back( sortedChunkFront.front() );
			else chunkPairs.push_back( mergeChunks( sortedChunkFront, input.front(), particles.front() ) );

			if( sortedChunkBack.size() == 1) chunkPairs.push_back( sortedChunkBack.front() );
			else chunkPairs.push_back( mergeChunks( sortedChunkBack, input.back(), particles.back() ) );

			readyChunkPairs.push_back( chunkPairs );
		}
		else
		{
			// For now we use basic matching.. later we should either split / merge
			for(auto v : distributeVectors(sortedChunkFront.size(), sortedChunkBack.size()))
			{
				QVector<SliceChunk> p;
				p << sortedChunkFront[v.first] << sortedChunkBack[v.second];
				readyChunkPairs.push_back( p );
			}
		}
	}
	else
	{
		// Same number of elements, simply match them up
		for(size_t i = 0; i < sortedChunk.front().size(); i++)
		{
			readyChunkPairs.push_back( QVector<SliceChunk>() << sortedChunkFront.at(i) << sortedChunkBack.at(i) );
		}
	}

	// Match each pair of chunks
	for(auto & pairChunk : readyChunkPairs)
		matchChunk(pairChunk, input, particles);
}

void PartCorresponder::matchChunk( QVector<SliceChunk> chunk, const QVector<ParticleMesh *> & input, QVector<Particles> & particles )
{
	for(size_t si = 0; si < input.size(); si++)
	{
		auto sj = (si+1) % input.size();

		if(chunk[si].vmap.empty() || chunk[sj].vmap.empty()) continue;

		std::vector<size_t> closestMap( chunk[si].vmap.size(), -1 );

		// Find closest particle
		for(int i = 0; i < (int)chunk[si].vmap.size(); i++)
		{
			auto vi = chunk[si].vmap[i];
			auto vj = chunk[sj].vmap[chunk[sj].tree->closest( particles[si][vi].relativePos )];

			closestMap[i] = vj;
		}

		// Match particles one-to-one
		for(int i = 0; i < (int)chunk[si].vmap.size(); i++)
		{
			auto vi = chunk[si].vmap[i];
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

SliceChunk PartCorresponder::mergeChunks(const QVector<SliceChunk> & chunks, ParticleMesh * input, Particles & particles)
{
	SliceChunk chunk;
	if(chunks.empty()) return chunk;
	if(chunks.size() == 1) return chunks.front();

	// Combine graphs - is it needed?
	if( false ){
		for(auto & c : chunks){
			for(auto v : c.g.vertices) chunk.g.AddVertex(v);
			for(auto e : c.g.GetEdgesSet()) chunk.g.AddEdge(e.index, e.target, 1);
		}
	}

	// Remove gaps between chunks
	QVector<Vector3> packedPoints;

	// Decide if we will allow moves along 'x' or 'y'
	Vector3 diagonal = chunks.back().box.center() - chunks.front().box.center();
	diagonal = diagonal.cwiseAbs().normalized();
	if(diagonal[0] > diagonal[1]) diagonal = Vector3(1,0,0);
	else diagonal = Vector3(0,1,0);

	Vector3 halfVoxel(input->grid.unitlength,input->grid.unitlength,input->grid.unitlength);
	halfVoxel *= 0.5;

	// Initial start
	Vector3 delta = chunks.front().box.min();

	// Snap chunks to each other
	for(size_t ci = 0; ci < chunks.size(); ci++)
	{
		for(auto v : chunks[ci].g.vertices)
		{
			auto p = particles[v].pos - delta;
			packedPoints.push_back( p );

			chunk.box.extend( p + halfVoxel );
			chunk.box.extend( p - halfVoxel );

			// Save mapped index
			chunk.vmap.push_back( v );
		}

		if(ci+1 == chunks.size()) continue;

		Vector3 d = chunk.box.diagonal().array() * diagonal.array();
		delta = chunks[ci+1].box.min() - d;
	}

	// Compute index of relative particle positions inside box
	chunk.tree = QSharedPointer<NanoKdTree>(new NanoKdTree);
	for(size_t i = 0; i < packedPoints.size(); i++)
	{
		Vector3 relativePos = (packedPoints[i] - chunk.box.min()).array() / chunk.box.sizes().array();
		chunk.tree->addPoint( relativePos );
		particles[chunk.vmap[i]].relativePos = relativePos;
	}
	chunk.tree->build();

	return chunk;
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
