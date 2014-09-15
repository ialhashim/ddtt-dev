#include "PartCorresponder.h"
#include "Bounds.h"

#include <NanoKdTree.h>

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

	std::vector< std::vector< std::vector<SegmentGraph> > > layer( input.size() );
	std::vector< std::vector< std::vector<NanoKdTree*> > > tree( input.size() );
	std::vector< std::vector< std::vector<Eigen::AlignedBox3d> > > lbox( input.size() );

	/// For each input shape:
	for(size_t s = 0; s < input.size(); s++)
	{
		Bounds<int> zcoords;
		std::map<size_t,bool> used;

		for(auto v : segment[s]->vertices)
		{
			auto & p = input[s]->particles[v];
			bbox[s].extend( p.pos );

			Vector3 pg = (p.pos - input[s]->grid.translation.cast<double>()) / input[s]->grid.unitlength;
			zcoords.extend( (int)pg.z() );
		}

		int total_range = 1 + zcoords.range();
		int fixedNumLayers = std::min(total_range, 10);
		int perLayer = std::floor((double)total_range / fixedNumLayers);
		int numLayers = std::ceil((double)total_range / perLayer);

		std::vector< std::vector<SegmentGraph> > layers( numLayers );
		std::vector< std::vector<NanoKdTree*> > trees( numLayers );
		std::vector< std::vector<Eigen::AlignedBox3d> > boxes( numLayers );

		auto bs = new starlab::BoxSoup;
		QColor color = starlab::qRandomColor2();

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

			layers[i] = layerGraph.toConnectedParts();

			std::vector<Eigen::AlignedBox3d> layerBoxes;

			for(auto g : layers[i])
			{
				Eigen::AlignedBox3d box;
				for(auto p : input[s]->particlesCorners(g.vertices)) box.extend(p);
				layerBoxes.push_back( box );

				auto t = new NanoKdTree;
				for(auto v : g.vertices) 
				{
					Vector3 p = (input[s]->particles[v].pos - box.min()).array() / box.sizes().array();
					t->addPoint( p );
				}
				t->build();
				trees[i].push_back(t);

				// DEBUG:
				//bs->addBox(box, color);
			}

			boxes[i] = layerBoxes;
		}

		debug << bs;

		layer[s] = layers;
		tree[s] = trees;
		lbox[s] = boxes;
	}

	for(size_t s = 0; s < input.size(); s++)
	{
		auto si = s, sj = (si+1) % input.size();

		for(size_t l = 0; l < layer[si].size(); l++)
		{
			double t = double(l) / std::max( size_t(1), (layer[si].size()-1) );

			int idxj = t * (layer[sj].size()-1);
			auto & targetLayer = layer[sj][idxj].front();
			auto & targetTree = tree[sj][idxj].front();

			std::vector<size_t> vjmap;
			for(auto v : targetLayer.vertices) vjmap.push_back( v );

			for(size_t k = 0; k < layer[si][l].size(); k++)
			{
				auto & curLayer = layer[si][l][k];
				auto & curBox = lbox[si][l][k];

				for( auto vi : curLayer.vertices )
				{
					Vector3 p = (input[si]->particles[vi].pos - curBox.min()).array() / curBox.sizes().array();
					auto vj = vjmap[targetTree->closest( p )];

					input[si]->particles[vi].correspondence = vj;
				}
			}
		}
	}
}
