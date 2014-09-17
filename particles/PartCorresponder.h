#pragma once
#include "ParticleMesh.h"

#include <NanoKdTree.h>
extern int slice_uid;
struct SliceChunk{
	Eigen::AlignedBox3d box;
	SegmentGraph g;
	NanoKdTree * tree;
	std::vector<size_t> vmap;
	SliceChunk(const SegmentGraph & graph) : g(graph), tree(NULL), uid(slice_uid++) {}
	~SliceChunk(){ if(tree) delete tree; }
	int uid;
	QMap<QString,QVariant> property;
};

struct Slice{
	std::vector<SliceChunk> chunks;
	void chunksFromGraphs(const std::vector<SegmentGraph> & graphs){
		for(auto & graph : graphs)
			chunks.push_back( SliceChunk(graph) );
	}
};

typedef QVector<Slice> Slices;

class PartCorresponder
{
public:
    PartCorresponder( ParticleMesh * pmeshA, SegmentGraph segA,
                      ParticleMesh * pmeshB, SegmentGraph segB );

    ParticleMesh *sA, *sB;
    SegmentGraph segA, segB;

    QVector<RenderObject::Base*> debug;

	QVector< QPair<int,int> > distributeVectors(int x, int y);
	void correspondChunks( SliceChunk & chunk_i, ParticleMesh* mesh_i, SliceChunk & chunk_j, ParticleMesh* mesh_j );
};
