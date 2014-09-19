#pragma once
#include <QSharedPointer>
#include "ParticleMesh.h"

#include <NanoKdTree.h>

extern int slice_uid;
struct SliceChunk{
	Eigen::AlignedBox3d box;
	SegmentGraph g;
	QSharedPointer<NanoKdTree> tree;
	std::vector<size_t> vmap;
	SliceChunk(const SegmentGraph & graph) : g(graph), uid(slice_uid++) {}
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

Q_DECLARE_METATYPE( Segments );
Q_DECLARE_METATYPE( Slices );

class PartCorresponder
{
public:
	static Slices computeSlices( ParticleMesh * input, const SegmentGraph & seg );

	static void correspondSegments( const QPair<size_t,size_t> & segmentsPair, 
		const QVector<ParticleMesh *> & input, QVector<Particles> & particles );

	static QVector< QPair<int,int> > distributeVectors(int x, int y);
};
