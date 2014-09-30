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
	SliceChunk(const SegmentGraph & graph = SegmentGraph()) : g(graph), uid(slice_uid++) {}
	int uid;
	QMap<QString,QVariant> property;
};

struct Slice{
	QVector<SliceChunk> chunks;
	void chunksFromGraphs(const std::vector<SegmentGraph> & graphs){
		for(auto & graph : graphs)
			chunks.push_back( SliceChunk(graph) );
	}
	Slice() : isGrid2D(false){}
	bool isGrid2D;
	Eigen::AlignedBox3d box(){ Eigen::AlignedBox3d b; for(auto & c : chunks) b.extend(c.box); return b; }
};
typedef QVector<Slice> Slices;

Q_DECLARE_METATYPE( Segments );
Q_DECLARE_METATYPE( Slices );
Q_DECLARE_METATYPE( Eigen::AlignedBox3d );

class PartCorresponder
{
public:
	static Slices computeSlices( ParticleMesh * input, const SegmentGraph & seg );

	static void correspondSegments( const QPair<size_t,size_t> & segmentsPair, 
		const QVector<ParticleMesh *> & input, QVector<Particles> & particles );

	static void matchChunk( QVector<SliceChunk> chunk, const QVector<ParticleMesh *> & input, 
		QVector<Particles> & particles );

	static void matchGridChunk(const QVector< QVector<SliceChunk> > & chunk, 
		bool isGrid2D, const QVector<ParticleMesh *> & input, QVector<Particles> & particles);

	static void match1DGridChunk( QVector< QVector<SliceChunk> > chunk, const QVector<ParticleMesh *> & input, 
		QVector<Particles> & particles );

	// helpers:
	static SliceChunk mergeChunks( const QVector<SliceChunk> & chunks, ParticleMesh * input, Particles & particles );
	static QVector< QPair<int,int> > distributeVectors(int x, int y);
};
