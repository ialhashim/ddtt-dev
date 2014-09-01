#pragma once

#include "ParticleMesh.h"

#include "convexhull.h"
#include "kmeans.h"

extern double solidity_threshold;

extern int rc;
extern QVector<RenderObject::Base*> globalDebug;

struct SplitOperation{
    SegmentGraph seg, neiGraph;
    double solidity;
    int level;
    std::vector<SplitOperation> children;
    ParticleMesh * s;
	SplitOperation * parent;

	ConvexHull<Vector3> hull;
	Eigen::AlignedBox3d bbox();

    typedef std::vector<float> FloatVec;
    typedef std::vector<FloatVec> FloatVecVec;
    typedef clustering::l2norm_squared<FloatVec> DestFn;

    SplitOperation( ParticleMesh * s, const SegmentGraph& seg, int level = 0);

    void split();

    void collectClusters( std::vector<SplitOperation*> & ops );

    void debugAllChildren(QVector<RenderObject::Base*> & debug);

    void report( QStringList & items );
};
