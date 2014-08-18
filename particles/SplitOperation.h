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
    std::vector<SplitOperation> children, undecided;
    SplitOperation * parent;
    ConvexHull<Vector3> hull;
    ParticleMesh * s;

    typedef std::vector<float> FloatVec;
    typedef std::vector<FloatVec> FloatVecVec;
    typedef clustering::lpnorm< VectorFloat > DestFn;

    SplitOperation( ParticleMesh * s, const SegmentGraph& seg, int level = 0, SplitOperation* parent = NULL);

    void split();

    void collectClusters( std::vector<SplitOperation*> & ops );

    void debugAllChildren(QVector<RenderObject::Base*> & debug);

    void report( QStringList & items );
};
