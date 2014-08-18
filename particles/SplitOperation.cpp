#include "SplitOperation.h"
#include "myglobals.h"

double solidity_threshold = 0.6;
int max_parts_threshold = 12;
int level_threshold = 3;

auto rcolors = rndColors2(200);
int rc = 0;
QVector<RenderObject::Base*> globalDebug;

SplitOperation::SplitOperation(ParticleMesh *s, const SegmentGraph &seg, int level, SplitOperation *parent) :
    s(s), seg(seg), level(level)
{
    // Compute solidity score
    hull = ConvexHull<Vector3>( s->particlesCorners( seg.vertices ), "FA Qt" );

    // Segment solidity (Volume / Convex Volume)
    solidity = hull.solidity( s->grid.unitlength );
}

void SplitOperation::split()
{
    // Limit depth of search
    if (level > level_threshold)
        return;

    // Check if it needs splinting
    if (solidity < solidity_threshold)
    {
        // Collect and map current particles
        FloatVecVec descs;
        QMap<size_t,size_t> pmap;

        for(auto v : seg.vertices){
            pmap[v] = pmap.size();
            descs.push_back( s->desc[v] );
        }

        // Perform binary split
		int K = 2;
        clustering::kmeans<FloatVecVec,DestFn> km( descs, K );
        km._centers.clear();
        for(auto pid : s->specialSeeding(ParticleMesh::DESCRIPTOR, K, seg.vertices)) km._centers.push_back( s->desc[pid] );
        km.run();

        // Assign found clusters
        for(auto v : seg.vertices)
            s->particles[v].segment = km.cluster( pmap[v] );

        auto newSegments = s->segmentToComponents( seg, neiGraph );

        // Create possible split operations
        std::vector<SplitOperation> newOps;
        for(auto & newSeg : newSegments)
            newOps.push_back( SplitOperation(s, newSeg, level + 1, this) );

        // Analyze clusters
        for(auto & newOp : newOps)
        {
            // Too small of a segment
            if(newOp.seg.vertices.size() < 10){
                undecided.push_back(newOp);
                continue;
            }

            // No further split should happen
            if(newOp.seg.vertices.size() == seg.vertices.size()){
                children.clear();
                break;
            }

            // Add it
            children.push_back(newOp);

            // Further split?
            if(newOp.solidity < solidity_threshold)
            {
                children.back().split();
                children.back().seg.pid = seg.uid;
            }
        }
    }
}

void SplitOperation::collectClusters(std::vector<SplitOperation *> &ops)
{
    if(children.empty()){
        ops.push_back(this);
        return;
    }

    for(auto & child : children){
        child.collectClusters(ops);
        child.parent = this;
    }
}

void SplitOperation::debugAllChildren(QVector<RenderObject::Base *> &debug)
{
    if(children.empty())
    {
        starlab::PolygonSoup * ps = new starlab::PolygonSoup;
        for(auto f : hull.faces)
        {
            QVector<starlab::QVector3> pnts;
            for(auto v : f) pnts << v;
            QColor c = rcolors[rc];
            //c.setAlphaF(0.5);
            ps->addPoly(pnts, c);
        }
        rc++;
        debug << ps;
    }

    for(auto & child : children)
        child.debugAllChildren(debug);
}

void SplitOperation::report(QStringList &items)
{
    QString r = QString("[%1 - level %2] solidity (%3) parent(%4) size(%5)").arg(
                seg.uid).arg(level).arg(solidity).arg(seg.pid).arg(seg.vertices.size());

    if(!children.empty())
        r += QString(" children(%1)").arg(children.size());

    items << r;

    for(auto & child : children)
        child.report( items );
}
