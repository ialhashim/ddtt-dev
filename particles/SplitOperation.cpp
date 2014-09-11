#include "SplitOperation.h"
#include "myglobals.h"
#include "BasicTable.h"
#include <queue>

double solidity_threshold = 0.65;
double isect_ratio_threshold = 0.9;
int max_parts_threshold = 12;
int level_threshold = 3;

auto rcolors = rndColors2(200);
int rc = 0;
QVector<RenderObject::Base*> globalDebug;

SplitOperation::SplitOperation(ParticleMesh *s, const SegmentGraph &seg, int level) :
    s(s), seg(seg), level(level), parent(NULL)
{
    // Compute solidity score
    hull = ConvexHull<Vector3>( s->particlesCorners( seg.vertices ), "FA Qt" );

    // Segment solidity (Volume / Convex Volume)
    solidity = hull.solidity( s->grid.unitlength );
}

void SplitOperation::split()
{
	// Check if it needs splitting
	if (solidity > solidity_threshold)
		return;

    // Limit depth of search
    if (level > level_threshold)
        return;

	int small_segment_threshold = 0.25 * s->grid.gridsize;

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
	auto seeds = s->specialSeeding(ParticleMesh::DESCRIPTOR, K, seg.vertices);

    for(auto pid : seeds) km._centers.push_back( s->desc[pid] );
    km.run();

    // Assign found clusters
    for(auto v : seg.vertices)
        s->particles[v].segment = km.cluster( pmap[v] );

	// Segment based on clustering
    auto newSegments = s->segmentToComponents( seg, neiGraph );

    // Create possible split operations
    std::vector<SplitOperation> newOps;
    for(auto & newSeg : newSegments)
        newOps.push_back( SplitOperation(s, newSeg, level + 1) );

	std::vector<SplitOperation> undecided;

    // Analyze clusters
    for(auto & newOp : newOps)
	{
		// No further split happened?
		if(newOp.seg.vertices.size() == seg.vertices.size())
			return;

        // Too small of a segment
        if(newOp.seg.vertices.size() < small_segment_threshold){
			undecided.push_back(newOp);
            continue;
		}

        // Add it
		children.push_back(newOp);

		children.back().split();
		children.back().seg.pid = seg.uid;
    }

	if(children.empty()) return;

	// Check if segment was split to arbitrary interconnected components
	if( level > 0 && children.size() > 1 )
	{	
		// Pick the two largest children
		std::priority_queue< std::pair<size_t,size_t> > pq;
		for(size_t i = 0; i < children.size(); i++)
			pq.push( std::make_pair(children[i].seg.vertices.size(), i) );
		auto i = pq.top().second; pq.pop();
		auto j = pq.top().second; pq.pop();

		auto bboxi = children[i].bbox(), bboxj = children[j].bbox();
		auto bboxivol = bboxi.volume(), bboxjvol = bboxj.volume();

		// Compute intersection ratio
		auto isect = bboxi.intersection(bboxj);
		auto isectvol = isect.volume();
		double isect_ratio = isectvol / std::min(bboxivol, bboxjvol);

		if( isect_ratio > isect_ratio_threshold ){
			//debugBox( QString("isect_ratio = %1, solidity = %2").arg(isect_ratio).arg(solidity) );
			children.clear();
			return;
		}
	}

	// Merge with neighbor
	for(auto & op : undecided)
	{
		auto & smallSeg = op.seg;

		// Decide which child to assign to
		SplitOperation * selected = &children.front();
		for(auto & child : children){
			if( neiGraph.IsEdgeExists(smallSeg.uid, child.seg.uid) ){
				selected = &child;
				break;
			}
		}

		// Bring back old edges
		for(auto vi : op.seg.vertices){
			for(auto edge : this->seg.adjacency_map[vi])
			{
				auto vj = edge.target;

				// Check if edge is valid
				bool isOutsideNeighbour = (op.seg.vertices.find(vj) == op.seg.vertices.end());

				// Only neighbors in selected child
				if( isOutsideNeighbour && selected->seg.vertices.find(vj) == selected->seg.vertices.end() )
					continue;

				// Add edge back
				selected->seg.AddEdge(vi, vj, edge.weight, edge.index);
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

Eigen::AlignedBox3d SplitOperation::bbox()
{
	Eigen::AlignedBox3d box;
	for(auto f : hull.faces) for(auto v : f) box.extend(v);
	return box;
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
