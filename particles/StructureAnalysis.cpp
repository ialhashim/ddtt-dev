#include <QMessageBox>
#include "StructureAnalysis.h"
#include "myglobals.h"
#include "Planes.h"
#include "Bounds.h"
#include "graph_helper.h"

#include "disjointset.h"

#include "convexhull.h"

// Declare property types
Q_DECLARE_METATYPE(Eigen::Vector3d)
Q_DECLARE_METATYPE(Eigen::Vector3f)
Q_DECLARE_METATYPE(Eigen::AlignedBox3d)
Q_DECLARE_METATYPE(Boundsd)

#include "kmeans.h"
double solidity_threshold = 0.7;
int max_parts_threshold = 12;
int level_threshold = 3;

auto rcolors = rndColors2(200);
int rc = 0;
QVector<RenderObject::Base*> globalDebug;

StructureAnalysis::StructureAnalysis(ParticleMesh * pmesh) : s(pmesh)
{
	rc = 0; // reset random color
	globalDebug.clear();

	struct SplitOperation{
		SegmentGraph seg, neiGraph;
		double solidity;
		int level;
		std::vector<SplitOperation> children;
		SplitOperation * parent;
		ConvexHull<Vector3> hull;
		ParticleMesh * s;

		typedef std::vector<float> FloatVec;
		typedef std::vector<FloatVec> FloatVecVec;
		typedef clustering::lpnorm< VectorFloat > DestFn;

		SplitOperation( ParticleMesh * s, const SegmentGraph& seg, int level = 0, SplitOperation* parent = NULL) : 
			s(s), seg(seg), level(level)
		{
			// Compute solidity score
			hull = ConvexHull<Vector3>( s->particlesCorners( seg.vertices ), "FA Qt" );

			// Segment solidity (Volume / Convex Volume)
			double seg_volume = pow(s->grid.unitlength,3) * seg.vertices.size();
			double convex_volume = hull.volume;
			solidity = seg_volume / convex_volume;
		}

		void split()
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
				clustering::kmeans<FloatVecVec,DestFn> km( descs, 2 );
				km._centers.clear();
				for(auto pid : s->specialSeeding(ParticleMesh::DESCRIPTOR, 2, seg.vertices)) km._centers.push_back( s->desc[pid] );
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
					// No further split should happen
					if(newOp.seg.vertices.size() == seg.vertices.size()){
						children.clear();
						break;
					}

					// Too small of a segment
					if(newOp.seg.vertices.size() < 4)
						continue;

					children.push_back(newOp);

					if(newOp.solidity < solidity_threshold)
					{
						children.back().split();
						children.back().seg.pid = seg.uid;
					}
				}
			}
		}

		void collectClusters( std::vector<SplitOperation*> & ops )
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

		void debugAllChildren(QVector<RenderObject::Base*> & debug)
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

		void report( QStringList & items )
		{
			QString r = QString("[%1 - level %2] solidity (%3) parent(%4) size(%5)").arg(
				seg.uid).arg(level).arg(solidity).arg(seg.pid).arg(seg.vertices.size());

			if(!children.empty())
				r += QString(" children(%1)").arg(children.size());

			items << r;

			for(auto & child : children)
				child.report( items );
		}
	};

	SplitOperation op( pmesh, s->toGraph() );
	op.split();
	
	std::vector<SplitOperation*> clusters;
	op.collectClusters( clusters );

	std::map<size_t,size_t> mappedClusters;
	for(auto c : clusters)
		mappedClusters[c->seg.uid] = mappedClusters.size();

	for(auto c : clusters)
		for(auto v : c->seg.vertices)
			s->particles[v].segment = mappedClusters[c->seg.uid];
		
	/*
	op.debugAllChildren(debug);
	QStringList rep; op.report(rep);
	saveToTextFile( "_split_report.txt", rep );
	for(auto d : globalDebug) debug << d;*/
}
