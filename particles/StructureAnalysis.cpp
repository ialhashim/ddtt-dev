#include "StructureAnalysis.h"
#include "myglobals.h"
#include "Planes.h"
#include "Bounds.h"
#include "graph_helper.h"

#include "disjointset.h"

#include "kmedoids.h"

// Declare property types
Q_DECLARE_METATYPE(Eigen::Vector3d)
Q_DECLARE_METATYPE(Eigen::Vector3f)
Q_DECLARE_METATYPE(Eigen::AlignedBox3d)
Q_DECLARE_METATYPE(Boundsd)

StructureAnalysis::StructureAnalysis(ParticleMesh * pmesh) : s(pmesh)
{

	// Debug points
	auto segsCentroid = new starlab::PointSoup(10);

	// Extract groups of similar segments
	std::map<int, std::vector<SegmentGraph*> > groupedSegments;
	SegmentGraph neiGraph;

	auto segments = s->segmentToComponents( neiGraph );

	for(auto & seg : segments)
	{
		int si = s->particles[seg.FirstVertex()].segment;

		auto c_seg = &seg;

		groupedSegments[ si ].push_back( c_seg );

		// Track segment
		allSegs[ c_seg->uid ] = c_seg;

		// Segment properties
		{
			Vector3 sum(0,0,0);
			Eigen::AlignedBox3d bbox;
			Boundsd interval;

			for(auto v : c_seg->vertices) 
			{
				Vector3 p = s->particles[v].pos;
				sum += p;
				bbox.extend( p );
				interval.extend( s->particles[v].measure );
			}

			Vector3 centroid = sum / c_seg->vertices.size();

			// Closest particle to this centroid
			uint centroid_id;
			double minDist = DBL_MAX;
			for(auto v : c_seg->vertices){
				double dist = (s->particles[v].pos - centroid).norm();
				if(dist < minDist){
					minDist = dist;
					centroid_id = v;
				}
			}

			c_seg->property["centroid"].setValue( centroid );
			c_seg->property["centroid_id"].setValue( centroid_id );
			c_seg->property["bbox"].setValue( bbox );
			c_seg->property["interval"].setValue( interval );

			// Visualize centroid
			segsCentroid->addPoint(centroid, Qt::black);
		}
	}

	debug << segsCentroid;

	// Experiment
	{
		std::vector<SegmentGraph*> allSegsVec;
		for(auto seg : allSegs)
			allSegsVec.push_back(seg);

		// Similarity matrix
		int N = allSegsVec.size();
		Eigen::MatrixXf M = Eigen::MatrixXf( N, N );

		for(size_t i = 0; i < N; i++){
			for(size_t j = 0; j < N; j++){
				auto segI = allSegsVec[i];
				auto segJ = allSegsVec[j];

				uint pi = segI->property["centroid_id"].value<uint>();
				uint pj = segJ->property["centroid_id"].value<uint>();

				auto di = Eigen::Map<Eigen::VectorXf>(&s->desc[pi][0], s->desc[pi].size());
				auto dj = Eigen::Map<Eigen::VectorXf>(&s->desc[pj][0], s->desc[pj].size());
				auto dist = (di-dj).norm();

				M(i,j) = dist;
			}
		}

		clustering::kmedoids km( M.rows() );
		km.pam(M, 5, 0);

		for(size_t i = 0; i < N; i++){
			auto newSeg = km.cluster_ids[i];
			auto seg = allSegsVec[i];

			for(auto v : seg->vertices)
				s->particles[v].segment = newSeg;
		}

		// For each level of clustering
		{

		}

		return;
	}


	// Merge similar isolated segments
	{
		segsCentroid->clear();

		std::vector<SegmentGraph*> isolated;
		for(auto & segGroup : groupedSegments){
			//if(segGroup.second.size() > 1) continue;
			//isolated.push_back( segGroup.second.front() );

			for(auto seg : segGroup.second)
				isolated.push_back(seg);
		}

		for(auto & seg : isolated)
			segsCentroid->addPoint(seg->property["centroid"].value<Vector3>(), Qt::black);

		DisjointSet disjointset(isolated.size());

		for(size_t i = 0; i < isolated.size(); i++){
			auto & segI = isolated[i];
			for(size_t j = i + 1; j < isolated.size(); j++){
				auto & segJ = isolated[j];

				// Consider neighbors only
				if( !neiGraph.IsEdgeExists( segI->uid, segJ->uid ) )
					continue;

				Vector3 centroid_i = segI->property["centroid"].value<Vector3>();
				Vector3 centroid_j = segJ->property["centroid"].value<Vector3>();

				uint centroid_id_i = segI->property["centroid_id"].value<uint>();
				uint centroid_id_j = segJ->property["centroid_id"].value<uint>();

				// Height difference
				//if( std::abs((centroid_i - centroid_j).z()) > s->grid.unitlength * 3 )
				//	continue;

				// Flatness
				auto flat_i = s->particles[centroid_id_i].flat;
				auto flat_j = s->particles[centroid_id_j].flat;
				auto flat_similiarty = 1 - abs(flat_i - flat_j);

				if( flat_i > 0.7 && flat_similiarty > 0.7 )
					disjointset.Union(i,j);
			}
		}

		for(size_t i = 0; i < isolated.size(); i++){
			auto & seg = isolated[i];
			for(auto & v : seg->vertices)
				s->particles[v].segment = isolated[disjointset.Parent[i]]->sid;
		}

		return;
	}

	// Relationships
	for(auto edge : neiGraph.GetEdgesSet())
	{
		Vector3 normal = edge.property["normal"].value<Vector3>();
		Vector3 center = edge.property["center"].value<Vector3>();

		starlab::PlaneSoup * ps = new starlab::PlaneSoup(0.05);
		ps->addPlane(center,normal);
		debug << ps;
	}
	
	// For visualization
	if(true)
	{
		QMap<int,QColor> nodeColors;
		QMap<int,QString> nodeLabels;
		for(auto e : neiGraph.GetEdgesSet())
		{
			auto si = allSegs[e.index];
			auto sj = allSegs[e.target];

			nodeColors[e.index] = s->rndcolors[si->sid];
			nodeColors[e.target] = s->rndcolors[sj->sid];

			nodeLabels[e.index] = QString("%1 (%2)").arg(si->sid).arg(si->uid);
			nodeLabels[e.target] = QString("%1 (%2)").arg(sj->sid).arg(sj->uid);
		}

		// Show segments graph
		static SvgView * viewer = new SvgView;
		viewer->show( buildSVG( toGraphvizFormat(neiGraph, nodeColors, nodeLabels) ) );
	}
	/*
	//double threshold = s->grid.unitlength * 3; // two voxels

	pcount = 0;

	std::vector<Plane> allPlanes;

	for(auto & segGroup : groupedSegments)
	{
		auto & group = segGroup.second;

		std::vector<Plane> groupPlanes;

		// Compare pair of parts in the same segment
		for(size_t i = 0; i < group.size(); i++)
		{
			for(size_t j = i+1; j < group.size(); j++)
			{
				auto groupI = group[i];
				auto groupJ = group[j];

				Vector3 centroid_i = groupI->property["centroid"].value<Vector3>();
				Vector3 centroid_j = groupJ->property["centroid"].value<Vector3>();

				uint centroid_id_i = groupI->property["centroid_id"].value<uint>();
				uint centroid_id_j = groupJ->property["centroid_id"].value<uint>();

				Vector3 planeCenter = (centroid_i+centroid_j) * 0.5;
				Vector3 planeNormal = (centroid_i-centroid_j).normalized();

				/// Filter:

				// by measure
				{
					uint closest_i, closest_j;

					// Compare with respect to measure
					double difference = abs(s->particles[centroid_id_i].measure - s->particles[centroid_id_j].measure);
					double similarity = 1-difference;

					if( similarity < 0.9 ) continue;
				}

				// by mass
				{
					Boundsd interval_i = groupI->property["interval"].value< Boundsd >();
					Boundsd interval_j = groupJ->property["interval"].value< Boundsd >();

					double minVol = std::min(interval_i.count, interval_j.count);
					double maxVol = std::max(interval_i.count, interval_j.count);
					double ratio = minVol / maxVol;

					if(ratio < 0.1) continue;

					groupPlanes.push_back( Plane(planeCenter, planeNormal, ratio) );
				}

				// by bounding volume
				{
					Eigen::AlignedBox3d bbox_i = groupI->property["bbox"].value< Eigen::AlignedBox3d >();
					Eigen::AlignedBox3d bbox_j = groupJ->property["bbox"].value< Eigen::AlignedBox3d >();
					double bbox_diff = (bbox_i.sizes() - bbox_j.sizes()).norm();
				}

				// DEBUG:
				{
					starlab::PlaneSoup * ps = new starlab::PlaneSoup(0.01);
					ps->addPlane( planeCenter, planeNormal );
					debug << ps;
					debug << new RenderObject::Segment(centroid_i, centroid_j, 1);
				}

				//Force segment:
				//for(auto v : groupI->vertices) s->particles[v].segment = 50;
				//for(auto v : groupJ->vertices) s->particles[v].segment = 50;

				pcount++;
			}
		}

		auto mergedPlanes = Plane::mergePlanes(groupPlanes, 0.99);
		for (auto & plane : mergedPlanes) allPlanes.push_back(plane);
	}

	auto rndColors = rndColors2(pcount);

	auto mergedPlanes = Plane::mergePlanes(allPlanes, 0.99);

	for(int i = 0; i < (int)mergedPlanes.size(); i++)
	{
		auto & plane = mergedPlanes[i];
		auto color = QColor::fromRgbF(rndColors[i].redF(),rndColors[i].greenF(),rndColors[i].blueF());
		starlab::PlaneSoup * ps = new starlab::PlaneSoup(0.1 * plane.weight, true, color);
		ps->addPlane( plane.pos, plane.n );

		//debug << ps;
	}
	*/
}
