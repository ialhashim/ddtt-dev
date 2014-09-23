#include "CorrespondencePrepare.h"
#include "PartCorresponder.h"
#include "Bounds.h"
#include "disjointset.h"
#include "myglobals.h"

CorrespondencePrepare::CorrespondencePrepare( std::vector<ParticleMesh*> meshes )
{
	QVector<ParticleMesh*> input;
	for(auto pm : meshes) input << pm;

	#pragma omp parallel for
	for(int i = 0; i < input.size(); i++)
	{
		/// Compute slices for each segment of shape:
		{
			SegmentGraph neiGraph;
			auto & segments = input[i]->segmentToComponents( input[i]->toGraph(), neiGraph );

			for(auto & seg : segments)
				seg.property["slices"].setValue( PartCorresponder::computeSlices( input[i], seg ) );

			for(auto & seg : segments)
				seg.property["bbox"].setValue( input[i]->segmentBoundingBox(seg) );

			input[i]->property["segments"].setValue( segments );
		}

		/// Compute groups of each shape:
		{
			input[i]->property["groups"].setValue( computeGroups(input[i]) );
		}
	}
}

std::vector< std::vector<size_t> > CorrespondencePrepare::computeGroups( ParticleMesh * input )
{
	auto inputbox = input->bbox().sizes();

	auto & segments = input->property["segments"].value<Segments>();

	auto segIDs = segments.keys();
	std::map<size_t,size_t> segMap;
	for(auto id : segIDs) segMap[id] = segMap.size();

	DisjointSet disjoint( segments.size() );
	QMap<size_t,bool> seen;

	// Group similar segments
	for(int coord = 2; coord >= 0; coord--)
	{
		NanoKdTree alongAxis;
		for( auto & seg : segments ){
			Vector3 p(0,0,seg.property["bbox"].value<Eigen::AlignedBox3d>().center()[coord] / inputbox[coord]);
			alongAxis.addPoint( p );
		}
		alongAxis.build();

		for(int sid = 0; sid < segIDs.size(); sid++)
		{
			if( seen[sid] ) continue;

			auto & si = segments[ segIDs[sid] ];
			auto box = si.property["bbox"].value<Eigen::AlignedBox3d>();
			Vector3 p(0,0,box.center()[coord] / inputbox[coord]);

			KDResults matches;
			alongAxis.ball_search(p, 0.01, matches);

			// Check for possible rotational symmetries
			if(coord == 2)
			{
				std::vector<size_t> group;

				// Check ground-based similarity
				Boundsd si_measure;
				for(auto v : si.vertices) si_measure.extend( input->particles[v].measure );

				for(auto match : matches)
				{
					if(seen[match.first]) continue;
					if(match.first == sid) continue;
					auto & sj = segments[ segIDs[match.first] ];

					Boundsd sj_measure;
					for(auto v : sj.vertices) sj_measure.extend( input->particles[v].measure );
					double max_range = std::max(si_measure.range(), sj_measure.range());
					double intersection = si_measure.intersection(sj_measure);

					if( intersection / max_range > 0.75 )
						group.push_back( segIDs[match.first] );
				}

				group.push_back( segIDs[sid] );

				if( group.size() > 2 && group.size() != 4 )
				{
					std::vector<Vector3> centroids;
					for(auto sid : group)
					{
						Vector3 center(0,0,0);
						for(auto v : segments[ sid ].vertices) center += input->particles[v].pos;
						center /= segments[ sid ].vertices.size();
						centroids.push_back( center );
					}

					Vector3 center(0,0,0);
					for(auto c : centroids) center += c;
					center /= centroids.size();

					// Check angle between closest two instances
					double smallestAngle = DBL_MAX;
					Vector3 v1 = (centroids[0] - center).normalized();
					for(size_t b = 1; b < centroids.size(); b++)
					{
						Vector3 v2 = (centroids[b] - center).normalized();
						smallestAngle = std::min(smallestAngle, acos( v1.dot(v2) ) );
					}

					double angle = rad_to_deg( smallestAngle );
					double expectedAngle = 360.0 / group.size();

					if(abs(angle - expectedAngle) < 10)
					{
						for(auto match : matches)
						{
							disjoint.Union(sid, match.first);
							seen[ match.first ] = true;
						}
						continue;
					}
				}
			}

			// Check for possible reflectional symmetries
			{
				double bbox_threshold = (box.sizes().norm() * 0.1);

				for(auto match : matches)
				{
					if(seen[match.first]) continue;
					if(match.first == sid) continue;
					auto & sj = segments[ segIDs[match.first] ];

					auto boxOther = input->segmentBoundingBox( sj );

					// Join if very similar bounding box-wise
					double diff = (box.sizes() - boxOther.sizes()).norm();
					if(diff > bbox_threshold) continue;

					disjoint.Union(sid, match.first);
					seen[ match.first ] = true;
				}
			}
		}
	}

	std::vector< std::vector<size_t> > foundGroups;
	for(auto group : disjoint.Groups())
	{
		std::vector<size_t> groupIDs;
		for(auto memeber : group) groupIDs.push_back( segIDs[memeber] );
		foundGroups.push_back(groupIDs);
	}

	// DEBUG:
	if( false )
	{	
		for(auto group : foundGroups)
		{
			auto bs = new starlab::BoxSoup;
			Eigen::AlignedBox3d groupBox;
			for(auto segID : group) groupBox.extend( input->segmentBoundingBox( segments[segID] ) );
			bs->addBox(groupBox, Qt::black);
			debug << bs;
		}
	}

	return foundGroups;
}
