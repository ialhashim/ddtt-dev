#include "ParticleCorresponder.h"
#include "PartCorresponder.h"
#include "Bounds.h"
#include "disjointset.h"
#include "myglobals.h"

ParticleCorresponder::ParticleCorresponder(ParticleMesh *pmeshA, ParticleMesh *pmeshB) : sA(pmeshA), sB(pmeshB)
{
	QVector<ParticleMesh*> input;
	input << sA << sB;

	/// Compute slices for each segment of shape:
	#pragma omp parallel for
	for(int i = 0; i < input.size(); i++)
	{
		SegmentGraph neiGraph;
		auto & segments = input[i]->segmentToComponents( input[i]->toGraph(), neiGraph );

		for(auto & seg : segments)
			seg.property["slices"].setValue( PartCorresponder::computeSlices( input[i], seg ) );
		
		input[i]->property["segments"].setValue( segments );
	}

	// Grouping Experiment
	{
		for(int i = 0; i < input.size(); i++)
		{
			double box_z = input[i]->bbox().sizes().z();
			auto & segments = input[i]->property["segments"].value<Segments>();

			auto segIDs = segments.keys();
			std::map<size_t,size_t> segMap;
			for(auto id : segIDs) segMap[id] = segMap.size();

			NanoKdTree alongZ;
			for( auto & s : segments ){
				Vector3 p(0,0,input[i]->segmentBoundingBox(s).center().z() / box_z);
				alongZ.addPoint( p );
			}
			alongZ.build();

			// Group similar segments
			DisjointSet disjoint( segments.size() );
			
			for(int sid = 0; sid < segIDs.size(); sid++)
			{
				auto & si = segments[ segIDs[sid] ];
				auto box = input[i]->segmentBoundingBox( si );
				Vector3 p(0,0,box.center().z() / box_z);

				KDResults matches;
				alongZ.ball_search(p, 0.01, matches);

				// Check for possible rotational symmetries
				{
					std::vector<size_t> group;

					// Check ground-based similarity
					Boundsd si_measure;
					for(auto v : si.vertices) si_measure.extend( input[i]->particles[v].measure );

					for(auto match : matches)
					{
						if(match.first == sid) continue;
						auto & sj = segments[ segIDs[match.first] ];

						Boundsd sj_measure;
						for(auto v : sj.vertices) sj_measure.extend( input[i]->particles[v].measure );
						double max_range = std::max(si_measure.range(), sj_measure.range());
						double intersection = si_measure.intersection(sj_measure);

						if( intersection / max_range > 0.75 )
							group.push_back( segIDs[match.first] );
					}

					group.push_back( segIDs[sid] );

					if( group.size() > 2 )
					{
						std::vector<Vector3> centroids;
						for(auto sid : group)
						{
							Vector3 center(0,0,0);
							for(auto v : segments[ sid ].vertices) center += input[i]->particles[v].pos;
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
								disjoint.Union(sid, match.first);
							continue;
						}
					}
				}

				// Check for possible reflectional symmetries
				{
					double bbox_threshold = (box.sizes().norm() * 0.1);

					for(auto match : matches)
					{
						if(match.first == sid) continue;
						auto & sj = segments[ segIDs[match.first] ];

						auto boxOther = input[i]->segmentBoundingBox( sj );

						// Join if very similar bounding box-wise
						double diff = (box.sizes() - boxOther.sizes()).norm();
						if(diff > bbox_threshold) continue;

						disjoint.Union(sid, match.first);
					}
				}
			}
			
			for(auto group : disjoint.Groups())
			{
				std::vector<size_t> groupIDs;
				for(auto memeber : group) groupIDs.push_back( segIDs[memeber] );

				auto bs = new starlab::BoxSoup;
				Eigen::AlignedBox3d groupBox;
				for(auto segID : groupIDs) groupBox.extend( input[i]->segmentBoundingBox( segments[segID] ) );
				bs->addBox(groupBox, Qt::black);
				debug << bs;
			}

			break;
		}
	}
}

Particles ParticleCorresponder::partToPartCorrespondence( const QVector< std::pair< size_t,size_t> > & partToPartAssignments )
{
	Particles result;



	return result;
}
