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
				auto box = input[i]->segmentBoundingBox( segments[ segIDs[sid] ] );
				Vector3 p(0,0,box.center().z() / box_z);

				KDResults matches;
				alongZ.ball_search(p, 0.1, matches);

				double threshold = (box.sizes().norm() * 0.1);

				for(auto match : matches)
				{
					if(match.first == sid) continue;
					auto boxOther = input[i]->segmentBoundingBox( segments[ segIDs[match.first] ] );

					// Join if very similar
					double diff = (box.sizes() - boxOther.sizes()).norm();
					if(diff > threshold) continue;

					disjoint.Union(sid, match.first);
				}
			}
			
			for(auto group : disjoint.Groups())
			{
				std::vector<size_t> groupIDs;
				for(auto memeber : group) groupIDs.push_back( segIDs[memeber] );

				auto bs = new starlab::BoxSoup;
				Eigen::AlignedBox3d groupBox;
				for(auto segID : groupIDs) groupBox.extend( input[i]->segmentBoundingBox( segments[segID] ) );
				bs->addBox(groupBox);
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
