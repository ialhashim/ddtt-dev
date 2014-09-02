#include "ParticleCorresponder.h"

ParticleCorresponder::ParticleCorresponder(ParticleMesh *pmeshA, ParticleMesh *pmeshB):
    sA(pmeshA), sB(pmeshB)
{
    /*

	// KD-tree
	{
		relativeKdtree = new NanoKdTree;

		Eigen::AlignedBox3d box = bbox();
		Vector3 sizes = box.sizes();

		for( auto & particle : particles )
		{
		Vector3 mapped = (particle.pos - box.min());
		for(int i = 0; i < 3; i++) mapped[i] /= sizes[i];

		particle.relativePos = mapped;
		relativeKdtree->addPoint( particle.relativePos );
		}

		relativeKdtree->build();
	}

	
	NanoKdTree * itree = sA->relativeKdtree;
    NanoKdTree * jtree = sB->relativeKdtree;

    for(auto & iparticle : sA->particles)
    {
        if( false )
        {
            iparticle.correspondence = jtree->closest( iparticle.relativePos );
        }
        else
        {
            // Experiment
            KDResults matches;
            jtree->k_closest( iparticle.relativePos, 20, matches );

            QMap<double, size_t> measures;
            for(auto p : matches)
            {
                double weight = p.second;
                //double desc_dist = dist_fn()( sA->desc[iparticle.id], sB->desc[p.first] );
                measures[ weight ] = p.first;
            }
            iparticle.correspondence = measures[ measures.keys().front() ];
        }
    }

    for(auto & jparticle : sB->particles)
    {
        if( true )
        {
            jparticle.correspondence = itree->closest( jparticle.relativePos );
        }
        else
        {
            // Experiment
            KDResults matches;
            itree->k_closest( jparticle.relativePos, 20, matches );

            QMap<double, size_t> measures;
            for(auto p : matches)
            {
                double weight = p.second;
                //double desc_dist = dist_fn()( sB->desc[jparticle.id], sA->desc[p.first] );
                measures[ weight ] = p.first;
            }
            jparticle.correspondence = measures[ measures.keys().front() ];
        }
    }*/
}
