#include "ParticleCorresponder.h"

ParticleCorresponder::ParticleCorresponder(ParticleMesh *pmeshA, ParticleMesh *pmeshB):
    sA(pmeshA), sB(pmeshB)
{
	for(auto & p : sA->particles) p.correspondence = 0;
	for(auto & p : sB->particles) p.correspondence = 0;

	basicCorrespondence();
}

void ParticleCorresponder::basicCorrespondence()
{
	QVector<ParticleMesh*> inputs;
	inputs << sA << sB;

	// KD-tree	
	QVector<NanoKdTree*> trees;
	for(auto & pmesh : inputs)
	{
		auto relativeKdtree = new NanoKdTree;

		Eigen::AlignedBox3d box = pmesh->bbox();
		Vector3 sizes = box.sizes();

		for( auto & particle : pmesh->particles )
		{
			particle.relativePos = pmesh->relativePos( particle.id );
			relativeKdtree->addPoint( particle.relativePos );
		}

		relativeKdtree->build();

		trees << relativeKdtree;
	}

	for(auto & pmesh : inputs)
	{
		// Select the kd-tree of the other
		auto otherTree = trees.back();
		if(pmesh == inputs.back()) otherTree = trees.front();

		// Match with closest based on euclidean distance
		for(auto & iparticle : pmesh->particles)
		{
			iparticle.correspondence = otherTree->closest( iparticle.relativePos );
		}
	}
}
