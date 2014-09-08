#include "ParticleCorresponder.h"
#include "Bounds.h"

#pragma warning( disable : 4005 4100 4267 4616 4700 ) // its not my code..
#include "flann/flann.hpp"
using namespace flann;

ParticleCorresponder::ParticleCorresponder(ParticleMesh *pmeshA, ParticleMesh *pmeshB):
    sA(pmeshA), sB(pmeshB)
{
	for(auto & p : sA->particles) p.correspondence = 0;
	for(auto & p : sB->particles) p.correspondence = 0;

	//basicCorrespondence();

	descriptorCorrespondence();
}

void ParticleCorresponder::descriptorCorrespondence()
{
	size_t selectedID = 0;
	auto & selected = sA->particles[selectedID];
	auto ps = new starlab::PointSoup(15);
	ps->addPoint( selected.pos );
	debug << ps;

	Boundsd bound;
	for(auto & p : sB->particles)
	{
		auto di = Eigen::Map<Eigen::VectorXf>(&sA->desc[selectedID][0], sA->desc[selectedID].size());
		auto dj = Eigen::Map<Eigen::VectorXf>(&sB->desc[p.id][0], sB->desc[p.id].size());

		bound.extend( (di - dj).norm() );
	}

	for(auto & p : sB->particles)
	{
		p.flag = VIZ_WEIGHT;

		auto di = Eigen::Map<Eigen::VectorXf>(&sA->desc[selectedID][0], sA->desc[selectedID].size());
		auto dj = Eigen::Map<Eigen::VectorXf>(&sB->desc[p.id][0], sB->desc[p.id].size());

		double dist = (di - dj).norm();

		p.weight = bound.normalized( dist );
	}
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
