#include "ParticleCorresponder.h"
#include "Bounds.h"

#pragma warning( disable : 4005 4100 4267 4616 4700 4291 4189 4267 4996 4267 ) // its not my code..
#include "flann/flann.hpp"
using namespace flann;

ParticleCorresponder::ParticleCorresponder(ParticleMesh *pmeshA, ParticleMesh *pmeshB):
    sA(pmeshA), sB(pmeshB)
{
	// Clear correspondences
	for(auto & p : sA->particles) p.correspondence = 0;
	for(auto & p : sB->particles) p.correspondence = 0;

	//basicCorrespondence();

	descriptorCorrespondence();

	// Compute relative positions
	for( auto & particle : sA->particles )	particle.relativePos = sA->relativePos( particle.id );
	for( auto & particle : sB->particles )	particle.relativePos = sB->relativePos( particle.id );
}

void ParticleCorresponder::descriptorCorrespondence()
{
	typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> MatrixXf;

	// Fill dataset
	MatrixXf m( sB->particles.size(), sB->desc.front().size() );
	for(auto & p : sB->particles) 
		m.row(p.id) = Eigen::Map<Eigen::VectorXf>(&sB->desc[p.id][0], sB->desc[p.id].size());
	Matrix<float> dataset( m.data(), m.rows(), m.cols() );

	// construct index using 4 kd-trees
	Index< L2<float> > index(dataset, flann::KDTreeIndexParams());
	index.buildIndex();

	int knn = 12;	
	flann::SearchParams params;
	params.cores = omp_get_num_procs();

	MatrixXf q( sA->particles.size(), m.cols() );
	for(auto & p : sA->particles) 
		q.row(p.id) = Eigen::Map<Eigen::VectorXf>(&sA->desc[p.id][0], sA->desc[p.id].size());
	Matrix<float> queries( q.data(), q.rows(), q.cols() );

	std::vector< std::vector<int> > indices;
	std::vector< std::vector<float> > dists;

	index.knnSearch( queries, indices, dists, knn, params );

	for(auto & p : sB->particles) p.weight = DBL_MAX;

	for(auto & p : sA->particles)
	{
		for(size_t j = 0; j < indices[p.id].size(); j++)
		{
			int idx = indices[p.id][j];
			double dist = dists[p.id][j];

			if( dist < sB->particles[idx].weight )
			{
				sB->particles[idx].correspondence = p.id;
				sA->particles[p.id].correspondence = idx;

				sB->particles[idx].weight = dist;
			}
		}
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
