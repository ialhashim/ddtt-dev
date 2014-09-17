#include "ParticleCorresponder.h"
#include "PartCorresponder.h"
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
	//descriptorCorrespondence();
	partToPartCorrespondence();

	// Compute relative positions
	for( auto & particle : sA->particles )	particle.relativePos = sA->relativePos( particle.id );
	for( auto & particle : sB->particles )	particle.relativePos = sB->relativePos( particle.id );
}

void ParticleCorresponder::partToPartCorrespondence()
{
	SegmentGraph neiGraphA, neiGraphB;

	auto segmentsA = sA->segmentToComponents(sA->toGraph(), neiGraphA);
	auto segmentsB = sB->segmentToComponents(sB->toGraph(), neiGraphB);

	int count = std::min(segmentsA.keys().size(), segmentsB.keys().size());

	for(int i = 0; i < count ; i++)
	{
		int segA = segmentsA.keys()[i];
		int segB = segmentsB.keys()[i];

		PartCorresponder pc( sA, segmentsA[segA], sB, segmentsB[segB] );
		for(auto d : pc.debug) debug << d;
	}
}

void ParticleCorresponder::descriptorCorrespondence()
{
	typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> MatrixXf;
	typedef Index< L2<float> > FlannIndex;

	int knn = 12;	
	flann::SearchParams params;
	params.cores = omp_get_num_procs();

	QVector<ParticleMesh*> inputs;
	inputs << sA << sB;

	QVector<MatrixXf*> matrices;
	QVector<FlannIndex*> trees;

	for(auto & pmesh : inputs)
	{
		// Fill dataset
		auto m = new MatrixXf( pmesh->particles.size(), pmesh->desc.front().size() );
		for(auto & p : pmesh->particles) 
			m->row(p.id) = Eigen::Map<Eigen::VectorXf>(&pmesh->desc[p.id][0], pmesh->desc[p.id].size());
		Matrix<float> dataset( m->data(), m->rows(), m->cols() );

		// construct index using 4 kd-trees
		auto index = new FlannIndex(dataset, flann::KDTreeIndexParams());
		index->buildIndex();

		trees << index;
		matrices << m;
	}

	for(size_t i = 0; i < inputs.size(); i++)
	{
		auto j = (i+1) % inputs.size();
		auto pi = inputs[i], pj = inputs[j];

		std::vector< std::vector<int> > indices;
		std::vector< std::vector<float> > dists;

		// Query the 'j's in the 'i's
		MatrixXf q( pj->particles.size(), pj->desc.front().size() );
		for(auto & p : pj->particles) q.row(p.id) = Eigen::Map<Eigen::VectorXf>(&pj->desc[p.id][0], pj->desc[p.id].size());
		trees[i]->knnSearch( Matrix<float>( q.data(), q.rows(), q.cols() ), indices, dists, knn, params );

		for(auto & p : pi->particles) p.weight = DBL_MAX;

		for(auto & p : pj->particles)
		{
			for(size_t w = 0; w < indices[p.id].size(); w++)
			{
				int idx = indices[p.id][w];
				double dist = dists[p.id][w];

				if( dist < pi->particles[idx].weight )
				{
					pi->particles[idx].correspondence = p.id;
					pi->particles[idx].weight = dist;
				}
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
