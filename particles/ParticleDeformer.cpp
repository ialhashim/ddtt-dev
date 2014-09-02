#include "ParticleDeformer.h"
#include "myglobals.h"

#include "Raytracing.h"

ParticleDeformer::ParticleDeformer(ParticleMesh *pmeshA, ParticleMesh *pmeshB):
    sA(pmeshA), sB(pmeshB)
{
	// Compute sum of errors
	int numSteps = 10;
	std::vector<double> errorVector;

	for(int i = 0; i < numSteps; i++)
	{
		double t = double(i) / (numSteps-1);

		// Mesh at time 't'
		std::vector<Eigen::Vector3f> movedPoints;

		for(auto & particle : sA->particles){
			Vector3 p = AlphaBlend(t, particle.pos, sB->particles[particle.correspondence].pos);
			movedPoints.push_back( p.cast<float>() );
		}

		auto curMesh = sA->meshPoints( movedPoints );

		// Accelerated raytracing
		raytracing::Raytracing<Eigen::Vector3d> rt( curMesh );

		std::vector<double> errors( sA->particles.size(), 0 );
		std::vector< std::vector<float> > newDescriptor = sA->desc;

		int numUsedParticles = 0;

		// Shoot rays around all particles
		#pragma omp parallel for
		for(int pi = 0; pi < (int)sA->particles.size(); pi++)
		{
			if( !sA->particles[pi].isMedial ) continue;

			auto & p = sA->particles[pi];

			int r = 0;
			for(auto d : sA->usedDirections)
			{
				raytracing::RayHit hit = rt.hit( p.pos, d );
				if(!hit.isHit) hit.distance = sA->grid.unitlength;
				
				newDescriptor[pi][r++] = hit.distance;
			}

			// Compute error
			auto di = Eigen::Map<Eigen::VectorXf>(&sA->desc[pi][0], sA->desc[pi].size());
			auto dj = Eigen::Map<Eigen::VectorXf>(&newDescriptor[pi][0], newDescriptor[pi].size());

			errors[pi] = (di-dj).norm();

			numUsedParticles++;
		}

		double errorSum = 0;
		for(auto e : errors) errorSum += e;
		errorVector.push_back( errorSum / numUsedParticles );
		
		delete curMesh;
	}

	debugBoxVec(errorVector);
}
