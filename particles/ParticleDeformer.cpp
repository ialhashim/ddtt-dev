#include "ParticleDeformer.h"
#include "myglobals.h"

#include "Raytracing.h"

ParticleDeformer::ParticleDeformer(ParticleMesh *pmeshA, ParticleMesh *pmeshB):
    sA(pmeshA), sB(pmeshB)
{
	// Compute sum of errors
	int numSteps = 10;
	std::vector<double> errorVector;

	Eigen::AlignedBox3d box = sA->bbox();

	// Overwrite
	{
		std::vector<Eigen::Vector3f> originalPoints;
		for(auto & particle : sA->particles)
			originalPoints.push_back( particle.pos.cast<float>() );
		auto originalMesh = sA->meshPoints( originalPoints );
		raytracing::Raytracing<Eigen::Vector3d> rtOriginal( originalMesh );
		#pragma omp parallel for
		for(int pi = 0; pi < (int)sA->particles.size(); pi++)
		{
			auto & p = sA->particles[pi];
			int r = 0;
			for(auto d : sA->usedDirections)
			{
				raytracing::RayHit hit = rtOriginal.hit( p.pos, d );
				sA->desc[pi][r++] = hit.distance;
			}
		}
	}

	for(int i = 0; i < numSteps; i++)
	{
		double t = double(i) / (numSteps-1);

		// Mesh at time 't'
		std::vector<Eigen::Vector3d> movedPoints;

		for(auto & particle : sA->particles){
			Vector3 relativePos = AlphaBlend(t, particle.relativePos, sB->particles[particle.correspondence].relativePos);
			movedPoints.push_back( sA->realPos( relativePos ) );
		}

		//auto curMesh = sA->meshPoints( movedPoints );
		auto curMesh = sA->meshPointsUsingGraph( movedPoints, Vector3(0,0,0) );

		// Accelerated raytracing
		raytracing::Raytracing<Eigen::Vector3d> rt( curMesh );

		std::vector<double> errors( sA->particles.size(), 0 );

		int misses = 0;

		// Shoot rays around all particles
		#pragma omp parallel for
		for(int pi = 0; pi < (int)sA->particles.size(); pi++)
		{
			auto & p = sA->particles[pi];

			std::vector<float> oldDescriptor = sA->desc[pi];
			std::vector<float> newDescriptor = sA->desc[pi];

			int r = 0;
			for(auto d : sA->usedDirections)
			{
				raytracing::RayHit hit = rt.hit( movedPoints[pi].cast<double>(), d );
				newDescriptor[r++] = hit.distance;


				#pragma omp critical
				if(!hit.isHit) misses++;
			}

			// Compute error
			auto di = Eigen::Map<Eigen::VectorXf>(&sA->desc[pi][0], sA->desc[pi].size());
			auto dj = Eigen::Map<Eigen::VectorXf>(&newDescriptor[0], newDescriptor.size());

			errors[pi] = (di-dj).norm();
		}

		debugBox(misses);
		curMesh->write((QString("%1_").arg(i) + "step.off").toStdString());

		double errorSum = 0;
		for(auto e : errors) errorSum += e;
		errorVector.push_back( errorSum );
		
		delete curMesh;
	}

	debugBoxVec(errorVector);
}
