#include "ParticleDeformer.h"
#include "myglobals.h"

#include "Raytracing.h"

ParticleDeformer::ParticleDeformer(ParticleMesh *pmeshA, ParticleMesh *pmeshB):
    sA(pmeshA), sB(pmeshB)
{
	// Compute sum of errors
	int numSteps = 10;
	std::vector<double> errorVector;

	Eigen::AlignedBox3d boxA( sA->bbox_min, sA->bbox_max );
	Eigen::AlignedBox3d boxB( sB->bbox_min, sB->bbox_max );
	Vector3 bbox_ratios(0,0,0);
	for(int i = 0; i < 3; i++){
		bbox_ratios[i] = boxA.sizes()[i] / boxB.sizes()[i];
		if(bbox_ratios[i] <= 1.49) bbox_ratios[i] = 0;
	}

	// Compute original descriptors
	{
		std::vector<Vector3> originalPoints;
		for(auto & particle : sA->particles)
			originalPoints.push_back( particle.pos );
		auto originalMesh = sA->meshPoints( originalPoints, Eigen::Vector3i(0,0,0) );
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

		// Clean up
		delete originalMesh;
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

		auto curMesh = sA->meshPoints( movedPoints, bbox_ratios.cast<int>() );

		// Accelerated raytracing
		raytracing::Raytracing<Eigen::Vector3d> rt( curMesh );

		std::vector<double> errors( sA->particles.size(), 0 );
		int misses = 0;

		double unitlength = sA->grid.unitlength;
		Vector3 center_delta = sA->grid.translation.cast<double>() + ( 0.5 * Vector3(unitlength,unitlength,unitlength) );

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
				// Snap point to its location on the grid
				Vector3 pg = (movedPoints[pi] - sA->grid.translation.cast<double>()) / unitlength;
				Eigen::Vector3i v( pg.x(), pg.y(), pg.z() );
				Vector3 gridpnt_pos = Vector3(v[0] * unitlength, v[1] * unitlength, v[2] * unitlength) + center_delta;
			
				raytracing::RayHit hit = rt.hit( gridpnt_pos, d );
				newDescriptor[r++] = hit.distance;

				if( !hit.isHit ) 
				{
					#pragma omp critical
					{
						misses++;
						auto vs = new starlab::VectorSoup;
						vs->addVector(movedPoints[pi].cast<double>(), d);
						debug << vs;
						//curMesh->write("missMesh.off");
					}
				}
			}

			// Compute error
			auto di = Eigen::Map<Eigen::VectorXf>(&sA->desc[pi][0], sA->desc[pi].size());
			auto dj = Eigen::Map<Eigen::VectorXf>(&newDescriptor[0], newDescriptor.size());

			errors[pi] = (di-dj).norm();
		}

		if( misses ){
			debugBox(QString("Missed rays = %1").arg(misses));
			curMesh->write((QString("%1_").arg(i) + "step.off").toStdString());
		}

		double errorSum = 0;
		for(auto e : errors) errorSum += e;
		errorVector.push_back( errorSum );
		
		delete curMesh;
	}

	debugBoxVec(errorVector);
}
