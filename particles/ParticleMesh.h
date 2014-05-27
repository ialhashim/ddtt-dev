#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

#include "Particle.h"
#include "NanoKdTree.h"

class ParticleMesh
{
public:
    ParticleMesh( const std::vector<Eigen::Vector3d> & fromPoints, double radius );

    std::vector<Eigen::Vector3d> extractSurface();
	void process();
	void drawParticles();

    std::vector<Particle> particles;
    double raidus;
	Eigen::Vector3d tranlsation;
	Eigen::AlignedBox3d bbox;

	NanoKdTree * kdtree;
};
