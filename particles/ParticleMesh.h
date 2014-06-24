#pragma once
#include <QWidget>
#include "qglviewer/qglviewer.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <random>

#include "Particle.h"
#include "NanoKdTree.h"
#include "RenderObjectExt.h"

#include "voxelization.h"

#include "GenericGraph.h"

class ParticleMesh
{
public:
	ParticleMesh(SurfaceMeshModel * mesh, int gridsize = 64, double particle_raidus = 0.1);
	~ParticleMesh();

	void process();

	std::vector< std::vector< std::vector<float> > > toGrid();
	GenericGraphs::Graph<uint,double> toGraph();

	std::vector< double > agd( int numStartPoints );

	void drawParticles( qglviewer::Camera * camera );
	void drawDebug(QGLWidget & widget);
	static QVector<QColor> rndcolors;

	typedef Eigen::Vector3f VoxelVector;
	VoxelContainer<VoxelVector> grid;
	std::map<uint64_t,size_t> mortonToParticleID;

	std::vector< std::vector<float> > desc, sig;

    std::vector<Particle> particles;
    double raidus;

	Eigen::AlignedBox3d bbox();

	SurfaceMeshModel * surface_mesh;
	NanoKdTree * relativeKdtree;

	QVector<RenderObject::Base*> debug;
	void distort();
};
