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

#include "SpatialHash.h"

class ParticleMesh
{
public:
	ParticleMesh(SurfaceMeshModel * mesh, int gridsize = 64, double particle_raidus = 0.1);
	~ParticleMesh();

	std::vector< Particle<Vector3> > particles;
	double raidus;

	void process();
	void computeDistanceToFloor();

	enum GraphEdgeWeight{ GEW_DISTANCE, GEW_DIAMETER };
	GenericGraphs::Graph<uint,double> & toGraph( GraphEdgeWeight wtype = GEW_DISTANCE );
	GenericGraphs::Graph<uint,double> cachedGraph;
	std::vector< GenericGraphs::Graph<uint,double> > segmentToComponents( GenericGraphs::Graph<uint,double> & neiGraph );

	std::vector< std::vector< std::vector<float> > > toGrid();
	SpatialHash< Vector3, Vector3::Scalar > spatialHash();
	std::vector<size_t> randomSamples( int numSamples, bool isSpread );
	std::vector<size_t> neighbourhood( const Particle<Vector3> & p, int step );
	std::vector< double > agd( int numStartPoints );

	void drawParticles( qglviewer::Camera * camera );
	void drawDebug(QGLWidget & widget);
	static QVector<QColor> rndcolors;

	typedef Eigen::Vector3f VoxelVector;
	VoxelContainer<VoxelVector> grid;
	std::map<uint64_t,size_t> mortonToParticleID;

	std::vector< std::vector<float> > desc, sig;

	Eigen::AlignedBox3d bbox();

	SurfaceMeshModel * surface_mesh;
	NanoKdTree * relativeKdtree;

	QVector<RenderObject::Base*> debug;
	void distort();
};
