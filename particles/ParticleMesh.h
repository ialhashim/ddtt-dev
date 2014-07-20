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

#include "SpatialHash.h"

#include "GenericGraph.h"
typedef GenericGraphs::Graph<uint,double> SegmentGraph;

class ParticleMesh
{
public:

	ParticleMesh(SurfaceMeshModel * mesh, int gridsize = 64, double particle_raidus = 0.1);
	~ParticleMesh();

	std::vector< Particle<Vector3> > particles;
	double raidus;
	std::vector< std::map< int, std::vector<size_t> > > cachedAdj;

	void process();
	void computeDistanceToFloor();
	std::vector<size_t> pathFromFloor;

	void cluster( int K, const std::set<size_t> & seeds, bool use_l1_norm );
	void shrinkSmallerClusters();

	enum GraphEdgeWeight{ GEW_DISTANCE, GEW_DIAMETER };
	SegmentGraph toGraph( GraphEdgeWeight wtype = GEW_DISTANCE );
	std::vector< SegmentGraph > segmentToComponents( SegmentGraph & neiGraph );

	std::vector< std::vector< std::vector<float> > > toGrid();
	SpatialHash< Vector3, Vector3::Scalar > spatialHash();
	std::vector<size_t> randomSamples( int numSamples, bool isSpread );
	std::vector<size_t> neighbourhood( Particle<Vector3> & p, int step );
	std::vector< double > agd( int numStartPoints );

	std::vector< std::pair< double, size_t > > closestParticles( const Vector3 & point, double threshold = 1e12 );

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
