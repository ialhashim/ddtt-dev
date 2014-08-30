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
typedef std::vector<float> VectorFloat;

class ParticleMesh
{
public:

	ParticleMesh(SurfaceMeshModel * mesh, int gridsize = 64, double particle_raidus = 0.1);
	~ParticleMesh();
	PropertyMap property;

	std::vector< Particle<Vector3> > particles;
	double raidus;
	std::vector< std::map< int, std::vector<size_t> > > cachedAdj;

	std::vector< Vector3 > usedDirections;
	std::vector< size_t > antiRays;
	Vector3 mainDirection( size_t particleID );

	void process();
	void computeDistanceToFloor();
	std::vector<size_t> pathFromFloor;

	std::vector<VectorFloat> cluster_centers;
	void cluster( int K, const std::set<size_t> & seeds, bool use_l1_norm, bool showSeeds );
	void shrinkSmallerClusters();

	enum SeedType{ RANDOM, PLUS_PLUS, GROUND, DESCRIPTOR };
	std::set<size_t> specialSeeding( SeedType seedType, int K, SegmentGraph::vertices_set selected = SegmentGraph::vertices_set() );

	SegmentGraph toGraph( SegmentGraph::vertices_set selected = SegmentGraph::vertices_set() );
	SegmentGraph cachedGraph;

	QMap< unsigned int, SegmentGraph > segmentToComponents( SegmentGraph fromGraph, SegmentGraph & neiGraph );
	std::vector< std::vector< std::vector<float> > > toGrid();
	SpatialHash< Vector3, Vector3::Scalar > spatialHash();
	std::vector<size_t> randomSamples( int numSamples, bool isSpread );
	std::vector< double > agd( int numStartPoints );
	std::vector<size_t> neighbourhood( Particle<Vector3> & p, int step = 2);
	Particle<Vector3> pointToParticle( const Vector3 & point );
	std::vector< Vector3 > particlesCorners( SegmentGraph::vertices_set selected = SegmentGraph::vertices_set() );

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

	SurfaceMeshModel * meshPoints( const std::vector<Eigen::Vector3f> & points ) const;
	QVector<RenderObject::Base*> debug;
	void distort();
};
