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

class ParticleMesh : public Serializable
{
public:

	ParticleMesh(SurfaceMeshModel * mesh = NULL, int gridsize = 64);
	~ParticleMesh();
	PropertyMap property;

	std::vector< Particle<Vector3> > particles;
	std::vector< std::vector<float> > desc, sig;

	SurfaceMeshModel * surface_mesh;

	int sphereResolutionUsed;
	std::vector< Vector3 > usedDirections;
	std::vector< size_t > antiRays;

	typedef Eigen::Vector3f VoxelVector;
	VoxelContainer<VoxelVector> grid;
	std::map<uint64_t,size_t> mortonToParticleID;

public:
	Vector3 mainDirection( size_t particleID );

	void process();
	void computeDistanceToFloor();
	std::vector<size_t> pathFromFloor;

	void cluster( int K, const std::vector<size_t> & seeds, bool use_l1_norm, bool showSeeds );
	std::vector<VectorFloat> cluster_centers;
	void shrinkSmallerClusters();

	enum SeedType{ RANDOM, PLUS_PLUS, GROUND, DESCRIPTOR };
	std::vector<size_t> specialSeeding( SeedType seedType, int K, SegmentGraph::vertices_set selected = SegmentGraph::vertices_set() );

	SegmentGraph toGraph( SegmentGraph::vertices_set selected = SegmentGraph::vertices_set() );

	QMap< unsigned int, SegmentGraph > segmentToComponents( SegmentGraph fromGraph, SegmentGraph & neiGraph );
	std::vector< std::vector< std::vector<float> > > toGrid();
	SpatialHash< Vector3, Vector3::Scalar > spatialHash();
	std::vector<size_t> randomSamples( int numSamples, bool isSpread );
	std::vector< double > agd( int numStartPoints );
	std::vector<size_t> neighbourhood( Particle<Vector3> & p, int step = 2);
	Particle<Vector3> pointToParticle( const Vector3 & point );
	std::vector< Vector3 > particlesCorners( SegmentGraph::vertices_set selected = SegmentGraph::vertices_set() );
	std::vector< Vector3 > particlesPositions(const std::set<unsigned int> & P);

	std::vector< std::pair< double, size_t > > closestParticles( const Vector3 & point, double threshold = 1e12 );

	Eigen::AlignedBox3d bbox();
	Vector3 bbox_min, bbox_max;

	Vector3 relativePos( size_t particleID );
	Vector3 realPos( Vector3 relative_pos );

	SurfaceMeshModel * meshPoints( const std::vector<Vector3> & points, Eigen::Vector3i ratios ) const;

	void distort();

	std::vector< std::map< int, std::vector<size_t> > > cachedAdj;
	SegmentGraph cachedGraph;

	static QVector<QColor> rndcolors;
	QVector<RenderObject::Base*> debug;
	void drawParticles( qglviewer::Camera * camera );
	void drawDebug(QGLWidget & widget);

	// Serialization:
	void serialize(QDataStream& os) const;	
	void deserialize(QDataStream&);
};
