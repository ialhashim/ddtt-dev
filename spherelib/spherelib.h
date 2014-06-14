#pragma once

#include <list>
#include <map>
#include "SurfaceMeshModel.h"

#define radToDegrees(angleRadians) ((angleRadians) * 180.0 / M_PI)
#define degreesToRad(angleDegrees) (angleDegrees * M_PI / 180.0)

namespace Spherelib
{

// Icosahedron subdivision from http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
struct SphereMaker
{
    typedef SurfaceMesh::SurfaceMeshModel::Vector3 Vector3;

    struct TriangleIndices{
        int v1,v2,v3;
        TriangleIndices(int v1, int v2, int v3) : v1(v1), v2(v2), v3(v3){}
    };

    int index;
    std::map<size_t, int> middlePointIndexCache;
    SurfaceMesh::SurfaceMeshModel geometry;

    // add vertex to mesh, fix position to be on unit sphere, return index
    int addVertex(const Vector3 & p);

    // return index of point in the middle of p1 and p2
    int getMiddlePoint(int p1, int p2);

    SphereMaker( int recursionLevel );
    static Surface_mesh makeSphere( int recursionLevel ){
        SphereMaker sm( recursionLevel );
        return sm.geometry;
    }
};


struct RadialGrid{    
	typedef SurfaceMesh::SurfaceMeshModel::Vector3 Vector3;

	int sectors, tracks;
	std::vector< std::vector<double> > front_values;
	std::vector< std::vector<double> > back_values;
	double average;

	RadialGrid(){ sectors = tracks = 0; average = 0; }
	static RadialGrid createGrid(const std::vector<Vector3> & projected, int numTracks = 10, int numSectors = 40 );

	void rotate( int steps = 1 );
	void align();

	// DEBUG:
	void draw(QPainter & painter);
};

struct Sphere{
    typedef SurfaceMesh::SurfaceMeshModel::Vector3 Vector3;

    SurfaceMesh::SurfaceMeshModel * geometry;
    SurfaceMesh::SurfaceMeshModel::Scalar radius;
    Vector3 center;
	int numPoints;
	Vector3 majorAxis;
	Spherelib::RadialGrid grid;

    Sphere(int resolution = 4, Vector3 center = Vector3(0,0,0), double radius = 1.0 );

	std::vector<Vector3> rays();
	std::vector<double> values();

	// Experimental / debug
	void fillRandom();	
	void fillPattern();
	void smoothValues( int iterations = 1 );
	void normalizeValues();

	// Grid related
	void assignLocalFrame( int tracks, int sectors );
	Spherelib::RadialGrid createGrid( const std::vector<double> & fromValues, int tracks, int sectors );
};

}
