#pragma once

#include <Qhull.h>
#include <QhullVertexSet.h>
#include <QhullFacetList.h>

template<typename Vector3>
struct ConvexHull{
	ConvexHull() { volume = area = 0; }
	ConvexHull( const std::vector<Vector3> & in_points, const char * options = "FA Qt") : volume(0), area(0), center(Vector3(0,0,0))
	{
		try{
			// Compute convex hull of points:
			orgQhull::Qhull q("", 3, (int)in_points.size(), in_points[0].data(), options);

			for(orgQhull::QhullFacet & f : q.facetList()){
				std::vector<Vector3> face;
				for(auto v : f.vertices()){
					auto p = v.point();
					face.push_back( Vector3(p[0],p[1],p[2]) );
					center += face.back();
				}
				if(f.isTopOrient()) std::reverse(face.begin(),face.end());
				faces.push_back( face );
			}

			center /= (3*faces.size());

			// Get volume and area:
			volume = q.volume();
			area = q.area();
		} catch (const std::exception &e) {
			qDebug() << e.what();
		}
	}
	std::vector< std::vector< Vector3 > > faces;
	Vector3 center;
	double volume, area;

	ConvexHull<Vector3> merged( const ConvexHull<Vector3> & otherHull ){
		std::vector< Vector3 > both_points;
		for(auto & f : faces) for(auto & v : f) both_points.push_back(v);
		for(auto & f : otherHull.faces) for(auto & v : f) both_points.push_back(v);
		return ConvexHull<Vector3>( both_points );
	}
};
