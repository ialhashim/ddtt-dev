#include "Raytracing.h"

#define USE_EMBREE

#ifdef USE_EMBREE

#include "embree2/rtcore.h"
#include "embree2/rtcore_ray.h"
#pragma comment(lib, "embree.lib")
namespace embree{
	struct Vertex   { float x, y, z, a; };
	struct Triangle { int v0, v1, v2; };
};
using namespace embree;

#else
#include "Octree.h"
#endif

#include <QElapsedTimer>

template<typename Vector3>
raytracing::Raytracing<Vector3>::Raytracing(SurfaceMesh::SurfaceMeshModel * mesh, std::vector<Vector3> & rayOrigins, 
											std::vector<Vector3> & rayDirections): time(-1)
{
	QElapsedTimer timer; timer.start();

	int rayCount = (int)rayOrigins.size();
	hits = std::vector<RayHit>( rayCount );

    #ifdef USE_EMBREE
    // Using Intel's accelerated ray tracer
    {
		/* initialize ray tracing core */
		rtcInit( "" );

		/* create scene */
		RTCScene scene = rtcNewScene( RTC_SCENE_STATIC, RTC_INTERSECT1 );

		/* Triangle mesh loading */
		int numTriangles = mesh->n_faces();
		int numVertices = mesh->n_vertices();
		unsigned int geomID = rtcNewTriangleMesh(scene, RTC_GEOMETRY_STATIC, numTriangles, numVertices, 1);

		/* set vertex */
		Vertex* vertices = (Vertex*) rtcMapBuffer(scene, geomID, RTC_VERTEX_BUFFER);
		for(int i = 0; i < numVertices; i++)
		{
			Eigen::Vector3d p = mesh->vertex_coordinates()[SurfaceMesh::Vertex(i)];
			vertices[i].x = p.x();
			vertices[i].y = p.y();
			vertices[i].z = p.z();
		}
		rtcUnmapBuffer(scene, geomID, RTC_VERTEX_BUFFER);

		/* set triangles */
		Triangle* triangles = (Triangle*) rtcMapBuffer(scene, geomID, RTC_INDEX_BUFFER);
		for(int i = 0; i < numTriangles; i++)
		{
			std::vector<int> mesh_triangle;
			for(auto v : mesh->vertices(SurfaceMesh::Face(i))) mesh_triangle.push_back(v.idx());

			triangles[i].v0 = mesh_triangle[0];
			triangles[i].v1 = mesh_triangle[1];
			triangles[i].v2 = mesh_triangle[2];
		}
		rtcUnmapBuffer(scene,geomID,RTC_INDEX_BUFFER);

		/* commit changes to scene */
		rtcCommit (scene);

		#pragma omp parallel for
		for(int i = 0; i < rayCount; i++)
		{
			Vector3 p(rayOrigins[i]), d(rayDirections[i]);

			/* initialize ray */
			RTCRay ray;
			ray.org[0] = p[0]; ray.org[1] = p[1]; ray.org[2] = p[2];
			ray.dir[0] = d[0]; ray.dir[1] = d[1]; ray.dir[2] = d[2];
			ray.tnear = 0.0f;
			ray.tfar = std::numeric_limits<float>::infinity();
			ray.geomID = RTC_INVALID_GEOMETRY_ID;
			ray.primID = RTC_INVALID_GEOMETRY_ID;
			ray.mask = -1;
			ray.time = 0;

			/* intersect ray with scene */
			rtcIntersect( scene, ray);

			hits[i] = RayHit( ray.tfar, ray.primID );
		}

		rtcDeleteScene( scene );
    }
    #else
    // Using an Octree
    {
        // Build octree
        Octree octree( mesh );

        #pragma omp parallel for
        for(int i = 0; i < rayCount; i++)
        {
            int findex = -1;
            Eigen::Vector3d isect = octree.closestIntersectionPoint(Ray(rayOrigins[i].cast<double>(), rayDirections[i].cast<double>()), &findex, true);
            double dist = (isect - rayOrigins[i]).norm();

            hits[i] = RayHit( dist, findex );
        }
    }
    #endif

	time = timer.elapsed();
}

// Explicit template instantiation
template class raytracing::Raytracing<Eigen::Vector3f>;
