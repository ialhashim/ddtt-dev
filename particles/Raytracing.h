#pragma once

#include <Eigen/Core>
#include "SurfaceMeshModel.h"

namespace raytracing{

struct RayHit{
    bool isHit;
    double distance;
    int faceid;
    RayHit(double distance = DBL_MAX, int faceid = -1) : distance(distance), faceid(faceid) {
        isHit = (faceid >= 0) ? true : false;
    }
};

template<typename Vector3>
class Raytracing
{
public:
    Raytracing( SurfaceMesh::SurfaceMeshModel * mesh, std::vector<Vector3> & rayOrigins, std::vector<Vector3> & rayDirections );
    std::vector<RayHit> hits;
    int time;
};

}
