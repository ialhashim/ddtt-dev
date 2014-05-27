#pragma once

#include <Eigen/Core>

class Particle
{
public:
    Particle(const Eigen::Vector3d& pos) : pos(pos), measure(0.0) {}

    size_t id, correspondence;
    Eigen::Vector3d pos;
    double measure;
};
