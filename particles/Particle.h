#pragma once

#include <vector>
#include <stdint.h>
#include <Eigen/Core>

template<typename Vector3>
struct Particle
{
	typedef double Scalar;

    Particle(const Vector3& pos) : pos(pos), measure(0.0), weight(1), 
		alpha(1.0), direction(Vector3(0,0,1)), flag(0), avgDiameter(0) {}

	size_t id, correspondence;
	uint64_t morton;
	int flag;
    Vector3 pos, direction, relativePos;
    Scalar measure;
	Scalar weight;
	Scalar alpha;
	Scalar avgDiameter;
};
