#pragma once

#include <QOpenGLShaderProgram>
#include <QWidget>
#include <QScrollArea>
#include <QScrollBar>
#include <QFileDialog>
#include <QElapsedTimer>

#include <vector>
#include <stdint.h>
#include <Eigen/Core>

enum ParticleFlags{ NONE, FLOOR, UNPROCESSED, VIZ_WEIGHT };

template<typename Vector3>
struct Particle
{
	typedef double Scalar;

    explicit Particle(const Vector3& pos) : pos(pos), measure(0.0), weight(1), 
		alpha(1.0), direction(Vector3(0,0,1)), flag(NONE), avgDiameter(0), segment(0) 
	{
		id = -1; // an invalid ID
	}

	size_t id, correspondence;
	uint64_t morton;
	ParticleFlags flag;
	int segment;
    Vector3 pos, direction, relativePos, axis;
    Scalar measure;
	Scalar weight;
	Scalar alpha;
	float avgDiameter, flat;
};
