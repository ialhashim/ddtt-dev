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

#include "Serializable.h"

template<typename Vector3>
struct Particle : public Serializable
{
	typedef double Scalar;

    explicit Particle(const Vector3& pos = Vector3(0,0,0)) : pos(pos), measure(0.0), weight(1), 
		alpha(1.0), direction(Vector3(0,0,1)), flag(NONE), avgDiameter(0), segment(0), isMedial(false)
	{
		id = -1; // an invalid ID
		correspondence = -1;
	}

	size_t id, correspondence;
	uint64_t morton;
	ParticleFlags flag;
	int segment;
    Vector3 pos, direction, axis;
    Scalar measure, weight, alpha;
	Scalar avgDiameter, flat;
	bool isMedial;

	// Serialization:
	void serialize(QDataStream& os) const {
		os << id << correspondence << morton << ((qint32)flag) << segment;
		os << pos << direction << axis;
		os << measure << weight << alpha;
		os << avgDiameter << flat;
		os << isMedial;
	}	
	void deserialize(QDataStream& is) {
		int flagItem;
		is >> id >> correspondence >> morton >> flagItem >> segment;
		is >> pos >> direction >> axis;
		is >> measure >> weight >> alpha;
		is >> avgDiameter >> flat;
		is >> isMedial;
		flag = (ParticleFlags)flagItem;
	}
};
