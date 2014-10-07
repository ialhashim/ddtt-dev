#pragma once

#include "SurfaceMeshHelper.h"
#include "qglviewer/manipulatedFrame.h"
using namespace qglviewer;
using namespace SurfaceMesh;

class DeformHandle : public ManipulatedFrame
{
public:
    DeformHandle(const Vector3 & start, double Radius){
        this->startPos = start;
        this->radius = Radius;
        this->setPosition(start.x(), start.y(), start.z());

		this->isActive = false;
    }

	Vector3 transformed(const Vector3 & originalPos, double gaussianWeight = 1.0){
        Vector3 d = originalPos - startPos;
        Vec delta(d[0], d[1], d[2]);
        Vec rotatedDelta = this->rotation() * delta;
        Vec r = this->position() + rotatedDelta;
        Vector3 newPos(r.x, r.y, r.z);

        double alpha = 1 - gaussianFunction(((originalPos - startPos).norm() / radius));
		alpha *= gaussianWeight;

        return (originalPos * (alpha)) + (newPos * (1-alpha));
    }

	Vector3 pos(){ auto p = this->position(); return Vector3(p[0], p[1], p[2]); }

public:
	QVector<size_t> element_id;
	QVector<Vector3> element_orig_pos;
	QVector<size_t> constraint_id;
	bool isActive;

private:
	Vector3 startPos;
    double radius;

    double inline gaussianFunction(double x, double mu = 0.0, double sigma = 2){
        //double a = 1.0 / (sigma * sqrt(2 * M_PI));
        double b = mu;
        double c = sigma;
        return exp( - (pow(x - b, 2) / (2 * pow(c, 2)) ) );
    }
};
