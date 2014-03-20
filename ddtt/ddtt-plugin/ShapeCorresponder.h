#pragma once
#include "RenderObjectExt.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace Structure{ struct Graph; }
typedef QMap<QString,QVariant> PropertyMap;

struct ShapeCorresponder{
	Structure::Graph * source;
	Structure::Graph * target;
	PropertyMap property;

	ShapeCorresponder(Structure::Graph * g1, Structure::Graph * g2);

	QVector<RenderObject::Base *> debug;
};
