#pragma warning(disable:4267)

#pragma once

#include "Hungarian.h"
#undef max
#undef min

#include <QVector>
#include <QVariant>
#include "RenderObjectExt.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace Structure{ struct Graph; }

#include "DeformationPath.h"

struct ShapeCorresponder{
	Structure::Graph * source;
	Structure::Graph * target;
	PropertyMap property;

	std::vector<DeformationPath> paths;

	ShapeCorresponder(Structure::Graph * g1, Structure::Graph * g2, bool isFullSet = true);

	QVector<RenderObject::Base *> debug;
};
