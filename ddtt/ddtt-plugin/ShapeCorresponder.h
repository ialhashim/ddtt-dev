#pragma once
#include <QVector>
#include <QVariant>
#include "RenderObjectExt.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace Structure{ struct Graph; }

typedef QMap<QString,QVariant> PropertyMap;
typedef QPair<QStringList,QStringList> Pairing;

struct ShapeCorresponder{
	Structure::Graph * source;
	Structure::Graph * target;
	PropertyMap property;

	ShapeCorresponder(Structure::Graph * g1, Structure::Graph * g2);

	QVector<Pairing> findPairing( QVector<Pairing> fixedPairs = QVector<Pairing>() );

	QVector<RenderObject::Base *> debug;
};
