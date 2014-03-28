#pragma warning(disable:4267)

#pragma once
#include <QVector>
#include <QVariant>
#include "RenderObjectExt.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace Structure{ struct Graph; }

typedef QMap<QString,QVariant> PropertyMap;
typedef QPair<QStringList,QStringList> Pairing;

#include "Hungarian.h"
#undef max
#undef min


#include "GraphCorresponder.h"
#include "Task.h"
#include "Scheduler.h"
#include "TopoBlender.h"

typedef QVector< QPair<QString,QString> > VectorPairStrings;
struct DeformationPath{
	double weight;
	VectorPairStrings pairsDebug;
	QVector<Pairing> pairs;
	QVector<double> errors;

	QMap<QString, QColor> scolors;
	QMap<QString, QColor> tcolors;

	int idx;
    int i;

	GraphCorresponder * gcorr;
	QSharedPointer<Scheduler> scheduler;
	QSharedPointer<TopoBlender> blender;
    DeformationPath(){ gcorr = NULL; weight = 0.0; idx = i = 0; }
};
static inline bool DeformationPathCompare (const DeformationPath & i, const DeformationPath & j) { return (i.weight < j.weight); }


struct ShapeCorresponder{
	Structure::Graph * source;
	Structure::Graph * target;
	PropertyMap property;

	std::vector<DeformationPath> paths;

	ShapeCorresponder(Structure::Graph * g1, Structure::Graph * g2);

	QVector<Pairing> findPairing( mat m, QVector<Pairing> fixedPairs = QVector<Pairing>() );

	QVector<RenderObject::Base *> debug;
};
