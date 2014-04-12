#pragma once
#include "GraphCorresponder.h"
#include "Task.h"
#include "Scheduler.h"
#include "TopoBlender.h"

typedef QMap<QString,QVariant> PropertyMap;
typedef QPair< QVector<QString>, QVector<QString> > Pairing;
typedef QVector< Pairing > VectorPairings;
typedef QVector< VectorPairings > Assignments;
typedef QVector< QPair<QString,QString> > VectorPairStrings;

class DeformationPath
{
public:
    DeformationPath();

    double weight;
    VectorPairStrings pairsDebug;
    QVector<Pairing> pairs;
    QVector<double> errors;

    QMap<QString, QColor> scolors;
    QMap<QString, QColor> tcolors;

    int idx;
    int i;
	int si;
    PropertyMap property;

    GraphCorresponder * gcorr;
    QSharedPointer<Scheduler> scheduler;
    QSharedPointer<TopoBlender> blender;

	void normalizeErrors();

    static double minWeight( const std::vector<DeformationPath> & paths );
    static double maxWeight( const std::vector<DeformationPath> & paths );
	static void normalize( std::vector<DeformationPath> & paths );

	void execute();
};

static inline bool DeformationPathCompare (const DeformationPath & i, const DeformationPath & j) { return (i.weight < j.weight); }

