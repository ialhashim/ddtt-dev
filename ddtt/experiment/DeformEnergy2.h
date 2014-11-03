#pragma once

#include "ShapeGraph.h"
#include "RenderObjectExt.h"

typedef QPair<QString, QString> PairParts;
typedef QVector<PairParts> ShapeEdges;

class DeformEnergy2
{
public:
    DeformEnergy2(Structure::ShapeGraph * shapeA, Structure::ShapeGraph * shapeB,
                 const QVector<QStringList> &a_landmarks,
                 const QVector<QStringList> &b_landmarks,
                 bool debugging);
    Structure::ShapeGraph *a, *b;

	static Array2D_Vector4d sideCoordinates;

	QVector<QStringList> mappedPartsA, mappedPartsB;
	QMap<QString, QStringList> mappingAB;

	PropertyMap energyTerms;
	double total_energy;

	// Energy terms:
	double computeAngles(const ShapeEdges & B, const ShapeEdges & B_missed);
	double computeEdges(const ShapeEdges & B, const ShapeEdges & B_missed);
	double computeContext(const QStringList & A);
	double computeSymmetry(const QStringList & A);

	// Regularization:
	double E_regularizer(const QVector<QStringList> & correspondedSet);

	// Utility:
	QPair< QVector<PairParts>, QVector<PairParts> > correspondedEdges(const PairParts & partsA);

	// Debug:
    bool debugging;
    QVector<RenderObject::Base*> debug;
};
