#pragma once

#include <QThread>
#include <QElapsedTimer>

#include "CorrespondenceGenerator.h"

class CorrespondenceSearch : public QThread
{
    Q_OBJECT
public:
	CorrespondenceSearch(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB, const Paths & paths);

    PropertyMap property;
    QVector<RenderObject::Base*> debug;

	Structure::ShapeGraph *shapeA, *shapeB;
	
	Paths paths;
	std::vector<double> pathScores;
	int bestCorrespondence;

public slots:
    void run();
    void increaseProgress();

signals:
    void done();
    void pathComputed();
};
