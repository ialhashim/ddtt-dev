#pragma once

#include <QThread>
#include <QElapsedTimer>

#include "PartCorresponder.h"
#include "CorrespondenceGenerator.h"

typedef QVector< QPair<size_t,size_t> > PartsCorrespondence;

class CorrespondenceSearch : public QThread
{
    Q_OBJECT
public:
    CorrespondenceSearch(CorrespondenceGenerator * generator);

    PropertyMap property;
    QVector<RenderObject::Base*> debug;

	ParticleMesh *sA, *sB;
    CorrespondenceGenerator * generator;

	PartsCorrespondence bestCorrespondence;
	void applyCorrespondence( const PartsCorrespondence & correspondence );

public slots:
    void run();
    void increaseProgress();

signals:
    void done();
    void pathComputed();
};
