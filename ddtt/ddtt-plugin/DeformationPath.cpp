#include <QApplication>
#include "DeformationPath.h"
#include "SynthesisManager.h"

SynthesisManager * smanager = NULL;
Q_DECLARE_METATYPE(SynthesisManager *)

DeformationPath::DeformationPath(){
    gcorr = NULL;
    weight = 0.0;
    idx = i = si = 0;
}

void DeformationPath::execute()
{
	scheduler = QSharedPointer<Scheduler>( new Scheduler );
	blender = QSharedPointer<TopoBlender>( new TopoBlender( gcorr, scheduler.data() ) );
	
	scheduler->executeAll();
	for(auto & g : scheduler->allGraphs) g->moveBottomCenterToOrigin( true );
	property["isReady"].setValue( true );
}

void DeformationPath::renderSamples()
{
	if(!gcorr) return;

	qApp->setOverrideCursor(Qt::WaitCursor);

	scheduler = QSharedPointer<Scheduler>( new Scheduler );
	blender = QSharedPointer<TopoBlender>( new TopoBlender( gcorr, scheduler.data() ) );

	// Surface sampling
	smanager = new SynthesisManager(gcorr, scheduler.data(), blender.data(), 4000);
	smanager->genSynData();

	// Execute path
	scheduler->executeAll();
	for(auto & g : scheduler->allGraphs) g->moveBottomCenterToOrigin( true );

	if( scheduler->allGraphs.size() )
	{
		int N = 6;

		for(int i = 0; i < N; i++)
		{
			double t = double(i) / (N-1);
			QString path = QFileInfo(gcorr->sg->property["name"].toString()).absolutePath() + "/";
			QString name = "output";
			QString filename = path + QString("%1-%2-%3").arg(name).arg(i).arg(t);

			//int idx = t * (scheduler->allGraphs.size()-1);
			//smanager->renderGraph(*scheduler->allGraphs[idx], filename, false, 6);
		}
	}

	property["synthManager"].setValue( smanager );

	qApp->restoreOverrideCursor();
	QCursor::setPos(QCursor::pos());

	property["isReady"].setValue( true );
}

void DeformationPath::renderProxies()
{
	if(!gcorr) return;

	qApp->setOverrideCursor(Qt::WaitCursor);

	scheduler = QSharedPointer<Scheduler>( new Scheduler );
	blender = QSharedPointer<TopoBlender>( new TopoBlender( gcorr, scheduler.data() ) );

	// Surface sampling via proxies
	smanager = new SynthesisManager(gcorr, scheduler.data(), blender.data());
	smanager->makeProxies(60,20);

	// Execute path
	scheduler->executeAll();
	for(auto & g : scheduler->allGraphs) g->moveBottomCenterToOrigin( true );

	property["synthManager"].setValue( smanager );

	qApp->restoreOverrideCursor();
	QCursor::setPos(QCursor::pos());

	property["isReady"].setValue( true );
}
