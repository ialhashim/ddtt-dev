#pragma warning(disable:4267)

#include <QFileDialog>

#include "ddtt-plugin.h"
#include "interfaces/ModePluginDockWidget.h"
#include "ui_ddtt_widget.h"
#include "ddtt_widget.h"

#include "StructureGraph.h"
#include "SynthesisManager.h"
#include "GraphExplorer.h"

#include "Corresponder.h"

#include "ShapeCorresponder.h"

#include "DeformScene.h"

#define BBOX_WIDTH(box) (box.max().x()-box.min().x())
#define PADDING_FACTOR 1.0

ddtt_widget * w = NULL;
ShapeCorresponder * sc = NULL;

void ddtt::create()
{
	if(!widget)
	{
		ModePluginDockWidget * dockwidget = new ModePluginDockWidget("TopoBlender", mainWindow());

		// Setup widget and signals
        w = new ddtt_widget();

		connect(w->ui->testButton, SIGNAL(clicked()), SLOT(execute()));
		connect(w->ui->clearButton, SIGNAL(clicked()), SLOT(clear()));

		connect(w->ui->loadGraphs, SIGNAL(clicked()), SLOT(loadGraphs()));
		connect(w->ui->correspondButton, SIGNAL(clicked()), SLOT(correspond()));

		dockwidget->setWidget(w);
		dockwidget->setWindowTitle(w->windowTitle());
		mainWindow()->addDockWidget(Qt::RightDockWidgetArea, dockwidget);
		this->widget = w;

		drawArea()->setShortcut(QGLViewer::DRAW_AXIS, Qt::Key_A);
		drawArea()->setShortcut(QGLViewer::DRAW_GRID, Qt::Key_G);

		w->ui->loadGraphs->click();
		//w->ui->correspondButton->click();
	}
}

SynthesisManager * sm = NULL;

void ddtt::decorate()
{
	// Cylinder experiment
	{
		// Test proxies
		if(graphs.size())
		{
			if(!sm)
			{
				GraphCorresponder * g = new GraphCorresponder(graphs.front(), graphs.back());
				Scheduler * s = new Scheduler;
				TopoBlender * t = new TopoBlender(g,s);

				sm = new SynthesisManager ( g, s, t );
				sm->makeProxies(20, 3);
			}

			if(sm)sm->drawWithProxies( sm->scheduler->activeGraph );
		}
	}

	double startX = bigbox.min().x();

	for(int g = 0; g < (int) graphs.size(); g++)
	{
		// Place and draw graph
		glPushMatrix();

		Eigen::AlignedBox3d curbox = graphs[g]->bbox();

		double curwidth = (curbox.max().x() - curbox.min().x());
		double deltaX = curwidth * 0.5;

		double padding = 0;
		if(g > 0) padding = curwidth * PADDING_FACTOR;

		double posX = startX + deltaX + padding;

		if(graphs.size() < 2) posX = 0;

		glTranslated(posX, 0, 0);

		// store for later use
		graphs[g]->property["posX"] = posX;
		graphs[g]->draw();
		//drawBBox( curbox );

		glPopMatrix();

		startX += curwidth + padding;
	}
}

void ddtt::setSceneBounds()
{
	if(!graphs.size())
	{
		drawArea()->setSceneRadius(2);
		drawArea()->setSceneCenter(qglviewer::Vec(0,0,0));
		drawArea()->setSceneBoundingBox(qglviewer::Vec(-1,-1,-1), qglviewer::Vec(1,1,1));
		drawArea()->camera()->setPosition(qglviewer::Vec(-1,-3,2));
		drawArea()->showEntireScene();
		drawArea()->updateGL();
		return;
	}

	// Set scene bounds
	bigbox = graphs.front()->bbox();
	double deltaX = BBOX_WIDTH(bigbox);
	bigbox.translate( Vector3(deltaX * 0.5, 0, 0) ); // start from zero

	for(int i = 1; i < (int)graphs.size(); i++)
	{
		Eigen::AlignedBox3d curbox = graphs[i]->bbox();

		double curWidth = BBOX_WIDTH(curbox);
		double padding = curWidth * PADDING_FACTOR;

		curbox.translate( Vector3(deltaX + (0.5 * curWidth) + padding, 0, 0) );
		bigbox = bigbox.merged( Eigen::AlignedBox3d(curbox) );

		deltaX += BBOX_WIDTH(curbox) + padding; 
	}

	// Move to center
	bigbox.translate( Vector3(-bigbox.center().x(), 0, 0) );

	// Setup camera
	{
		Vector3 a = bigbox.min();
		Vector3 b = bigbox.max();

		qglviewer::Vec vecA(a.x(), a.y(), a.z());
		qglviewer::Vec vecB(b.x(), b.y(), b.z());

		drawArea()->camera()->setUpVector(qglviewer::Vec(0,0,1));
		drawArea()->camera()->setPosition(qglviewer::Vec(-2,-3,1));
		drawArea()->camera()->lookAt(qglviewer::Vec(0,0,0));

		drawArea()->setSceneCenter((vecA + vecB) * 0.5);
		drawArea()->setSceneBoundingBox(vecA, vecB);
		drawArea()->showEntireScene();
		drawArea()->updateGL();
	}
}

void ddtt::loadModels(QStringList fileNames)
{
	if(fileNames.isEmpty()) return;

	foreach(QString file, fileNames){
		graphs.push_back( new Structure::Graph(file) );
	}

	QFileInfo fileInfo(fileNames.back());
	mainWindow()->settings()->set( "lastUsedDirectory", fileInfo.absolutePath() );

	setSceneBounds();
}

void ddtt::loadGraphs()
{
	if( graphs.size() < 2 )
	{
		loadModels(QFileDialog::getOpenFileNames(0, "Open Model", 
			mainWindow()->settings()->getString("lastUsedDirectory"), "Model Files (*.xml)"));
	}
}

void ddtt::clear()
{
	drawArea()->clear();
	graphs.clear();
	drawArea()->updateGL();
}

void ddtt::execute()
{
	loadGraphs();

	int nidx = w->ui->nodeIndex->value();
	int metric = w->ui->metricBox->currentIndex();
	int splits = w->ui->splits->value();

	PropertyMap prop;
	prop["forceBest"] = w->ui->bestAssign->isChecked();

	Corresponder c(graphs.front(), graphs.back(), nidx, (Corresponder::Metric)metric, splits);
	c.compute( prop );

	drawArea()->clear();

	for(auto r : c.debug)
		drawArea()->addRenderObject(r);

	drawArea()->updateGL();
}

void ddtt::correspond()
{
	loadGraphs();

	QElapsedTimer timer; timer.start();

	// Make first has smaller number of nodes
	if(graphs.front()->nodes.size() > graphs.back()->nodes.size())
		std::swap( graphs.front(), graphs.back() );

	sc = new ShapeCorresponder( graphs.front(), graphs.back(), !w->ui->bestAssign->isChecked() );

	drawArea()->clear();
	for(auto r : sc->debug) drawArea()->addRenderObject(r);
	drawArea()->update();

	// Stats
	QStringList message;
	message << QString("Number of paths ( %1 ).").arg( sc->property["pathsCount"].toInt() );
	message << QString("Prepare time ( %1 ms ).").arg( sc->property["prepareTime"].toInt() );
	message << QString("Compute time ( %1 ms ).").arg( sc->property["computeTime"].toInt() );
	message << QString("Evaluate time ( %1 ms ).").arg( sc->property["evaluateTime"].toInt() );
	message << QString("Total correspondence search time ( %1 ms ).").arg(timer.elapsed());
	mainWindow()->setStatusBarMessage( message.join("\n") );

	// View correspondences
	DeformScene * ds = new DeformScene;

	int limit = 100;

	int i = 0;
	
	for(auto & path : sc->paths)
	{
		path.idx = i++;
		ds->addDeformationPath( &path );

		if(limit > 0 && i > limit) break; 
	}

	ds->pack();
}

bool ddtt::keyPressEvent(QKeyEvent* event)
{
	if(event->key() == Qt::Key_I)
	{
		GraphExplorer * ge = new GraphExplorer;
		ge->update( graphs.back() );
		ge->show();

		return true;
	}

	return false;
}
