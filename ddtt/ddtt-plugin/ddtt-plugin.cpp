#include <QFileDialog>

#include "ddtt-plugin.h"
#include "interfaces/ModePluginDockWidget.h"
#include "ui_ddtt_widget.h"
#include "ddtt_widget.h"

#include "StructureGraph.h"

#include "Corresponder.h"

#define BBOX_WIDTH(box) (box.max().x()-box.min().x())
#define PADDING_FACTOR 1.0

ddtt_widget * w = NULL;

void ddtt::create()
{
	if(!widget)
	{
		ModePluginDockWidget * dockwidget = new ModePluginDockWidget("TopoBlender", mainWindow());

		// Setup widget and signals
        w = new ddtt_widget();

		connect(w->ui->testButton, SIGNAL(clicked()), SLOT(execute()));
		connect(w->ui->clearButton, SIGNAL(clicked()), SLOT(clear()));

		dockwidget->setWidget(w);
		dockwidget->setWindowTitle(w->windowTitle());
		mainWindow()->addDockWidget(Qt::RightDockWidgetArea, dockwidget);
		this->widget = w;

		drawArea()->setShortcut(QGLViewer::DRAW_AXIS, Qt::Key_A);
		drawArea()->setShortcut(QGLViewer::DRAW_GRID, Qt::Key_G);
	}
}

void ddtt::decorate()
{
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
		//graphs[g]->draw();
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

void ddtt::execute()
{
	if( graphs.isEmpty() )
	{
		loadModels(QFileDialog::getOpenFileNames(0, "Open Model", 
			mainWindow()->settings()->getString("lastUsedDirectory"), "Model Files (*.xml)"));
	}

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

void ddtt::clear()
{
	drawArea()->clear();
	graphs.clear();
	drawArea()->updateGL();
}
