#include "experiment.h"
#include "experiment-widget.h"
#include "ui_experiment-widget.h"
#include <QFileDialog>

#include "SurfaceMeshHelper.h"
#include "RenderObjectExt.h"

#include "GraphCorresponder.h"
#include "Scheduler.h"
#include "TopoBlender.h"

#include "myglobals.h"

static QVector<QColor> colors = rndColors2(100);

#include "CorrespondenceSearch.h"
void experiment::doCorrespondSearch()
{
	auto search = new CorrespondenceSearch(graphs.front(), graphs.back(), CorrespondenceGenerator(graphs.front(), graphs.back()).generate());
	
	search->start();
	//mainWindow()->setStatusBarMessage(QString("Search done in (%1 ms)").arg(search->property["allSearchTime"].toInt()));
}

#include "DeformEnergy.h"
void experiment::doCorrespond()
{
    QElapsedTimer timer; timer.start();

    QVector<QStringList> landmarks_front;
    QVector<QStringList> landmarks_back;

    for(auto landmark : graphs.front()->landmarks){
        QStringList m;
        for(auto l : landmark) m << l.partid;
        landmarks_front << m;
    }

    for(auto landmark : graphs.back()->landmarks){
        QStringList m;
        for(auto l : landmark) m << l.partid;
        landmarks_back << m;
    }

    DeformEnergy d( graphs.front(), graphs.back(), landmarks_front, landmarks_back );

    for (auto debug : d.debug)drawArea()->addRenderObject(debug);

    mainWindow()->setStatusBarMessage(QString("%1 ms").arg(timer.elapsed()));
}

#include "Deformer.h"
void experiment::doCorrespond2()
{
	QElapsedTimer timer; timer.start();

	int num_solver_iterations = ((ExperimentWidget*)widget)->ui->numIterations->value();

	Deformer d( graphs.front(), graphs.back(), num_solver_iterations );

	for(auto debug : d.debug)drawArea()->addRenderObject( debug );

	mainWindow()->setStatusBarMessage(QString("%1 ms").arg(timer.elapsed()));
}

void experiment::create()
{
    // Prepare UI
	if (widget) return;

	graphs << new Structure::ShapeGraph("C:/Temp/dataset/ChairBasic1/SimpleChair1.xml");
	graphs << new Structure::ShapeGraph("C:/Temp/dataset/ChairBasic2/shortChair01.xml");

	graphs.front()->setColorAll(Qt::blue);
	graphs.back()->setColorAll(Qt::green);

	graphs.front()->loadLandmarks("front.landmarks");
	graphs.back()->loadLandmarks("back.landmarks");

	//GraphCorresponder gcorr( graphs.front(), graphs.back() );
	//QSharedPointer<Scheduler> scheduler ( new Scheduler );
	//QSharedPointer<TopoBlender> blender = QSharedPointer<TopoBlender>( new TopoBlender(&gcorr, scheduler.data()) );
	//scheduler->executeAll();

    // Setup viewer
    {
        //drawArea()->setAxisIsDrawn(true);
        drawArea()->camera()->setType(qglviewer::Camera::PERSPECTIVE);

        double worldRadius = 1;
        drawArea()->camera()->setUpVector(qglviewer::Vec(0,0,1));
        drawArea()->camera()->setPosition(qglviewer::Vec(2,-2,1.5));
        drawArea()->camera()->lookAt(qglviewer::Vec());
        drawArea()->camera()->setSceneRadius( worldRadius );
        drawArea()->camera()->showEntireScene();
    }

    ModePluginDockWidget * dockwidget = new ModePluginDockWidget("Particles", mainWindow());
    auto * pw = new ExperimentWidget();
    widget = pw;

    dockwidget->setWidget( widget );
    mainWindow()->addDockWidget(Qt::RightDockWidgetArea, dockwidget);
    
	// UI:
	connect(pw->ui->saveLandmarks, &QPushButton::released, [=]{
		graphs.front()->saveLandmarks("front.landmarks");
		graphs.back()->saveLandmarks("back.landmarks");
		mainWindow()->setStatusBarMessage("Landmarks saved.");
	});
	connect(pw->ui->loadLandmarks, &QPushButton::released, [=]{
		graphs.front()->loadLandmarks("front.landmarks");
		graphs.back()->loadLandmarks("back.landmarks");
		drawArea()->update();
	});
	connect(pw->ui->clearLandmarks, &QPushButton::released, [=]{
		for (auto g : graphs) g->landmarks.clear();
		drawArea()->update();
	});
	connect(pw->ui->executeButton, &QPushButton::released, [=]{
		drawArea()->clear();
		this->doCorrespond();
		drawArea()->update();
	});
	connect(pw->ui->clearShapes, &QPushButton::released, [=]{
		graphs.clear();
		drawArea()->clear();
		drawArea()->update();
	});
	connect(pw->ui->loadShapes, &QPushButton::released, [=]{
		QString filename = QFileDialog::getOpenFileName(mainWindow(), tr("Load Shape"), "", tr("Shape File (*.xml)"));
		if (!filename.size()) return;
		graphs << new Structure::ShapeGraph(filename);
		drawArea()->update();
	});
	connect(pw->ui->searchBest, &QPushButton::released, [=]{
		drawArea()->clear();
		this->doCorrespondSearch();
		drawArea()->update();
	});
}

void experiment::decorate()
{
    int startX = 0;
    for(auto g : graphs)
    {
        glPushMatrix();
        glTranslated(startX,0,0);

        g->draw();
		g->property["startX"].setValue(startX);

		if (((ExperimentWidget*)widget)->ui->isShowLandmarks->isChecked())
		{
			//starlab::SphereSoup ss;
			glDisable(GL_LIGHTING);

			int r = 0;
			for (auto & l : g->landmarks)
			{
				glBegin(GL_POINTS);
				auto c = colors[r++];
				glColor3d(c.redF(), c.greenF(), c.blueF());
				for(auto & landmark : l) glVertex3dv(landmark.data());
				glEnd();

				//for (auto & landmark : l) ss.addSphere(landmark, 0.02, colors[r]);
				r++;
			}
			//ss.draw();
		}

        glPopMatrix();
        startX += 1;
    }
}

bool experiment::keyPressEvent(QKeyEvent *)
{
    return false;
}

bool experiment::mouseMoveEvent(QMouseEvent *)
{
    return false;
}

bool experiment::mousePressEvent(QMouseEvent *event)
{
	if (event->modifiers())
	{
		bool found = false;
		auto pos = drawArea()->camera()->pointUnderPixel(event->pos(), found);
		Vector3 worldPos(pos[0], pos[1], pos[2]);
		if (!found) return false;

		for (auto g : graphs)
		{
			Vector3 p = worldPos - Vector3(g->property["startX"].toDouble(), 0, 0);
			if (!g->cached_bbox().contains(p)) continue;

			QMap<QString, double> distMap;
			for (auto n : g->nodes){
				auto nodeBox = n->bbox(1.01);
				Vector3 projection = p;
				for (int i = 0; i < 3; i++){
					projection[i] = std::min(p[i], nodeBox.max()[i]);
					projection[i] = std::max(projection[i], nodeBox.min()[i]);
				}
				distMap[n->id] = (p - projection).norm();
			}
			auto dists = distMap.values();
			size_t minidx = std::min_element(dists.begin(), dists.end()) - dists.begin();
			auto closeNode = g->getNode(distMap.keys().at(minidx));

			auto coord = closeNode->approxCoordinates(p);
			auto closestPoint = closeNode->position(coord);

			Structure::Landmark landmark(g->landmarks.size(), closestPoint);
			landmark.u = coord[0];
			landmark.v = coord[1];
			landmark.partid = closeNode->id;

			if (event->modifiers() & Qt::SHIFT)	g->landmarks.push_back(Structure::Landmarks());
			if (g->landmarks.size()) g->landmarks.back().push_back(landmark);
		}

		drawArea()->update();
	}
	
	return false;
}
