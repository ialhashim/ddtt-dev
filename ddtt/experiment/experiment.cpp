#include "experiment.h"
#include "experiment-widget.h"
#include "ui_experiment-widget.h"

#include "SurfaceMeshHelper.h"
#include "RenderObjectExt.h"

#include "GraphCorresponder.h"
#include "Scheduler.h"
#include "TopoBlender.h"

#include "myglobals.h"

static QVector<QColor> colors = rndColors2(100);

#include "Deformer.h"
void experiment::doCorrespond()
{
	int num_solver_iterations = ((ExperimentWidget*)widget)->ui->numIterations->value();

	Deformer d( graphs.front(), graphs.back(), num_solver_iterations );

	for(auto debug : d.debug)drawArea()->addRenderObject( debug );
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
			starlab::SphereSoup ss;
			//glDisable(GL_LIGHTING);

			int r = 0;
			for (auto & landmark : g->landmarks)
			{
				//glBegin(GL_POINTS);
				//auto c = colors[r++];
				//glColor3d(c.redF(), c.greenF(), c.blueF());
				//glVertex3dv(landmark.data());
				//glEnd();

				ss.addSphere(landmark, 0.02, colors[r++]);
			}
			ss.draw();
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
	if (event->modifiers() & Qt::SHIFT)
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
			size_t minidx = std::min_element(dists.begin(),dists.end()) - dists.begin();
			auto closeNode = g->getNode(distMap.keys().at(minidx));

			auto coord = closeNode->approxCoordinates(p);
			auto closestPoint = closeNode->position(coord);

			Structure::Landmark landmark(g->landmarks.size(), closestPoint);
			landmark.u = coord[0];
			landmark.v = coord[1];
			landmark.partid = closeNode->id;

			g->landmarks.push_back( landmark );
		}

		drawArea()->update();
	}

	return false;
}