#include "experiment.h"
#include "experiment-widget.h"
#include "ui_experiment-widget.h"

#include "SurfaceMeshHelper.h"
#include "RenderObjectExt.h"

#include "GraphCorresponder.h"
#include "Scheduler.h"
#include "TopoBlender.h"

void experiment::create()
{
    // Prepare UI
    if( !widget )
    {
        // Setup viewer
        {
            drawArea()->setAxisIsDrawn(true);
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
    }

    graphs << new Structure::ShapeGraph("C:/Temp/dataset/ChairBasic1/SimpleChair1.xml");
    graphs << new Structure::ShapeGraph("C:/Temp/dataset/ChairBasic2/shortChair01.xml");

    GraphCorresponder gcorr( graphs.front(), graphs.back() );

    QSharedPointer<Scheduler> scheduler ( new Scheduler );
    QSharedPointer<TopoBlender> blender = QSharedPointer<TopoBlender>( new TopoBlender(&gcorr, scheduler.data()) );
    //scheduler->executeAll();
}

void experiment::decorate()
{
    int startX = -0.5;
    for(auto g : graphs)
    {
        glPushMatrix();
        glTranslated(startX,0,0);

        g->draw();

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
		Vector3 p(pos[0], pos[1], pos[2]);
		if (!found) return false;

		for (auto g : graphs)
		{
			QMap<double, QString> distMap;
			for (auto & n : g->nodes){
				auto nodeBox = n->bbox();
				Vector3 projection = p;
				for (int i = 0; i < 3; i++){ 
					projection[i] = std::min(p[i], nodeBox.max()[i]);
					projection[i] = std::max(p[i], nodeBox.min()[i]);
				}
				distMap[(p - projection).norm()] = n->id;
			}
			auto closeNode = g->getNode(distMap.values().front());
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