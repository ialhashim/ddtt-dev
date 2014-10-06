#include "deform.h"
#include "deform-widget.h"
#include "ui_deform-widget.h"
#include "RenderObjectExt.h"

void deform::create()
{
    if(widget) return;

    // Viewer
    double worldRadius = mesh()->bbox().diagonal().norm();
    {
        drawArea()->setAxisIsDrawn(true);
        drawArea()->camera()->setType(qglviewer::Camera::PERSPECTIVE);

        drawArea()->camera()->setUpVector(qglviewer::Vec(0,0,1));
        drawArea()->camera()->setPosition(qglviewer::Vec(worldRadius * 2,-worldRadius * 2, worldRadius * 1.5));
        drawArea()->camera()->lookAt(qglviewer::Vec());
        drawArea()->camera()->setSceneRadius( worldRadius );
        drawArea()->camera()->showEntireScene();
    }

    // Rendering
    {
        drawArea()->setRenderer(mesh(), "Flat Wire");
    }

    ModePluginDockWidget * dockwidget = new ModePluginDockWidget("Deform", mainWindow());
    DeformWidget * dw = new DeformWidget();
    widget = dw;

    dockwidget->setWidget( widget );
    mainWindow()->addDockWidget(Qt::RightDockWidgetArea, dockwidget);

    // UI
    connect(dw->ui->createAnchor, &QPushButton::released, [&]{
		DeformWidget * dw = (DeformWidget *)widget;
        auto vid = Vertex( dw->ui->vertexID->value() );
        Vector3 p = mesh()->vertex_coordinates()[vid];

        auto handleRadius = worldRadius * 0.2;

        QSharedPointer<DeformHandle> handle( new DeformHandle(p, handleRadius) );

		this->connect(handle.data(), SIGNAL(manipulated()), SLOT(apply_deformation()));

		drawArea()->setManipulatedFrame( handle.data() );

		handles << handle;
        drawArea()->update();
    });
}

void deform::apply_deformation()
{

}

void deform::decorate()
{
	if (handles.isEmpty()) return;
	
	double worldRadius = mesh()->bbox().diagonal().norm();

	// Draw handles
	auto handleRadius = worldRadius * 0.2;
	starlab::FrameSoup fs(handleRadius);
	starlab::PointSoup ps (20);
	for (auto & handle : handles)
	{
		auto handlepos = handle->position();
		auto p = Vector3(handlepos[0], handlepos[1], handlepos[2]);
		ps.addPoint(p);
		fs.addFrame(Vector3(Vector3::UnitX()), Vector3(Vector3::UnitY()), Vector3(Vector3::UnitZ()), p);
	}
	ps.draw();
	fs.draw();
}

void deform::drawWithNames()
{

}

bool deform::postSelection(const QPoint &)
{
    return true;
}

bool deform::keyPressEvent(QKeyEvent *)
{
    return false;
}

bool deform::mouseMoveEvent(QMouseEvent *)
{
    return false;
}

bool deform::mousePressEvent(QMouseEvent * event)
{
    if (event->modifiers() & Qt::SHIFT)
    {
        drawArea()->select(event->pos());
        return true;
    }

    return false;
}
