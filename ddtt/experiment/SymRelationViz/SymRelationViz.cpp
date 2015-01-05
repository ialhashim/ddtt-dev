#include "SymRelationViz.h"
#include "SymRelationViz-widget.h"
#include "ui_SymRelationViz-widget.h"
#include <QMessageBox>
#include <QFileDialog>

SymRelationVizWidget * pw;

void SymRelationViz::create()
{
	// Prepare UI
	if (widget) return;

	// Setup viewer
	{
        drawArea()->setAxisIsDrawn(true);
		drawArea()->camera()->setType(qglviewer::Camera::ORTHOGRAPHIC);

		double worldRadius = 1;
		drawArea()->camera()->setUpVector(qglviewer::Vec(0, 0, 1));
		drawArea()->camera()->setPosition(qglviewer::Vec(-0.36, -2.2, 1.3));
		auto center = qglviewer::Vec(0.5, 0, 0.5);
		drawArea()->setSceneCenter(center);
		drawArea()->camera()->lookAt(center);
		drawArea()->camera()->setSceneRadius(worldRadius);
		drawArea()->camera()->showEntireScene();
	}

    ModePluginDockWidget * dockwidget = new ModePluginDockWidget("SymRelationViz", mainWindow());
	pw = new SymRelationVizWidget();
	widget = pw;

	dockwidget->setWidget(widget);
	mainWindow()->addDockWidget(Qt::RightDockWidgetArea, dockwidget);

	// UI:
    connect(pw->ui->doButton, &QPushButton::released, [&]{
        QMessageBox::information(0,"hello","hi");
	});

    connect(pw->ui->loadButton, &QPushButton::released, [&]{
        auto files = QFileDialog::getOpenFileNames(mainWindow(), "");
    });
}

void SymRelationViz::decorate()
{
    // Draw stuff here:
}

bool SymRelationViz::keyPressEvent(QKeyEvent * event)
{
	if (event->key() == Qt::Key_P)
	{
	}

	drawArea()->update();
	return false;
}

bool SymRelationViz::mouseMoveEvent(QMouseEvent *)
{
	return false;
}

bool SymRelationViz::mousePressEvent(QMouseEvent *event)
{
	if (event->modifiers())
	{
		
		drawArea()->update();
	}

	return false;
}
