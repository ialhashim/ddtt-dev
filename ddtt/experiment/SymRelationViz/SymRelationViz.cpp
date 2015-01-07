#include "SymRelationViz.h"
#include "SymRelationViz-widget.h"
#include "ui_SymRelationViz-widget.h"
#include <QMessageBox>
#include <QFileDialog>

SymRelationVizWidget * pw;

void SymRelationViz::changeScale(double ds)
{
	m_scene.m_modelScale = ds;
	m_scene.layout();
	m_scene.buildModelDislayList();
}

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
	pw->ui->scaleSpinBox->setValue(m_scene.m_modelScale);
	//m_scene.m_modelScale
	connect(pw->ui->scaleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(changeScale(double)));

    connect(pw->ui->loadButton, &QPushButton::released, [&]{
        //auto files = QFileDialog::getOpenFileNames(mainWindow(), "");
		QString dirname = QFileDialog::getExistingDirectory(mainWindow(), tr("Open Directory"),
			".",
			QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
		
		if (dirname.isNull() || dirname.isEmpty()) return;
		m_scene.clearScene();
		m_scene.loadScene(dirname);
		m_scene.layout();
		m_scene.buildModelDislayList();
    });
}

void SymRelationViz::decorate()
{
	m_scene.draw();
    // Draw stuff here:
	//glEnable(GL_LIGHTING);
	//glEnable(GL_POLYGON_OFFSET_FILL);
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
