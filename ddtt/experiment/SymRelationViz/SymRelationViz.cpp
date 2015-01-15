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
void SymRelationViz::changePtSize(double ds)
{
	m_scene.m_ptSize = ds;
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
	pw->ui->ptSizeSpinBox->setValue(m_scene.m_ptSize);
	connect(pw->ui->scaleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(changeScale(double)));
	connect(pw->ui->ptSizeSpinBox, SIGNAL(valueChanged(double)), this, SLOT(changePtSize(double)));
	connect(pw->ui->modelTreeWidget, SIGNAL(itemClicked(QTreeWidgetItem *, int)),
		this, SLOT(changeModelOrPart(QTreeWidgetItem *, int)));

	connect(pw->ui->loadButton, &QPushButton::released, [&]
	{
		const QString DEFAULT_DIR_KEY("default_dir");
		QSettings MySettings;

		QString dirname = QFileDialog::getExistingDirectory(mainWindow(), tr("Open Directory"),
			MySettings.value(DEFAULT_DIR_KEY).toString(),
			QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

		if (dirname.isNull() || dirname.isEmpty()) return;
		QDir dir(dirname);
		MySettings.setValue(DEFAULT_DIR_KEY, dir.absolutePath());

		drawArea()->setAxisIsDrawn(false);
		m_scene.clearScene();
		m_scene.loadScene(dirname);
		buildModelTree();
		m_scene.layout();
		m_scene.buildModelDislayList();
	});
}
void SymRelationViz::buildModelTree()
{
	pw->ui->modelTreeWidget->setColumnCount(2);
	for (int i = 0; i < m_scene.m_modelNum; ++i)
	{
		MeshModel* mm = m_scene.m_modelList[i];
		QTreeWidgetItem *treeItem = new QTreeWidgetItem(pw->ui->modelTreeWidget);
		treeItem->setText(0, mm->m_name);
		for (int j = 0; j < mm->m_parts.size(); ++j)
		{
			QVector<int> &tmp = mm->m_matches[0];
			if ( j+1 >= tmp.size()) continue;
			QTreeWidgetItem *treeItemPart = new QTreeWidgetItem();
			treeItemPart->setText(0, QString::number(j));
			QString str;
			for (int k = 0; k < mm->m_matches.size(); ++k)
			{
				str.append(QString::number(mm->m_matches[k][j]));
				str.append(" ");
			}
			treeItemPart->setText(1, str);
			treeItem->addChild(treeItemPart);
		}
	}
}

void SymRelationViz::changeModelOrPart(QTreeWidgetItem *item, int colNo)
{
	int modelNo(-1), partNo(-1);
	if (0 == item->childCount()) // a part is selected
	{
		QTreeWidgetItem *parent = item->parent();
		modelNo = pw->ui->modelTreeWidget->indexOfTopLevelItem(parent);
		partNo = parent->indexOfChild(item);
	}
	else // a model is selected
	{
		modelNo = pw->ui->modelTreeWidget->indexOfTopLevelItem(item);		
	}
	m_scene.setCurrentModelAndPart(modelNo, partNo);
	m_scene.buildModelDislayList();
	drawArea()->update();
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
