#include "experiment.h"
#include "experiment-widget.h"
#include "ui_experiment-widget.h"
#include <QFileDialog>
#include <QListWidget>
#include <QMatrix4x4>

#include "SurfaceMeshHelper.h"
#include "RenderObjectExt.h"

#include "GraphCorresponder.h"
#include "Scheduler.h"
#include "TopoBlender.h"

#include "myglobals.h"

static QVector<QColor> colors = rndColors2(100);
ExperimentWidget * pw = NULL;

#include "DeformEnergy.h"
#include "DeformEnergy2.h"

#include "Deformer.h"
#include "CorrespondenceSearch.h"
CorrespondenceSearch * search = NULL;

#include "EnergyGuidedDeformation.h"
void experiment::doEnergySearch()
{
	QElapsedTimer timer; timer.start();

	auto shapeA = new Structure::ShapeGraph(*graphs.front());
	auto shapeB = new Structure::ShapeGraph(*graphs.back());

	QMatrix4x4 mat;

	if (pw->ui->isAnisotropy->isChecked())
	{
		auto bboxA = shapeA->bbox();
		auto bboxB = shapeB->bbox();

		Vector3 s = bboxA.diagonal().array() / bboxB.diagonal().array();
		mat.scale(s.x(), s.y(), s.z());
		shapeB->transform(mat, true);

		// Show
		//graphs.removeLast();
		//graphs.push_back(shapeB);
	}

	QVector<QStringList> landmarks_front;
	QVector<QStringList> landmarks_back;

	for (auto landmark : graphs.front()->landmarks){
		QStringList m;
		for (auto l : landmark) m << l.partid;
		landmarks_front << m;
	}

	for (auto landmark : graphs.back()->landmarks){
		QStringList m;
		for (auto l : landmark) m << l.partid;
		landmarks_back << m;
	}

	EnergyGuidedDeformation egd(shapeA, shapeB, landmarks_front, landmarks_back, true);

	for (auto debug : egd.debug) drawArea()->addRenderObject(debug);
	for (auto debug : shapeA->debug) drawArea()->addRenderObject(debug);
	for (auto debug : shapeB->debug) drawArea()->addRenderObject(debug);

	// Show deformed
	graphs.clear();
	graphs << shapeA << shapeB;
	//for (auto g : graphs) g->property["showCtrlPts"].setValue(true);

	mainWindow()->setStatusBarMessage(QString("%1 ms").arg(timer.elapsed()));
}

void experiment::doCorrespondSearch()
{
	pw->ui->searchBest->setEnabled(false);

	auto shapeA = new Structure::ShapeGraph(*graphs.front());
	auto shapeB = new Structure::ShapeGraph(*graphs.back());
	
	QMatrix4x4 mat;

	if (pw->ui->isAnisotropy->isChecked())
	{
		auto bboxA = shapeA->bbox();
		auto bboxB = shapeB->bbox();

		Vector3 s = bboxA.diagonal().array() / bboxB.diagonal().array();
		mat.scale(s.x(), s.y(), s.z());
		shapeB->transform(mat,true);

		// Show
		//graphs.removeLast();
		//graphs.push_back(shapeB);
	}

	search = new CorrespondenceSearch(shapeA, shapeB, CorrespondenceGenerator(shapeA, shapeB).generate(), pw->ui->isOtherEnergy->isChecked());

	connect(search, SIGNAL(done()), SLOT(postCorrespond()));

	search->start();
}

void experiment::showCorrespond(int idx)
{
	if (graphs.empty() || !search || search->paths.empty() || idx > search->paths.size()-1) return;

	drawArea()->clear();

	for (auto n : graphs.front()->nodes){
		graphs.front()->setColorFor(n->id, QColor(255, 255, 255, 10));
		n->vis_property["meshSolid"].setValue(false);
	}
	for (auto n : graphs.back()->nodes){
		graphs.back()->setColorFor(n->id, QColor(255, 255, 255, 10));
		n->vis_property["meshSolid"].setValue(false);
	}

	auto nidA = search->paths[idx].first;
	auto nidB = search->paths[idx].second;

	QMap<QString, QStringList> mapback;

	for (size_t i = 0; i < nidA.size(); i++)
	{
		auto nidsA = nidA[i], nidsB = nidB[i];
		QColor color = colors[i];
		for (auto nid : nidsA) {
			graphs.front()->setColorFor(nid, color);
			graphs.front()->getNode(nid)->vis_property["meshSolid"].setValue(true);
		}
		for (auto nid : nidsB) {
			graphs.back()->setColorFor(nid, color);
			graphs.back()->getNode(nid)->vis_property["meshSolid"].setValue(true);

			if (mapback.contains(nid)){
				for (auto nid : mapback[nid]) {
					graphs.front()->setColorFor(nid, color);
					graphs.front()->getNode(nid)->vis_property["meshSolid"].setValue(true);
				}
			}
			else
				mapback[nid] = nidsA;
		}
	}

	// Apply for debugging:
	{
		// anisotropy:		
		auto shapeA = new Structure::ShapeGraph(*graphs.front());
		auto shapeB = new Structure::ShapeGraph(*graphs.back());
		{
			QMatrix4x4 mat;

			if (pw->ui->isAnisotropy->isChecked())
			{
				auto bboxA = shapeA->bbox();
				auto bboxB = shapeB->bbox();

				Vector3 s = bboxA.diagonal().array() / bboxB.diagonal().array();
				mat.scale(s.x(), s.y(), s.z());
				shapeB->transform(mat, true);
			}
		}

		if (!pw->ui->isOtherEnergy->isChecked())
		{
			DeformEnergy2 d(shapeA, shapeB, nidA, nidB, pw->ui->isVisualize->isChecked());
			for (auto debug : d.debug) drawArea()->addRenderObject(debug);
		}
		else
		{
			DeformEnergy d(shapeA, shapeB, nidA, nidB, pw->ui->isVisualize->isChecked());
			for (auto debug : d.debug) drawArea()->addRenderObject(debug);
		}
	}
}

void experiment::postCorrespond()
{
	pw->ui->searchBest->setEnabled(true);

	mainWindow()->setStatusBarMessage(QString("Search done in (%1 ms)").arg(search->property["allSearchTime"].toInt()));

	bool isVisualize = pw->ui->isVisualize->isChecked();

	if (pw->ui->isShowParts->isChecked())
		showCorrespond( search->bestCorrespondence );

	pw->ui->pathsList->clear();

	// List scores
	{
		QSet<double> scoreSet;

		for (size_t pi = 0; pi < search->pathScores.size(); pi++)
		{
			//if (scoreSet.contains(search->pathScores[pi])) continue; // avoid duplicated scores
			scoreSet.insert(search->pathScores[pi]);

			//Arg1: the number, Arg2: how many 0 you want?, Arg3: i don't know but only 10 can take negative numbers
			QString number;
			number.sprintf("%09.3f ", search->pathScores[pi]);

			for (auto key : search->pathDetails[pi].keys()) number += QString(",%1=%2").arg(key.left(4)).arg(search->pathDetails[pi][key].toDouble());

			auto item = new QListWidgetItem(number);
			item->setData(Qt::UserRole, pi);
			pw->ui->pathsList->addItem(item);
		}
	}

	drawArea()->update();
}

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

    DeformEnergy2 d( graphs.front(), graphs.back(), landmarks_front, landmarks_back, true );

	for (auto debug : d.debug)drawArea()->addRenderObject(debug);

	// Display score
	{
		pw->ui->pathsList->clear();
		QSet<double> scoreSet;

		QString number;
		number.sprintf("%09.3f ", d.total_energy);
		for (auto key : d.energyTerms.keys()) number += QString(",%1=%2").arg(key.left(4)).arg(d.energyTerms[key].toDouble());

		auto item = new QListWidgetItem(number);
		item->setData(Qt::UserRole, 0);
		pw->ui->pathsList->addItem(item);
	}

    mainWindow()->setStatusBarMessage(QString("%1 ms").arg(timer.elapsed()));
}

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

	// Load last used shapes:
	QSettings settingsFile(QSettings::IniFormat, QSettings::UserScope, "GrUVi", "experiment");
	QString shapeLeft = settingsFile.value("shapeLeft", "C:/Temp/dataset/ChairBasic1/SimpleChair1.xml").toString();
	QString shapeRight = settingsFile.value("shapeRight", "C:/Temp/dataset/ChairBasic2/shortChair01.xml").toString();
	settingsFile.sync();

	graphs << new Structure::ShapeGraph(shapeLeft);
	graphs << new Structure::ShapeGraph(shapeRight);

	//graphs.front()->loadLandmarks("front.landmarks");
	//graphs.back()->loadLandmarks("back.landmarks");

	graphs.front()->setColorAll(Qt::blue);
	graphs.back()->setColorAll(Qt::green);

    // Setup viewer
    {
        //drawArea()->setAxisIsDrawn(true);
        drawArea()->camera()->setType(qglviewer::Camera::ORTHOGRAPHIC);

        double worldRadius = 1;
        drawArea()->camera()->setUpVector(qglviewer::Vec(0,0,1));
        drawArea()->camera()->setPosition(qglviewer::Vec(-0.36,-2.2,1.3));
		auto center = qglviewer::Vec(0.5, 0, 0.5);
		drawArea()->setSceneCenter(center);
        drawArea()->camera()->lookAt(center);
        drawArea()->camera()->setSceneRadius( worldRadius );
        drawArea()->camera()->showEntireScene();
    }

    ModePluginDockWidget * dockwidget = new ModePluginDockWidget("Particles", mainWindow());
    pw = new ExperimentWidget();
    widget = pw;

    dockwidget->setWidget( widget );
    mainWindow()->addDockWidget(Qt::RightDockWidgetArea, dockwidget);
    
	// UI:
	connect(pw->ui->saveLandmarks, &QPushButton::released, [&]{
		graphs.front()->saveLandmarks("front.landmarks");
		graphs.back()->saveLandmarks("back.landmarks");
		mainWindow()->setStatusBarMessage("Landmarks saved.");
	});
	connect(pw->ui->loadLandmarks, &QPushButton::released, [&]{
		graphs.front()->loadLandmarks("front.landmarks");
		graphs.back()->loadLandmarks("back.landmarks");
		drawArea()->update();
	});
	connect(pw->ui->clearLandmarks, &QPushButton::released, [&]{
		for (auto g : graphs) g->landmarks.clear();
		drawArea()->update();
	});
	connect(pw->ui->executeButton, &QPushButton::released, [&]{
		drawArea()->clear();
		this->doCorrespond();
		drawArea()->update();
	});
	connect(pw->ui->clearShapes, &QPushButton::released, [&]{
		graphs.clear();
		drawArea()->clear();
		pw->ui->pathsList->clear();
		drawArea()->update();
	});
	connect(pw->ui->swapButton, &QPushButton::released, [&]{
		std::swap(graphs.front(), graphs.back());
		drawArea()->update();
	});
	connect(pw->ui->loadShapes, &QPushButton::released, [&]{
		QString filename = QFileDialog::getOpenFileName(mainWindow(), tr("Load Shape"), "", tr("Shape File (*.xml)"));
		if (!filename.size()) return;
		graphs << new Structure::ShapeGraph(filename);
		drawArea()->update();

		// Save last used shapes:
		if (graphs.size() == 2){
			QSettings settingsFile(QSettings::IniFormat, QSettings::UserScope, "GrUVi", "experiment");
			settingsFile.setValue("shapeLeft", graphs.front()->property["name"].toString());
			settingsFile.setValue("shapeRight", graphs.back()->property["name"].toString());
			settingsFile.sync();
		}
	});
	connect(pw->ui->resetShapes, &QPushButton::released, [&]{
		drawArea()->clear();
		graphs.clear();
		QSettings settingsFile(QSettings::IniFormat, QSettings::UserScope, "GrUVi", "experiment");
		QString shapeLeft = settingsFile.value("shapeLeft").toString();
		QString shapeRight = settingsFile.value("shapeRight").toString();
		graphs << new Structure::ShapeGraph(shapeLeft);
		graphs << new Structure::ShapeGraph(shapeRight);
		drawArea()->update();
	});
	connect(pw->ui->searchBest, &QPushButton::released, [&]{
		drawArea()->clear();
		this->doCorrespondSearch();
		drawArea()->update();
	});
	connect(pw->ui->pathsList, &QListWidget::itemSelectionChanged, [&]{
		if (search){
			int pid = pw->ui->pathsList->currentItem()->data(Qt::UserRole).toInt();
			showCorrespond(pid);
		}
		drawArea()->update();
	});

	connect(pw->ui->doEnergyGuided, &QPushButton::released, [&]{
		drawArea()->clear();
		this->doEnergySearch();
		drawArea()->update();
	});

	this->isReady = true;
}

void experiment::decorate()
{
	if (!isReady) return;

    double startX = 0;
	for (size_t i = 0; i < graphs.size(); i++)
    {
		auto & g = graphs[i];

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

		// Spacing between models:
		if (i + 1 < graphs.size()){
			auto & g2 = graphs[i + 1];
			Eigen::AlignedBox3d bbox = g2->bbox();
			if (!g2->property.contains("width"))
			{
				for (auto n : g2->nodes){
					auto m = g2->getMesh(n->id);
					if (!m) g2->property["width"].setValue(bbox.sizes().x());
					m->updateBoundingBox();
					bbox.extend(m->bbox());
				}
				g2->property["width"].setValue(bbox.sizes().x());

				drawArea()->setSceneRadius(g->property["width"].toDouble() * 4);
			}
			double g1_width = g->property["width"].toDouble();
			double g2_width = g2->property["width"].toDouble();

			startX += (g1_width * 0.3) + g2_width + 0.05;
		}else{
			Eigen::AlignedBox3d bbox = g->bbox();
			for (auto n : g->nodes){
				auto m = g->getMesh(n->id);
				if (!m) g->property["width"].setValue(bbox.sizes().x());
				m->updateBoundingBox();
				bbox.extend(m->bbox());
			}
			g->property["width"].setValue(bbox.sizes().x());
		}
    }
}

bool experiment::keyPressEvent(QKeyEvent * event)
{
	if (event->key() == Qt::Key_P) 
	{
		drawArea()->setStateFileName("cameraSettings.xml");
		drawArea()->saveStateToFile();
	}

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
			auto graphBox = g->bbox();

			// Expand for very flat
			double diagonal_len = graphBox.diagonal().norm();
			double s = diagonal_len * 0.01;
			graphBox.extend(graphBox.max() + Vector3(s, s, s));
			graphBox.extend(graphBox.min() - Vector3(s, s, s));

			if (!graphBox.contains(p)) continue;

			QMap<QString, double> distMap;
			for (auto n : g->nodes){
				auto nodeBox = n->bbox();

				// Expand for very flat
				if (nodeBox.diagonal().minCoeff() == 0){
					double diagonal_len = nodeBox.diagonal().norm();
					double s = diagonal_len * 0.01;
					nodeBox.extend(nodeBox.max() + Vector3(s, s, s));
					nodeBox.extend(nodeBox.min() - Vector3(s, s, s));
				}

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
