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
#include "CorrespondenceSearch.h"
CorrespondenceSearch * mysearch = NULL;
#include "Deformer.h"
#include "Synthesizer.h"

#include "EvaluateCorrespondence.h"

#include "Viewer.h"

#include "BatchProcess.h"

Energy::GuidedDeformation * egd = NULL;
Viewer * v = NULL;
Q_DECLARE_METATYPE(Structure::ShapeGraph*);
Q_DECLARE_METATYPE(Energy::SearchNode*);
experiment * exprmnt = NULL;
Energy::SearchNode * selected_path;

QString pathToHtml(Energy::SearchNode & p)
{
	QStringList html;

	for (auto a : p.assignments)
	{
		html << QString("<span class=source>%1</span>").arg(a.first.join("-"));
		html << "<br/>";
		html << QString("<span class=target>%1</span>").arg(a.second.join("-"));
		if(p.assignments.size() > 1 && a != p.assignments.back()) html << "<br/>";
	}

	html << "<br/>";
	html << "<span class=cost>" + QString::number(p.cost) + "</span>";
	html << "<span class=children>(" + QString::number(p.num_children) + ")</span>";

	//html << "<span>" + p.current.join("-") + "</span>";
	//html << "<span>fixed:" + p.fixed.join("-") + "</span><br/>";
	//html << "<span>unassigned:" + p.unassigned.join("-") + "</span><br/>";
	return html.join("");
};

void experiment::setSearchPath(Energy::SearchNode * path)
{
	// Show deformed
	graphs.clear();
	graphs << path->shapeA << path->shapeB;

	// Grey out
	for (auto n : graphs.front()->nodes){
		graphs.front()->setColorFor(n->id, QColor(255, 255, 255, 10));
		n->vis_property["meshSolid"].setValue(false);
		if (n->type() == Structure::SHEET) ((Structure::Sheet*)n)->surface.quads.clear();
	}
	for (auto n : graphs.back()->nodes){
		graphs.back()->setColorFor(n->id, QColor(255, 255, 255, 10));
		n->vis_property["meshSolid"].setValue(false);
		if (n->type() == Structure::SHEET) ((Structure::Sheet*)n)->surface.quads.clear();
	}

	// Assign colors based on target
	int ci = 0;
	for (auto & relation : graphs.back()->relations)
	{
		QColor color = colors[ci++];
		for (auto nid : relation.parts)
		{
			graphs.back()->setColorFor(nid, color);
			graphs.back()->getNode(nid)->vis_property["meshSolid"].setValue(true);
		}
	}

	// Color matching source
	for (auto spart : selected_path->mapping.keys())
	{
		auto tpart = selected_path->mapping[spart];
		if (tpart == Structure::null_part) continue;

		auto color = graphs.back()->getNode(tpart)->vis_property["meshColor"].value<QColor>();

		graphs.front()->setColorFor(spart, color);
		graphs.front()->getNode(spart)->vis_property["meshSolid"].setValue(true);
	}

	// For visualization
	if (!graphs.front()->animation.empty())
	{
		graphs.front()->setAllControlPoints(graphs.front()->animation.front());
		encodeGeometry();
		graphs.front()->setAllControlPoints(graphs.front()->animation.back());
	}
}

void experiment::doEnergySearch()
{
	QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));

	//for (auto g : graphs) g->property["showCtrlPts"].setValue(true);
	for (auto g : graphs) g->property["showMeshes"].setValue(false);

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

	QElapsedTimer timer; timer.start();

	egd = new Energy::GuidedDeformation;

	// Create a search path
	Energy::Assignments assignments;
	for (size_t i = 0; i < landmarks_front.size(); i++) assignments << qMakePair(landmarks_front[i], landmarks_back[i]);

	QVector<Energy::SearchNode> search_roots;
	Energy::SearchNode path(shapeA, shapeB, QStringList(), assignments);
	search_roots << path;

	// Explore path
	egd->searchAll(graphs.front(), graphs.back(), search_roots);

	auto timeElapsed = timer.elapsed();

	// Visualize search graph
	{
		exprmnt = this;

		v = new Viewer;
		v->show();

		connect(v->wv, &QWebView::loadFinished, [&](){
			v->addCSS(".source{color:red} .target{color:blue} .cost{color:gray} .children{color:gray}");

			auto & path = *egd->searchTrees.front().begin();
			int idx = v->addNode(pathToHtml(path));
			Energy::SearchNode* path_ptr = &path;
			v->nodeProperties[idx]["path"].setValue(path_ptr);

			v->updateGraph();

			connect(v, &Viewer::nodeSelected, [&](int nid){
				auto node = v->nodeProperties[nid]["path"].value<Energy::SearchNode*>();
				if (!node) return;

				egd->applySearchPath( egd->getEntirePath(node) );

				selected_path = node;
				exprmnt->setSearchPath(selected_path);

				auto shape = node->shapeA;
				if (!shape) return;
				node->property["nid"].setValue(nid);

				// Expand children
				if (!v->nodeProperties[nid]["expanded"].toBool())
				{
					// Expand children
					for (auto child : egd->childrenOf(selected_path))
					{
						int child_idx = v->addNode(pathToHtml(*child));
						v->addEdge(nid, child_idx);
						v->nodeProperties[child_idx]["path"].setValue(child);
					}

					v->updateGraph();
					v->nodeProperties[nid]["expanded"].setValue(true);
				}

				exprmnt->drawArea()->update();
			});
		});

		connect(v, &Viewer::goingToExpand, [&](){
			auto toExpand = egd->getEntirePath(selected_path);

			int nid = toExpand.front()->property["nid"].value<int>();

			for (auto element : toExpand)
			{
				int child_nid = v->addNode(pathToHtml(*element));
				v->addEdge(nid, child_nid);
				v->nodeProperties[child_nid]["path"].setValue(element);

				nid = child_nid;
			}

			v->updateGraph();
		});
	}

	// Select least cost path
	auto all_solutions = egd->solutions();
	QList < QPair<double, Energy::SearchNode*> > solutions;
	for (auto s : all_solutions) solutions << qMakePair(s->cost, s);
	qSort(solutions.begin(), solutions.end());
	
	auto entire_path = egd->getEntirePath(solutions.front().second);
	egd->applySearchPath( entire_path );
	selected_path = entire_path.back();
	setSearchPath( selected_path );

	QApplication::restoreOverrideCursor();

	double cost = EvaluateCorrespondence::evaluate(graphs.front());
	mainWindow()->setStatusBarMessage(QString("%1 ms - cost = %2 - solutions %3").arg(timeElapsed).arg(cost).arg(solutions.size()));
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
		shapeB->transform(mat, true);

		// Show
		//graphs.removeLast();
		//graphs.push_back(shapeB);
	}

	mysearch = new CorrespondenceSearch(shapeA, shapeB, CorrespondenceGenerator(shapeA, shapeB).generate(), pw->ui->isOtherEnergy->isChecked());

	connect(mysearch, SIGNAL(done()), SLOT(postCorrespond()));

	mysearch->start();
}

void experiment::showCorrespond(int idx)
{
	if (graphs.empty() || !mysearch || mysearch->paths.empty() || idx > mysearch->paths.size() - 1) return;

	drawArea()->clear();

	for (auto n : graphs.front()->nodes){
		graphs.front()->setColorFor(n->id, QColor(255, 255, 255, 10));
		n->vis_property["meshSolid"].setValue(false);
	}
	for (auto n : graphs.back()->nodes){
		graphs.back()->setColorFor(n->id, QColor(255, 255, 255, 10));
		n->vis_property["meshSolid"].setValue(false);
	}

	auto nidA = mysearch->paths[idx].first;
	auto nidB = mysearch->paths[idx].second;

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

	mainWindow()->setStatusBarMessage(QString("Search done in (%1 ms)").arg(mysearch->property["allSearchTime"].toInt()));

    //bool isVisualize = pw->ui->isVisualize->isChecked();

	if (pw->ui->isShowParts->isChecked())
		showCorrespond(mysearch->bestCorrespondence);

	pw->ui->pathsList->clear();

	// List scores
	{
		QSet<double> scoreSet;

		for (size_t pi = 0; pi < mysearch->pathScores.size(); pi++)
		{
			//if (scoreSet.contains(search->pathScores[pi])) continue; // avoid duplicated scores
			scoreSet.insert(mysearch->pathScores[pi]);

			//Arg1: the number, Arg2: how many 0 you want?, Arg3: i don't know but only 10 can take negative numbers
			QString number;
			number.sprintf("%09.3f ", mysearch->pathScores[pi]);

			for (auto key : mysearch->pathDetails[pi].keys()) number += QString(",%1=%2").arg(key.left(4)).arg(mysearch->pathDetails[pi][key].toDouble());

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

	DeformEnergy2 d(graphs.front(), graphs.back(), landmarks_front, landmarks_back, true);

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

	Deformer d(graphs.front(), graphs.back(), num_solver_iterations);

	for (auto debug : d.debug)drawArea()->addRenderObject(debug);

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
		drawArea()->camera()->setUpVector(qglviewer::Vec(0, 0, 1));
		drawArea()->camera()->setPosition(qglviewer::Vec(-0.36, -2.2, 1.3));
		auto center = qglviewer::Vec(0.5, 0, 0.5);
		drawArea()->setSceneCenter(center);
		drawArea()->camera()->lookAt(center);
		drawArea()->camera()->setSceneRadius(worldRadius);
		drawArea()->camera()->showEntireScene();
	}

	ModePluginDockWidget * dockwidget = new ModePluginDockWidget("Particles", mainWindow());
	pw = new ExperimentWidget();
	widget = pw;

	dockwidget->setWidget(widget);
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

		// Save last used shapes:
		if (graphs.size() == 2){
			QSettings settingsFile(QSettings::IniFormat, QSettings::UserScope, "GrUVi", "experiment");
			settingsFile.setValue("shapeLeft", graphs.front()->property["name"].toString());
			settingsFile.setValue("shapeRight", graphs.back()->property["name"].toString());
			settingsFile.sync();
		}
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
		if (mysearch){
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
	connect(pw->ui->loadJobs, &QPushButton::released, [&]{
		QString filename = QFileDialog::getOpenFileName(mainWindow(), tr("Load Jobs"), "", tr("Jobs File (*.json)"));
		QTimer::singleShot(0, [=] { 
			BatchProcess * bp = new BatchProcess(filename);
			bp->start();
		});
	});
	connect(pw->ui->saveJob, &QPushButton::released, [&]{
		if (graphs.size() < 2 || graphs.front()->landmarks.isEmpty()) return;
		auto g1 = graphs.front(), g2 = graphs.back();

		QString filename = QFileDialog::getSaveFileName(mainWindow(), tr("Load Jobs"), "", tr("Jobs File (*.json)"));
		QVariantMap job;
		job["source"].setValue(g1->property["name"].toString());
		job["target"].setValue(g2->property["name"].toString());
		job["title"].setValue(g1->name() + "-" + g2->name());

		QVariantList assignments;
		for (int i = 0; i < g1->landmarks.size(); i++)
		{
			QVariantMap a;
			QStringList a1; for (auto l : g1->landmarks[i]) a1 << l.partid;
			QStringList a2; for (auto l : g2->landmarks[i]) a2 << l.partid;
			a["source"].setValue(a1);
			a["target"].setValue(a2);

			assignments.push_back(a);
		}
		job["assignments"].setValue(assignments);

		BatchProcess::appendJob(job, filename);
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
		glTranslated(startX, 0, 0);

		g->draw(drawArea());
		g->property["startX"].setValue(startX);
		g->property["posX"].setValue(startX);

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
				for (auto & landmark : l) glVertex3dv(landmark.data());
				glEnd();

				//for (auto & landmark : l) ss.addSphere(landmark, 0.02, colors[r]);
				r++;
			}
			//ss.draw();
		}

		for (auto debug : g->debug) debug->draw(*drawArea());

		// Debug deformation:
		if (g == graphs.front())
		{
			if (!g->animation_debug.isEmpty())
			{
				int idx = std::min(g->animation_debug.size()-1, g->animation_index);
				auto dbg = g->animation_debug[idx];
				for (auto d : dbg) d->draw(*drawArea());
			}

			if (g->property["showNames"].toBool())
			{
				for (auto n : g->nodes)
				{
					Vector3 diagonal = n->diagonal(), start = n->startPoint(), end = n->endPoint();
					Vector3 midPoint = 0.5 * (start + end);

					Vector3 orig_diagonal = n->property["orig_diagonal"].value<Vector3>();
					Vector3 orig_diagonal_start = midPoint - (orig_diagonal * 0.5);

					starlab::LineSegments ls(4);
					double v = (diagonal.normalized().dot(orig_diagonal.normalized()));
					ls.addLine(start, Vector3(n->startPoint() + diagonal), starlab::qtJetColor(v, -1, 1));
					ls.draw();
				}
			}

			// Spokes
			if (g->property["showSpokes"].toBool())
			{
				// Original spokes
				Array1D_Vector3 orig_spokes;
				for (auto l : g->edges)
					for (auto s : l->property["orig_spokes"].value<Array1D_Vector3>())
						orig_spokes.push_back(s);

				int si = 0;
				starlab::LineSegments ls;
				Array1D_Vector3 allspokes;
				for (auto l : g->edges)
				{
					auto samples1 = l->n1->property["samples_coords"].value<Array2D_Vector4d>();
					auto samples2 = l->n2->property["samples_coords"].value<Array2D_Vector4d>();

					// all combinations
					for (auto rowi : samples1){
						for (auto ci : rowi){
							for (auto rowj : samples2) {
								for (auto cj : rowj){
									auto a = l->n1->position(ci), b = l->n2->position(cj);
									Vector3 cur_spoke = a - b;
									double v = orig_spokes[si++].normalized().dot(cur_spoke.normalized());
									if (v < 0) v = 0;
									ls.addLine(a, b, starlab::qtJetColor(1.0 - v));
								}
							}
						}
					}
				}

				ls.draw();
			}
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
		}
		else{
			Eigen::AlignedBox3d bbox = g->bbox();
			for (auto n : g->nodes){
				auto m = g->getMesh(n->id);
				if (!m)
					g->property["width"].setValue(bbox.sizes().x());
				else
				{
					m->updateBoundingBox();
					bbox.extend(m->bbox());
				}
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

	if (event->key() == Qt::Key_N)
	{
		for (auto g : graphs) g->property["showNames"].setValue(!g->property["showNames"].toBool());

		for (auto g : graphs) g->property["showEdgeLines"].setValue(!g->property["showEdgeLines"].toBool());
		for (auto g : graphs) g->property["showEdges"].setValue(!g->property["showEdges"].toBool());
	}

	if (event->key() == Qt::Key_K)
	{
		for (auto g : graphs) g->property["showSpokes"].setValue(!g->property["showSpokes"].toBool());
	}

	if (event->key() == Qt::Key_M)
	{
		for (auto g : graphs)
		{
			for (auto n : g->nodes)
			{
				n->vis_property["meshColor"].setValue(n->vis_property["color"].value<QColor>());
				n->vis_property["meshSolid"].setValue(true);
			}

			g->property["showMeshes"].setValue(!g->property["showMeshes"].toBool());
			g->property["showNodes"].setValue(!g->property["showMeshes"].toBool());
		}
	}

	if (event->key() == Qt::Key_E)
	{
		encodeGeometry();
	}

	if (event->key() == Qt::Key_D)
	{
		decodeGeometry();
	}

	if (event->key() == Qt::Key_Space)
	{
		if (graphs.empty() || graphs.front()->animation.empty()) return false;

		QTimer * timer = new QTimer;
		connect(timer, &QTimer::timeout, [&]()
		{
			if (graphs.empty() || graphs.front()->animation.empty()) return;

			static bool isForward = true;
			static double t = 0;
			static int index = 0;

			static Structure::ShapeGraph * prev = NULL;
			if (prev != graphs.front()){
				t = 0;
				index = 0;
				prev = graphs.front();
			}

			auto source = graphs.front();

			double delta = pw->ui->speed->value();
			t += delta;

			if (t >= 1.0 || t <= 0.0)
			{
				t = (t >= 1.0) ? 0.0 : 1.0;

				index += 1;

				if (index >= source->animation.size() + 1)
				{
					t = 0;
					index = 0;
				}
			}

			auto ptsBefore = graphs.front()->animation[std::max(0, index - 1)];
			auto ptsAfter = graphs.front()->animation[std::min(index, graphs.front()->animation.size() - 1)];

			int activeNodes = graphs.front()->nodes.size();
			Array2D_Vector3 ptsCurrent;
			for (size_t i = 0; i < activeNodes; i++)
			{
				Array1D_Vector3 pts_cur;

				if (i < ptsBefore.size())
				{
					for (size_t j = 0; j < ptsBefore[i].size(); j++)
					{
						if (ptsBefore[i].size() != ptsAfter[i].size())
							ptsBefore[i] = ptsAfter[i];

						pts_cur.push_back(AlphaBlend(t, ptsBefore[i][j], ptsAfter[i][j]));
					}

					ptsCurrent.push_back(pts_cur);
				}
				else
				{
					auto n = graphs.front()->nodes[i];
					pts_cur.resize(n->numCtrlPnts(), Vector3(0,0,0));

					ptsCurrent.push_back(pts_cur);
				}
			}
			
			graphs.front()->setAllControlPoints(ptsCurrent);

			for (auto n : graphs.front()->nodes) if (n->type() == Structure::SHEET) ((Structure::Sheet*)n)->surface.quads.clear();

			if (graphs.front()->property["showMeshes"].toBool())
			{
				decodeGeometry();
			}

			graphs.front()->animation_index = index;

			drawArea()->update();

			//mainWindow()->setStatusBarMessage(QString("index [%1], t=%2").arg(index).arg(t));
		});
		timer->start(25);
	}

	drawArea()->update();
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

void experiment::encodeGeometry()
{
	for (auto g : graphs)
	{
		if (g->property["isGeometryEncoded"].toBool()) continue;

		for (auto n : g->nodes)
		{
			auto mesh = g->getMesh(n->id);
			if (!mesh) continue;
			std::vector < Vector3f > points, normals;
			auto mesh_points = mesh->vertex_coordinates(), mesh_normals = mesh->vertex_normals();
			for (auto v : mesh->vertices()){
				points.push_back(mesh_points[v].cast<float>());
				normals.push_back(mesh_normals[v].cast<float>());
			}

			QVector<ParameterCoord> encoding;
			if (n->type() == Structure::CURVE) encoding = Synthesizer::genPointCoordsCurve((Structure::Curve*)n, points, normals);
			if (n->type() == Structure::SHEET) encoding = Synthesizer::genPointCoordsSheet((Structure::Sheet*)n, points, normals);
			n->property["encoding"].setValue(encoding);
		}

		g->property["isGeometryEncoded"].setValue(true);

		break; // only encode source
	}
}

void experiment::decodeGeometry()
{
	if (!graphs.front()->property["isGeometryEncoded"].toBool()) return;

	for (auto g : graphs)
	{
		for (auto n : g->nodes)
		{
			auto mesh = g->getMesh(n->id);
			if (!mesh) continue;
			auto mesh_points = mesh->vertex_coordinates();
			QVector<ParameterCoord> encoding = n->property["encoding"].value< QVector<ParameterCoord> >();

			if (encoding.empty()) continue;

			// Generate consistent frames along curve
			Array1D_Vector4d coords;
			RMF rmf;
			if (n->type() == Structure::CURVE) rmf = Synthesizer::consistentFrame((Structure::Curve*)n, coords);
			int rmfCount = rmf.count();

			for (int i = 0; i < encoding.size(); i++){
				auto & sample = encoding[i];

				Vector3d startPoint;
				std::vector<Vector3d> frame;
				n->get(Eigen::Vector4d(sample.u, sample.v, 0, 0), startPoint, frame);
				Vector3f rayPos = Vector3f(startPoint[0], startPoint[1], startPoint[2]);

				Vector3d _X, _Y, _Z;
				if (n->type() == Structure::CURVE)
				{
					int idx = sample.u * (rmfCount - 1);
					_X = rmf.U[idx].r; _Y = rmf.U[idx].s; _Z = rmf.U[idx].t;
				}
				else
				{
					_X = frame[0]; _Z = frame[2]; _Y = cross(_Z, _X);
				}

				// double float
				Vector3f X(_X[0], _X[1], _X[2]), Y(_Y[0], _Y[1], _Y[2]), Z(_Z[0], _Z[1], _Z[2]);
				Vector3f rayDir;
				localSphericalToGlobal(X, Y, Z, sample.theta, sample.psi, rayDir);

				// Reconstructed point
				Vector3f isect = rayPos + (rayDir * sample.origOffset);
				mesh_points[SurfaceMesh::Vertex(i)] = isect.cast<double>();
			}

			mesh->update_face_normals();
			mesh->update_vertex_normals();
		}

		break; // only decode source
	}
}
