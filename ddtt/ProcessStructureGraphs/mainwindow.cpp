#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "globals.h"

#include <QDesktopWidget>
#include <QFileDialog>
#include <QTemporaryFile>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QSettings>
#include <QKeyEvent>
#include <QThread>
#include <QListWidget>
#include "Evaluator.h"

#include "OBB_Volume.h"

#include "StructureGraph.h"

#include "DBSCAN.hpp"
#include "optics.h"
#include "kmeans.h"

#include <string>

QString default_folder = "";
auto paired_colors = {"#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"};

int num_classes = 128;
QVector<QColor> class_color;

Q_DECLARE_METATYPE(OBB_Volume)

class MyViewer : public QGLViewer{
public:
	QMap <QString, QSharedPointer<Structure::Graph> > graphs;
	QStringList shape_names;
	QVector < QSharedPointer<RenderObject::Base> > debugs;

	QVector <QSharedPointer<SurfaceMeshModel> > models;

	QVector < QSharedPointer<QWidget> > widgets;

	void makeShapePairsWidget()
	{
		auto w = QSharedPointer<QWidget>(new QWidget);
		w->setMinimumSize(150,300);
		auto l1 = new QHBoxLayout;

		auto list1 = QSharedPointer<QListWidget>(new QListWidget());
		auto list2 = QSharedPointer<QListWidget>(new QListWidget());

		widgets << list1;
		widgets << list2;

		for (auto shape : shape_names)
		{
			list1->addItem(shape);
			list2->addItem(shape);
		}

		connect(list1.data(), &QListWidget::itemSelectionChanged, [&](){this->showPairs();});
		connect(list2.data(), &QListWidget::itemSelectionChanged, [&](){this->showPairs();});

		l1->addWidget(list1.data());
		l1->addWidget(list2.data());

		w->setLayout(l1);
		widgets << w;
		w->show();
	}

	void showPairs()
	{
		auto list1 = (QListWidget*)widgets.at(0).data();
		auto list2 = (QListWidget*)widgets.at(1).data();
		if (list1->selectedItems().empty() || list2->selectedItems().empty()) return;

		auto CalcHausdorffDist = [](const std::vector<Vector3> &vp, const std::vector<Vector3> &bbvp){
			double d1(0);
			for (unsigned int i = 0; i < vp.size(); i++) {
				double dd(std::numeric_limits<double>::max());
				for (unsigned int j = 0; j < bbvp.size(); j++) {
					dd = std::min(dd, (vp[i]-bbvp[j]).norm());
				}
				d1 = std::max(d1, dd);
			}
			double d2(0);
			for (unsigned int i = 0; i < bbvp.size(); i++) {
				double dd(std::numeric_limits<double>::max());
				for (unsigned int j = 0; j < vp.size(); j++) {
					dd = std::min(dd, (bbvp[i]-vp[j]).norm());
				}
				d2 = std::max(d2, dd);
			}
			return std::max(d1, d2);
		};

		for (auto g : graphs)
		{
			for (auto n : g->nodes)
				g->setColorFor(n->id, QColor(200, 200, 200));

			g->property["showMeshes"].setValue(true);
		}

		auto g1 = graphs[list1->selectedItems().front()->text()];

		for (auto g2 : graphs)
		{
			for (auto ni : g1->nodes){
				auto mesh = g1->getMesh(ni->id);

				auto obbi = OBB_Volume(mesh);

				ni->property["obb"].setValue(obbi);
				//for (auto p : obbi.corners<Vector3>()) g1->debug << starlab::PointSoup::drawPoint(p, 4);
			}

			for (auto nj : g2->nodes){
				auto meshj = g2->getMesh(nj->id);

				auto obbj = OBB_Volume(meshj);

				nj->property["obb"].setValue(obbj);
				//for (auto p : obbj.corners<Vector3>()) g2->debug << starlab::PointSoup::drawPoint(p, 4, Qt::green);
			}

			int ci = 0;

			QVector < QPair<QString,QString> > paring;
			QMap < QString, QSet<QString> > target_mapped;
		
			for (auto ni : g1->nodes)
			{
				auto obbi = ni->property["obb"].value<OBB_Volume>();
			
				double minDist = DBL_MAX;
				QString correspond = "";

				for (auto nj : g2->nodes)
				{
					auto obbj = nj->property["obb"].value<OBB_Volume>();

					double dist = CalcHausdorffDist(obbi.corners(), obbj.corners());

					if (dist < minDist)
					{
						minDist = dist;
						correspond = nj->id;
					}

					nj->vis_property["meshSolid"].setValue(true);
				}

				g1->property["showMeshes"].setValue(true);	g1->property["showNodes"].setValue(false);
				g2->property["showMeshes"].setValue(true);	g2->property["showNodes"].setValue(false);

				ni->vis_property["meshSolid"].setValue(true);

				g1->setColorFor(ni->id, class_color[ci]);
				g2->setColorFor(correspond, class_color[ci]);

				for (auto nid : target_mapped[correspond]) g1->setColorFor(nid, class_color[ci]);

				target_mapped[correspond].insert(ni->id);

				ci++;
			}

		}

		updateGL();
	}

	void keyPressEvent(QKeyEvent *event)
	{
		if (event->key() == Qt::Key_S)
		{
			int curWidth = std::min(this->width(), this->height());

			int thumbSize = 512;

			int numObjects = std::max(models.size(), graphs.size());
			int count = ceil(sqrt(double(numObjects)));

			int newWidth = thumbSize * count;

			this->setGeometry(this->pos().x(), this->pos().y(), newWidth, newWidth);
			qApp->processEvents();
			this->updateGL();
			this->updateGL();
			qApp->processEvents();
			this->grabFrameBuffer(false).save(this->windowTitle() + ".png");

			this->setGeometry(this->pos().x(), this->pos().y(), curWidth, curWidth);
		}

		QGLViewer::keyPressEvent(event);
	}

	MyViewer()
	{  
		QGLFormat glf = QGLFormat::defaultFormat();
		glf.setSamples(8);
		QGLFormat::setDefaultFormat(glf);

		int width = 800;
		int height = 800;

		QDesktopWidget wid;
		int screenWidth = wid.screen()->width();
		int screenHeight = wid.screen()->height();
		this->setGeometry((screenWidth / 2) - (width / 2), (screenHeight / 2) - (height / 2), width, height);
		this->show();

		//this->setGridIsDrawn();
		//this->setAxisIsDrawn();
	}

	void draw()
	{
		this->setBackgroundColor(Qt::white);

		//return;

		int numObjects = std::max(models.size(), graphs.size());
		if (numObjects == 0) return;

		int w = std::min(this->width(), this->height());
		int count = ceil(sqrt(double(numObjects)));
		int size = w / count;

		QVector < QRect > rects;
		for (int i = 0; i < count; i++){
			for (int j = 0; j < count; j++){
				rects << QRect(i*size, j*size, size, size);
			}
		}

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glEnable(GL_MULTISAMPLE);
		glEnable(GL_DEPTH_TEST);

		int ridx = 0;

		qglviewer::Vec startCameraPos(-1.75, -1.5, 0.8);
		double scaling = 0.08;

		for (auto m : models)
		{
			QRect r = rects.at(ridx++); int x = r.left(), y = r.top();
			glViewport(x, y, size, size);
			glScissor(x, y, size, size);

			// Setup camera
			auto s = m->bbox().diagonal().maxCoeff();
			auto cameraPos = startCameraPos;
			if (s > 1.0){
				auto delta = cameraPos;
				cameraPos += (delta * s) * scaling;
			}
			qglviewer::Camera cam;
			cam.setType(qglviewer::Camera::ORTHOGRAPHIC);
			cam.setScreenWidthAndHeight(size,size);
			cam.setSceneRadius(20.0f);
			cam.setSceneCenter(qglviewer::Vec(0, 0, 0.5));
			cam.lookAt(cam.sceneCenter());
			cam.setUpVector(qglviewer::Vec(0, 0, 1));
			cam.setPosition(cameraPos);
			cam.setViewDirection((cam.sceneCenter() - cameraPos).unit());
			cam.loadProjectionMatrix();
			cam.loadModelViewMatrix();

			// Draw segmented mesh
			glShadeModel(GL_FLAT);
			auto fcolor = m->face_property<QColor>("f:color");
			auto fnormal = m->face_property<Vector3>("f:normal");
			auto vpos = m->vertex_coordinates();
			glBegin(GL_TRIANGLES);
			for (auto f : m->faces())
			{
				auto c = fcolor[f];

				glColor3f(c.redF(), c.greenF(), c.blueF());
				glNormal3(fnormal[f]);
				for (auto v : m->vertices(f)) glVector3(vpos[v]);
			}
			glEnd();
		}

		for (auto gname : graphs.keys())
		{
			QRect r = rects.at(ridx++); int x = r.left(), y = r.top();
			glViewport(x, y, size, size);
			glScissor(x, y, size, size);

			auto cur_shape = graphs[gname];

			auto bbox = cur_shape->bbox();

			// Setup camera
			auto s = bbox.diagonal().maxCoeff(); 
			auto cameraPos = startCameraPos;
			if (s > 1.0){
				auto delta = cameraPos;
				cameraPos += (delta * s) * scaling;
			}
			qglviewer::Camera cam;
			cam.setType(qglviewer::Camera::ORTHOGRAPHIC);
			cam.setScreenWidthAndHeight(size,size);
			cam.setSceneRadius(20.0f);
			cam.setSceneCenter(qglviewer::Vec(0, 0, 0.5));
			cam.lookAt(cam.sceneCenter());
			cam.setUpVector(qglviewer::Vec(0, 0, 1));
			cam.setPosition(cameraPos);
			cam.setViewDirection((cam.sceneCenter() - cameraPos).unit());
			cam.loadProjectionMatrix();
			cam.loadModelViewMatrix();

			// Draw shape
			cur_shape->draw();

			for (auto d : cur_shape->debug) d->draw(*this);

			cur_shape->setColorAll(Qt::gray);
		}

		for (auto d : debugs) d->draw(*this);
	}
};

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
	ui->setupUi(this);

	// Default folder
	{
		const QString DEFAULT_DIR_KEY("default_dir");
		QString path = QDir(qApp->applicationDirPath()).absoluteFilePath("settings.ini");
		QSettings MySettings(path, QSettings::IniFormat);
		QString oldPath = MySettings.value(DEFAULT_DIR_KEY).toString();
		QString dirpath = QFileDialog::getExistingDirectory(this, tr("Open Directory"),
			oldPath, QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
		if (dirpath.isNull() || dirpath.isEmpty()) return;
		MySettings.setValue(DEFAULT_DIR_KEY, dirpath);
		MySettings.sync();
		default_folder = dirpath;
	}

	// Random colors
	{
		for (int i = 0; i < num_classes; i++){
			if (i >= paired_colors.size()) class_color.push_back(starlab::qRandomColor3());
			else{
				QColor c;
				c.setNamedColor(*(paired_colors.begin() + i));
				class_color.push_back(c);
			}
		}
	}

	connect(ui->fuzzyButton, &QPushButton::clicked, [&](){
		MyViewer * viewer = new MyViewer;

		viewer->setWindowTitle("fuzzy");

		auto folders = shapesInDataset(default_folder);

		for (auto folderName : folders.keys())
		{
			auto folder = folders[folderName];
			Structure::Graph g(folder["graphFile"].toString());

			viewer->shape_names << folderName;
			viewer->graphs[folderName] = QSharedPointer<Structure::Graph>(new Structure::Graph(folder["graphFile"].toString()));
		}

		viewer->makeShapePairsWidget();

	});

	connect(ui->baselineButton, &QPushButton::clicked, [&](){
		QMap <QString, QSharedPointer<Structure::Graph> > graphs;
		auto folders = shapesInDataset(default_folder);

		for (auto folderName : folders.keys())
		{
			auto folder = folders[folderName];
			graphs[folderName] = QSharedPointer<Structure::Graph>(new Structure::Graph(folder["graphFile"].toString()));
		}

		auto CalcHausdorffDist = [](const std::vector<Vector3> &vp, const std::vector<Vector3> &bbvp){
			double d1(0);
			for (unsigned int i = 0; i < vp.size(); i++) {
				double dd(std::numeric_limits<double>::max());
				for (unsigned int j = 0; j < bbvp.size(); j++) {
					dd = std::min(dd, (vp[i] - bbvp[j]).norm());
				}
				d1 = std::max(d1, dd);
			}
			double d2(0);
			for (unsigned int i = 0; i < bbvp.size(); i++) {
				double dd(std::numeric_limits<double>::max());
				for (unsigned int j = 0; j < vp.size(); j++) {
					dd = std::min(dd, (bbvp[i] - vp[j]).norm());
				}
				d2 = std::max(d2, dd);
			}
			return std::max(d1, d2);
		};

		std::vector<std::vector<std::pair<QString, QString>>> allMaps, allMapsLabel;
		for (auto g1 : graphs)
		{
			for (auto g2 : graphs)
			{
				if (g1->name() == g2->name()) continue;
// 				cout << "g1:" << g1->name().toStdString() << "  g2:" << g2->name().toStdString() << endl;

				std::vector<std::pair<QString, QString>> gmap, gmapl;

				for (auto ni : g1->nodes){
					auto mesh = g1->getMesh(ni->id);

					auto obbi = OBB_Volume(mesh);

					ni->property["obb"].setValue(obbi);
					//for (auto p : obbi.corners<Vector3>()) g1->debug << starlab::PointSoup::drawPoint(p, 4);
				}

				for (auto nj : g2->nodes){
					auto meshj = g2->getMesh(nj->id);

					auto obbj = OBB_Volume(meshj);

					nj->property["obb"].setValue(obbj);
					//for (auto p : obbj.corners<Vector3>()) g2->debug << starlab::PointSoup::drawPoint(p, 4, Qt::green);
				}


				QVector < QPair<QString, QString> > paring;
				QMap < QString, QSet<QString> > target_mapped;

				for (auto ni : g1->nodes)
				{
					auto obbi = ni->property["obb"].value<OBB_Volume>();

					double minDist = DBL_MAX;
					QString correspond = "";
					auto corrn = g2->nodes.front();

					for (auto nj : g2->nodes)
					{
						auto obbj = nj->property["obb"].value<OBB_Volume>();

						double dist = CalcHausdorffDist(obbi.corners(), obbj.corners());
						if (dist < minDist || correspond.isEmpty())
						{
							minDist = dist;
							correspond = nj->id;
							corrn = nj;
						}

						nj->vis_property["meshSolid"].setValue(true);
					}

						gmap.push_back(std::pair<QString, QString>(correspond, ni->id));
						gmapl.push_back(std::pair<QString, QString>(corrn->meta["label"].toString(), ni->meta["label"].toString()));

				}

				allMaps.push_back(gmap);
				allMapsLabel.push_back(gmapl);
			}
		}
		new Evaluator(default_folder, allMaps, allMapsLabel);
	});

	connect(ui->exportLabelsButton, &QPushButton::clicked, [&](){
		auto shapesDirName = default_folder;
		auto shapesDir = QDir(shapesDirName);

		auto folders = shapesInDataset(shapesDirName);
		auto outputFolder = QFileDialog::getExistingDirectory();

		QMap<QString, int> labelIndex, coarseLabelIndex;
		labelIndex["nullLabel"] = labelIndex.size();
		coarseLabelIndex["nullLabel"] = coarseLabelIndex.size();

		for (auto folderName : folders.keys())
		{
			auto & folder = folders[folderName];
			auto g = QSharedPointer<Structure::Graph>(new Structure::Graph(folder["graphFile"].toString()));

			int voffset = 0;

			SurfaceMeshModel sm;
			std::vector < std::pair<SurfaceMeshModel::Face, int> > faceLabels;
			std::vector < std::pair<SurfaceMeshModel::Face, int> > faceLabelsCoarse;

			for (auto n : g->nodes)
			{
				int label_idx = 0;
				int coarse_label_idx = 0;

				if (n->meta.contains("label"))
				{
					auto l = n->meta["label"].toString();
					if(!labelIndex.contains(l)) labelIndex[l] = labelIndex.size();
					label_idx = labelIndex[l];

					auto lc = l.split("-").front();
					if (!coarseLabelIndex.contains(lc)) coarseLabelIndex[lc] = coarseLabelIndex.size();
					coarse_label_idx = coarseLabelIndex[lc];
				}

				auto m = g->getMesh(n->id);

				for (auto v : m->vertices())
				{
					sm.add_vertex(m->vertex_coordinates()[v]);
				}

				for (auto f : m->faces())
				{
					std::vector<SurfaceMeshModel::Vertex> verts;

					for (auto v : m->vertices(f))
					{
						verts.push_back(SurfaceMeshModel::Vertex(v.idx() + voffset));
					}

					auto currFace = sm.add_face(verts);

					faceLabels.push_back(std::make_pair(currFace, label_idx));
					faceLabelsCoarse.push_back(std::make_pair(currFace, coarse_label_idx));
				}

				voffset += m->n_vertices();
			}

			auto f_label = sm.add_face_property<int>("f:label", 0);
			auto f_label_coarse = sm.add_face_property<int>("f:label_coarse", 0);

			for (auto fl : faceLabels) f_label[fl.first] = fl.second;
			for (auto fl : faceLabelsCoarse) f_label_coarse[fl.first] = fl.second;

			// Output meshes and labels
			QFile file(outputFolder + "/" + folderName + ".obj");
			QFile fileLabel(outputFolder + "/" + folderName + ".labels.txt");

			// Output mesh
			{
				if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return;
				QTextStream out(&file);

				out << QString("# Extracted from Structure Graph \n");

				for (auto v : sm.vertices())
				{
					auto p = sm.vertex_coordinates()[v];
					out << QString("v %1 %2 %3\n").arg(p.x()).arg(p.y()).arg(p.z());
				}

				for (auto f : sm.faces())
				{
					std::vector<SurfaceMeshModel::Vertex> verts;
					for (auto v : sm.vertices(f)) verts.push_back(v);

					out << QString("f %1 %2 %3\n").arg(verts[0].idx() + 1).arg(verts[1].idx() + 1).arg(verts[2].idx() + 1);
				}
			}

			// Output labels
			{
				if (!fileLabel.open(QIODevice::WriteOnly | QIODevice::Text)) return;
				QTextStream out(&fileLabel);

				for (auto f : sm.faces())
				{
					out << QString("%1 %2\n").arg(f_label_coarse[f]).arg(f_label[f]);
				}
			}
		}

		// output labels file
		{
			QFile file(outputFolder + "/labels.txt");
			if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return;
			QTextStream out(&file);

			for (auto key : labelIndex.keys())
			{
				QString l = key;
				QString cl = l.split("-").front().trimmed();
				out << QString("%1 %2 %3\n").arg(coarseLabelIndex[cl]).arg(labelIndex[l]).arg(l);
			}
		}
	});

    connect(ui->extractPly, &QPushButton::clicked,[&](){
		auto folders = shapesInDataset(default_folder);
		auto outputFolder = QFileDialog::getExistingDirectory();

        for(auto folderName : folders.keys())
        {
            auto folderInfo = folders[folderName];
            //debugBox(folderName);

            Structure::Graph g(folderInfo["graphFile"].toString());

			QTemporaryFile file("XXXXXX.obj");

			if (file.open()){
				QString objFilename = file.fileName();
				g.exportAsOBJ(objFilename);

				auto destination = QString("%1/%2").arg(outputFolder).arg(folderName /*+ ".ply"*/);

				// Convert 'obj' to 'ply' using meshlabserver
				//auto cmd = QString("\"C:/Program Files/VCG/MeshLab/meshlabserver.exe\" -i %1 -o %2").arg(objFilename).arg(destination);
				auto cmd = QString("\"C:/Program Files/VCG/MeshLab/meshconv.exe\" %1 -c ply -kz -kd -ascii -o %2").arg(objFilename).arg(destination);

				system(cmd.toLatin1());
			}
        }
    });

	connect(ui->clusterButton, &QPushButton::clicked, [&](){
		auto corrFilename = QFileDialog::getOpenFileName(this, "Open correspondence", "", "*.json");

		// Open JSON file
		QFile file;
		file.setFileName(corrFilename);
		if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;
		QJsonDocument jdoc = QJsonDocument::fromJson(file.readAll());
		auto corrArray = jdoc.toVariant().value<QVariantList>();

		QStringList edges;
		
		for (auto c : corrArray)
		{
			auto obj = c.toMap();
			auto i = obj["i"].toInt();
			auto j = obj["j"].toInt();
			auto cost = obj["cost"].toDouble();
			auto corr = obj["correspondence"].value<QVariantList>();

			for (auto match : corr)
			{
				auto matching = match.value<QVariantList>();
				auto nid1 = QString("%1:%2").arg(i).arg(matching.front().toString());
				auto nid2 = QString("%1:%2").arg(j).arg(matching.back().toString());

				auto edge = QString("%1 %2 %3").arg(nid1).arg(nid2).arg(cost);

				edges << edge;
			}
		}

		// Write graph file
		{
			QFile gfile(corrFilename + ".tgf");
			QTextStream out(&gfile);
			if (!gfile.open(QIODevice::WriteOnly | QIODevice::Text)) return;
			out << edges.join("\n");
		}
	});

	connect(ui->visualizeButton, &QPushButton::clicked, [&](){
		MyViewer * viewer = new MyViewer;

		viewer->setWindowTitle("ours");

		auto shapesDirName = default_folder;
		auto shapesDir = QDir(shapesDirName);

		auto folders = shapesInDataset(shapesDirName);
		auto shapesFilename = shapesDir.entryList(QStringList() << "*.txt", QDir::Files).join("");
		auto classFilename = shapesDir.entryList(QStringList() << "*.classsets", QDir::Files).join("");
		auto corrFilename = shapesDir.entryList(QStringList() << "*.json", QDir::Files).join("");

		// Get shapes
		{
			QFile file;
			file.setFileName(shapesDir.filePath(shapesFilename));
			if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;
			for (auto line : QString(file.readAll()).split("\n", QString::SkipEmptyParts))
				viewer->shape_names << line.split(" ").back();
		}

		// Load shapes
		{
			for (auto folderName : folders.keys())
			{
				auto & folder = folders[folderName];
				viewer->graphs[folderName] = QSharedPointer<Structure::Graph>(new Structure::Graph(folder["graphFile"].toString()));
			}
		}

		// Load part classification
		QMap < QString, QStringList > classes;
		QMap < QString, QString > partClassMap;
		{
			QFile file;
			file.setFileName(shapesDir.filePath(classFilename));
			if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;
			for (auto line : QString(file.readAll()).split("\n", QString::SkipEmptyParts))
			{
				auto lineTokens = line.split("\t", QString::SkipEmptyParts);
				if (lineTokens.at(0) != "//NODECLASS") continue;

				QString part = QString(lineTokens.at(1)).replace("\"", "");
				QString className = QString(lineTokens.at(2)).replace("\"", "");
				classes[className] << part;
				partClassMap[part] = className;
			}
		}

		// Assign colors
		{
			// Generate colors
			QMap<QString, QColor> classColor;
			for (int i = 0; i < classes.keys().size(); i++)
			{
				auto c = classes.keys().at(i);
				classColor[c] = class_color[i];
			}

			// Color the parts
			for (auto gname : viewer->graphs.keys())
			{
				int gidx = viewer->shape_names.indexOf(gname);
				auto g = viewer->graphs[gname];

				for (auto n : g->nodes)
				{
					QString partName = QString("%1:%2").arg(gidx).arg(n->id);
					
					g->setColorFor(n->id, classColor[partClassMap[partName]]);
					n->vis_property["meshSolid"].setValue(true);
				}

				g->property["showMeshes"].setValue(true);
				g->property["showNodes"].setValue(false);
			}

			viewer->update();
		}
	});

	connect(ui->greedyButton, &QPushButton::clicked, [&](){
		MyViewer * viewer = new MyViewer;

		viewer->setWindowTitle("greedy_basic");

		auto shapesDirName = default_folder;
		auto shapesDir = QDir(shapesDirName);

		auto folders = shapesInDataset(shapesDirName);

		// Load shapes
		{
			for (auto folderName : folders.keys())
			{
				auto & folder = folders[folderName];
				viewer->graphs[folderName] = QSharedPointer<Structure::Graph>(new Structure::Graph(folder["graphFile"].toString()));
			}
		}

		auto ps = QSharedPointer<starlab::PointSoup>(new starlab::PointSoup(12));

		// Cluster parts based on features (e.g. spatial, bounds, etc.)
		typedef std::vector<double> FeatureVector;
		std::vector<FeatureVector> features;
		QMap<Structure::Node *,int> part_ptrs;
		int pi = 0;

		for (auto gname : viewer->graphs.keys())
		{
			auto g = viewer->graphs[gname];

			auto gbox = g->robustBBox();

			auto graphColor = starlab::qRandomColor3();

			for (auto n : g->nodes)
			{
				FeatureVector fv;

				Vector3d center = (n->center() - gbox.min()).array() / gbox.sizes().array();
				//ps->addPoint(center, graphColor); // debug
				
				for (int i = 0; i < 3; i++) fv.push_back(center(i));

				features.push_back(fv);

				// Keep track of indices
				part_ptrs[n] = pi++;
			}
		}

		// DEBUG:
		//viewer->debugs << ps;

		if (features.empty()) return;

		// Clustering

		std::vector<int> labels(features.size(), 0);

		Eigen::MatrixXd M(features.size(), features.front().size());
		for (int i = 0; i < M.rows(); i++)
			for (int j = 0; j < M.cols(); j++)
				M(i, j) = features[i][j];

		// DBSCAN
		if (ui->checkbox1->isChecked())
		{
			double eps = ui->greedyParam->value();
			int min_elems = ui->greedyParam2->value();

			clustering::DBSCAN<Eigen::VectorXd, Eigen::MatrixXd> dbscan(eps, min_elems);

			dbscan.fit(M);

			auto labels = dbscan.get_labels();
			//num_classes = dbscan.num_classes();
		}

		// K-means
		else if (ui->checkbox2->isChecked())
		{
			int k = ui->greedyParam2->value();
			//num_classes = k;

			auto clusters = kmeansFast::kmeans(M, k);

			labels = clusters.cluster;
		}

		if (labels.empty()) return;

		for (auto gname : viewer->graphs.keys())
		{
			auto g = viewer->graphs[gname];
			for (auto n : g->nodes)
			{
				int class_id = labels[part_ptrs[n]];

				g->setColorFor(n->id, class_id < 0 ? Qt::black : class_color[class_id]);
				n->vis_property["meshSolid"].setValue(true);
			}

			g->property["showMeshes"].setValue(true);
			g->property["showNodes"].setValue(false);
		}

		// OPTICS
		{/*
			Optics optics(eps, min_elems);
			for (auto fv : features) optics.addSample(fv);
			optics.tick();
		*/}

	});
	
	connect(ui->conSegButton, &QPushButton::clicked, [&](){
		MyViewer * viewer = new MyViewer;

		viewer->setWindowTitle("consistent_seg");

		auto shapesDir = QDir(default_folder + "_consistent_seg");

		auto meshFiles = shapesDir.entryList(QStringList() << "*.ply", QDir::Files);

		QString p = shapesDir.absolutePath() + "/";

		for (int i = 0; i < meshFiles.size(); i++)
		{
			std::vector<int> face_class;

			// Read geometry and class per face
			auto mesh = QSharedPointer<SurfaceMeshModel>(new SurfaceMeshModel("mesh"));
			//mesh->read((p + meshFiles.at(i)).toStdString());
			{
				QFile plyFile(p + meshFiles.at(i));
				if (!plyFile.open(QIODevice::ReadOnly | QIODevice::Text)) continue;
				QTextStream in(&plyFile);
				auto lines = in.readAll().split("\n", QString::SkipEmptyParts);

				int nv = lines.at(2).split(" ").at(2).toInt();
				int nf = lines.at(6).split(" ").at(2).toInt();

				face_class.resize(nf, 0);

				int voffset = 11;

				for (int vi = 0; vi < nv; vi++)
				{
					auto vline = lines.at(vi + voffset).split(" ");
					mesh->add_vertex(Vector3(vline[0].toDouble(), vline[1].toDouble(), vline[2].toDouble()));
				}

				int foffset = voffset + nv;

				for (int fi = 0; fi < nf; fi++)
				{
					auto fline = lines.at(fi + foffset).split(" ");
					mesh->add_triangle(
						SurfaceMeshModel::Vertex( fline[1].toInt() ), 
						SurfaceMeshModel::Vertex( fline[2].toInt() ),
						SurfaceMeshModel::Vertex( fline[3].toInt() ));

					face_class[fi] = fline[4].toInt();
				}
			}

			auto fcolor = mesh->add_face_property<QColor>("f:color");
			for (auto f : mesh->faces())
			{
				fcolor[f] = class_color[face_class[f.idx()]];
			}

			// post-process
			mesh->updateBoundingBox();
			mesh->update_face_normals();
			mesh->update_vertex_normals();

			viewer->models << mesh;
		}
	});

	connect(ui->vkButton, &QPushButton::clicked, [&](){
		MyViewer * viewer = new MyViewer;

		viewer->setWindowTitle("vk_template");

		auto shapesDir = QDir(default_folder + "_vk_template");

		auto meshFiles = shapesDir.entryList(QStringList() << "*.off", QDir::Files);
		auto segFiles = shapesDir.entryList(QStringList() << "*.seg", QDir::Files);

		QString p = shapesDir.absolutePath() + "/";

		for (int i = 0; i < meshFiles.size(); i++)
		{
			// Read geometry
			auto mesh = QSharedPointer<SurfaceMeshModel>(new SurfaceMeshModel("mesh"));
			mesh->read((p + meshFiles.at(i)).toStdString());

			// Read segmentation
			auto filename = (p + segFiles.at(i)).toStdString();
			std::vector<double> m_aValues;
			{
				FILE * file = fopen(filename.c_str(), "rb");
				enum ShapeValueCacheType{
					SVTYPE_Unknown,
					SVTYPE_PerFace,			// single value for the whole face
					SVTYPE_PerIndexedFaces,	// for each face store indices 
					SVTYPE_PerVertex,		// defined at every vertex of a mesh
					SVTYPE_PerPointSet		// defined at every point in a set
				};
				auto ReadString = [](FILE * file){
					int l = 0;
					fread(&l, sizeof(int), 1, file);
					char * buffer = new char[l];
					fread(buffer, sizeof(char), l, file);
					std::string str(buffer);
					delete[] buffer;
					return str;
				};
				auto m_sValueName = ReadString(file);
				int m_iNDims, m_iNValues;
				ShapeValueCacheType m_CacheType;
				fread(&m_CacheType, sizeof(ShapeValueCacheType), 1, file);
				fread(&m_iNDims, sizeof(int), 1, file);
				fread(&m_iNValues, sizeof(int), 1, file);
				m_aValues.resize(m_iNDims*m_iNValues, 0);
				fread(&m_aValues[0], sizeof(double), m_iNDims*m_iNValues, file);
			}

			// Count number of classes
			//auto vec = m_aValues;
			//std::sort(vec.begin(), vec.end());
			//vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
			//QMap < int, int > imap;
			//for (auto i : vec) imap[i] = imap.size();
			//int num_classes = *std::max_element(m_aValues.begin(), m_aValues.end());

			auto fcolor = mesh->add_face_property<QColor>("f:color");
			for (auto f : mesh->faces())
			{
				fcolor[f] = class_color[m_aValues[f.idx()]];
			}

			// post-process
			mesh->updateBoundingBox();
			mesh->update_face_normals();
			mesh->update_vertex_normals();

			viewer->models << mesh;
		}
	});


	connect(ui->pairWiseButton, &QPushButton::clicked, [&](){
		new Evaluator(default_folder, false);
	});

	connect(ui->setWiseButton, &QPushButton::clicked, [&](){
		new Evaluator(default_folder, true);
	});
}

MainWindow::~MainWindow()
{
    delete ui;
}
