#pragma once

// Simplify runtime debugging
#include <QMessageBox>
template<typename DataType>
static inline void debugBox( DataType message ){
    QMessageBox msgBox;
    msgBox.setTextFormat(Qt::RichText);
    msgBox.setText( QString("%1").arg(message) );
    msgBox.setStandardButtons(QMessageBox::Ok);
    msgBox.exec();
}
static inline void debugBoxList( QStringList messages ){
    debugBox( messages.join("<br>") );
}
template<typename Container>
static inline void debugBoxVec( Container data ){
    QStringList l;
    for(auto d : data) l << QString("%1").arg(d);
    debugBoxList(l);
}
template<typename Container2D>
static inline void debugBoxVec2( Container2D data, int limit = -1 ){
    QStringList l;
    for(auto row : data){
        QStringList line;
        for(auto d : row) line << QString("%1").arg( d );
        l << QString("%1").arg( line.join(", ") );
        if(limit > 0 && l.size() == limit) break;
    }
    if(limit > 0 && data.size() - l.size() > 0) l << QString("... (%1) more").arg(data.size() - l.size());
    debugBoxList(l);
}

// Loading XML files
#include <QDir>
typedef QMap<QString, QVariantMap> DatasetMap;
static inline DatasetMap shapesInDataset(QString datasetPath)
{
    DatasetMap dataset;

    QDir datasetDir(datasetPath);
    QStringList subdirs = datasetDir.entryList(QDir::Dirs | QDir::NoSymLinks | QDir::NoDotAndDotDot);

    for(auto subdir : subdirs)
    {
        // Special folders
        if (subdir == "corr") continue;

        QDir d(datasetPath + "/" + subdir);

        // Check if no graph is in this folder
        if (d.entryList(QStringList() << "*.xml", QDir::Files).isEmpty()) continue;

        auto xml_files = d.entryList(QStringList() << "*.xml", QDir::Files);

        dataset[subdir]["Name"] = subdir;
        dataset[subdir]["graphFile"] = d.absolutePath() + "/" + (xml_files.empty() ? "" : xml_files.front());
        dataset[subdir]["thumbFile"] = d.absolutePath() + "/" + d.entryList(QStringList() << "*.png", QDir::Files).join("");
        dataset[subdir]["objFile"] = d.absolutePath() + "/" + d.entryList(QStringList() << "*.obj", QDir::Files).join("");
    }

    return dataset;
}


// My customized viewer
#include <qglviewer/qglviewer.h>
#include <QDesktopWidget>
#include <QListWidget>
#include <QHBoxLayout>
#include <QKeyEvent>
#include "StructureGraph.h"

#include "OBB_Volume.h"
Q_DECLARE_METATYPE(OBB_Volume)

class MyViewer : public QGLViewer{
public:
	QStringList shape_names;
	QMap <QString, QSharedPointer<Structure::Graph> > graphs;
	QVector < QSharedPointer<RenderObject::Base> > debugs;
	QVector <QSharedPointer<SurfaceMeshModel> > models;
	QVector < QSharedPointer<QWidget> > widgets;

	static int num_classes;
	static QVector<QColor> class_color;

	void makeShapePairsWidget()
	{
		auto w = QSharedPointer<QWidget>(new QWidget);
		w->setMinimumSize(150, 300);
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

		connect(list1.data(), &QListWidget::itemSelectionChanged, [&](){this->showPairs(); });
		connect(list2.data(), &QListWidget::itemSelectionChanged, [&](){this->showPairs(); });

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

		for (auto g : graphs)
		{
			for (auto n : g->nodes)
				g->setColorFor(n->id, QColor(200, 200, 200));

			g->property["showMeshes"].setValue(true);
		}

		auto g1 = graphs[list1->selectedItems().front()->text()];
		//auto g2 = graphs[list2->selectedItems().front()->text()];

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

			QVector < QPair<QString, QString> > paring;
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
			cam.setScreenWidthAndHeight(size, size);
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
			cam.setScreenWidthAndHeight(size, size);
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
