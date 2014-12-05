#include "BatchProcess.h"
#include "myglobals.h"
#include <QCoreApplication>

Structure::ShapeGraph *shapeA, *shapeB;
static QVector<QColor> myrndcolors = rndColors2(100);

BatchProcess::BatchProcess(QString filename) : filename(filename)
{
	// Clean up my self
	connect(this, SIGNAL(finished()), this, SLOT(deleteLater()));

	// Progress
	pd = new QProgressDialog("Searching..", "Cancel", 0, 0);
	pd->setValue(0);
	pd->show();
	pd->connect(this, SIGNAL(jobFinished(int)), SLOT(setValue(int)));
	pd->connect(this, SIGNAL(allJobsFinished(int)), SLOT(deleteLater()));

	// Rendering	
	renderer = new RenderingWidget(256, NULL);
	renderer->moveToThread(this);
	renderer->show();
}

void BatchProcess::run()
{
	QFile file;
	file.setFileName(filename);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;
	QJsonDocument jdoc = QJsonDocument::fromJson(file.readAll());
	auto json = jdoc.object();

	int resultsCount = json["resultsCount"].toInt();
	QString outputPath = json["outputPath"].toString();
	auto jobsArray = json["jobs"].toArray();

	// Progress
	pd->setMaximum(jobsArray.size());

	for (int idx = 0; idx < jobsArray.size(); idx++ )
	{
		auto & j = jobsArray[idx];

		int allTime = 0, searchTime = 0;
		QElapsedTimer allTimer; allTimer.start();

		/// Input Shapes:
		auto job = j.toObject(); if (job.isEmpty()) continue;
		auto source = job["source"].toString();
		auto target = job["target"].toString();

		/// Initial Assignments:
		Energy::Assignments assignments;
		for (auto a : job["assignments"].toArray())
		{
			QStringList la, lb;
			for (auto part : a.toObject()["source"].toArray().toVariantList()) la << part.toString();
			for (auto part : a.toObject()["target"].toArray().toVariantList()) lb << part.toString();
			assignments.push_back(qMakePair(la, lb));
		}

		/// Search solutions:
		Energy::GuidedDeformation egd;			
		
		// Load shapes
		shapeA = new Structure::ShapeGraph(source);
		shapeB = new Structure::ShapeGraph(target);

		// Set initial correspondence
		Energy::SearchPath path(shapeA, shapeB, QStringList(), assignments);
		egd.search_paths << path;

		// Search for all solutions
		QElapsedTimer searchTimer; searchTimer.start();
		egd.searchAll();
		emit(jobFinished(idx));
		searchTime = searchTimer.elapsed();
		
		/// Rank solutions:
		QMap <double, Energy::SearchPath*> sorted_solutions;
		{
			auto all_solutions = egd.solutions();
			for (auto s : all_solutions) sorted_solutions[s->cost] = s;
		}

		/// Draw top solutions:
		QVector<QImage> images;
		for (int r = 0; r < resultsCount; r++)
		{
			auto cost = sorted_solutions.keys().at(r);

			auto entire_path = Energy::SearchPath::getEntirePath(sorted_solutions[cost], egd.search_paths);
			egd.applySearchPath(entire_path);
			auto selected_path = entire_path.back();

			shapeA = selected_path->shapeA;
			shapeB = selected_path->shapeB;

			// Color corresponded nodes
			{
				// Grey out
				for (auto n : shapeA->nodes){
					shapeA->setColorFor(n->id, QColor(255, 255, 255, 10));
					n->vis_property["meshSolid"].setValue(false);
					if (n->type() == Structure::SHEET) ((Structure::Sheet*)n)->surface.quads.clear();
				}
				for (auto n : shapeB->nodes){
					shapeB->setColorFor(n->id, QColor(255, 255, 255, 10));
					n->vis_property["meshSolid"].setValue(false);
					if (n->type() == Structure::SHEET) ((Structure::Sheet*)n)->surface.quads.clear();
				}

				shapeA->property["showNodes"].setValue(false);
				shapeB->property["showNodes"].setValue(false);

				// Assign colors based on target
				int ci = 0;
				for (auto & relation : shapeB->relations)
				{
					QColor color = myrndcolors[ci++];
					for (auto nid : relation.parts)
					{
						shapeB->setColorFor(nid, color);
						shapeB->getNode(nid)->vis_property["meshSolid"].setValue(true);
					}
				}

				// Color matching source
				for (auto spart : selected_path->mapping.keys())
				{
					auto tpart = selected_path->mapping[spart];
					if (tpart == Structure::null_part) continue;

					auto color = shapeB->getNode(tpart)->vis_property["meshColor"].value<QColor>();

					shapeA->setColorFor(spart, color);
					shapeA->getNode(spart)->vis_property["meshSolid"].setValue(true);
				}
			}

			auto simg = renderer->render(shapeA);
			auto timg = renderer->render(shapeB);

			QImage img(simg.width() + timg.width(), simg.height(), QImage::Format_ARGB32_Premultiplied);
			QPainter painter;
			painter.begin(&img);
			painter.drawImage(0, 0, simg);
			painter.drawImage(simg.width(), 0, timg);
			painter.setBrush(Qt::black);
			painter.drawText(10, 10, QString::number(cost));
			painter.end();

			QString number; number.sprintf("%04d", r);
			img.save(QString("result_%1.png").arg(number));
			images << img;
		}

		allTime = allTimer.elapsed();
	}

	emit(jobFinished(jobsArray.size()));
	emit(allJobsFinished());

	renderer->hide();
	renderer->deleteLater();
}

void BatchProcess::appendJob(QVariantMap job, QString filename)
{
	QFile file;
	file.setFileName(filename);
	if (!file.open(QIODevice::ReadWrite | QIODevice::Text)) return;
	QJsonDocument jdoc = QJsonDocument::fromJson(file.readAll());
	file.close();

	auto json = jdoc.object();

	// Default values if needed
	if (!json.contains("outputPath"))
	{
		json["outputPath"] = QString("outputPath");
		json["resultsCount"] = 10;
	}

	auto jobs = json["jobs"].toArray();

	auto jj = QJsonObject::fromVariantMap(job);
	auto tt = jj.keys();
	jobs.push_back(jj);

	json["jobs"] = jobs;

	QJsonDocument saveDoc(json);
	QFile saveFile(filename);
	if (!saveFile.open(QIODevice::WriteOnly)) return;
	saveFile.write(saveDoc.toJson());
}

QImage RenderingWidget::render(Structure::ShapeGraph * shape)
{
	this->cur_shape = shape;
	buffer = QImage();
	this->update();
	while (buffer.isNull());
	return buffer;
}

void RenderingWidget::initializeGL()
{
	initializeOpenGLFunctions();
	glClearColor(1, 1, 1, 1);

	glDisable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);

	// Setup lights and material
	GLfloat ambientLightColor[] = { 0.2f, 0.2f, 0.2f, 1 };
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLightColor);

	GLfloat diffuseLightColor[] = { 0.9f, 0.9f, 0.9f, 1 };
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLightColor);

	GLfloat specularLightColor[] = { 0.95f, 0.95f, 0.95f, 1 };
	glLightfv(GL_LIGHT0, GL_SPECULAR, specularLightColor);

	float posLight0[] = { 3, 3, 3, 0 };
	glLightfv(GL_LIGHT0, GL_POSITION, posLight0);

	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);

	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

	// Specular lighting
	float specReflection[] = { 0.8f, 0.8f, 0.8f, 1.0f };
	glMaterialfv(GL_FRONT, GL_SPECULAR, specReflection);
	glMateriali(GL_FRONT, GL_SHININESS, 56);
}

void RenderingWidget::paintGL()
{
	if (!cur_shape) return;
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_MULTISAMPLE);

	// Setup camera
	qglviewer::Vec cameraPos(-2,-2,1.5);

	qglviewer::Camera cam;
	cam.setType(qglviewer::Camera::ORTHOGRAPHIC);
	cam.setScreenWidthAndHeight(this->width(), this->height());
	cam.setSceneRadius(20.0f);	
	cam.setSceneCenter(qglviewer::Vec(0, 0, 0.5));
	cam.lookAt(cam.sceneCenter());
	cam.setUpVector(qglviewer::Vec(0, 0, 1));
	cam.setPosition(cameraPos);
	cam.setViewDirection((cam.sceneCenter()-cameraPos).unit());
	cam.loadProjectionMatrix();
	cam.loadModelViewMatrix();

	// Draw shape
	cur_shape->draw();

	// Return QImage of buffer
	auto size = this->size();
	bool alpha_format = true, include_alpha = true;
	QImage img(size, (alpha_format && include_alpha) ? QImage::Format_ARGB32_Premultiplied : QImage::Format_RGB32);
	int w = size.width();
	int h = size.height();
	glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, img.bits());
	convertFromGLImage(img, w, h, alpha_format, include_alpha);

	buffer = img;

	// Reset
	cur_shape = NULL;
}

RenderingWidget::RenderingWidget(int width, QWidget * parent) : cur_shape(NULL), QOpenGLWidget(parent)
{
	setFixedSize(QSize(width, width));
}
