﻿#include "BatchProcess.h"
#include "myglobals.h"
#include <QCoreApplication>

Structure::ShapeGraph *shapeA, *shapeB;
static QVector<QColor> myrndcolors = rndColors2(100);

BatchProcess::BatchProcess(QString filename) : filename(filename)
{
	// Clean up my self
	connect(this, SIGNAL(finished()), this, SLOT(deleteLater()));

	// Rendering	
	renderer = new RenderingWidget(512, NULL);
	renderer->move(0, 0);
	renderer->show();

	// Progress
	pd = new QProgressDialog("Searching..", "Cancel", 0, 0);
	pd->setValue(0);
	pd->show();
	pd->connect(this, SIGNAL(jobFinished(int)), SLOT(setValue(int)));
	pd->connect(this, SIGNAL(allJobsFinished(int)), SLOT(deleteLater()));
}

QImage stitchImages(const QImage & a, const QImage & b, bool isVertical = false, int padding = 2, QColor background = Qt::white)
{
	int newWidth = isVertical ? (2 * padding) + std::max(a.width(), b.width()) : (3 * padding) + a.width() + b.width();
	int newHeight = isVertical ? (3 * padding) + a.height() + b.height() : (2 * padding) + std::max(a.height(), b.height());
	QImage img(newWidth, newHeight, QImage::Format_ARGB32_Premultiplied);
	QPainter painter;
	painter.begin(&img);
	painter.fillRect(img.rect(), background);
	if (isVertical)
	{
		painter.drawImage(padding, padding, a);
		painter.drawImage(padding, padding + a.height() + padding, b);
	}
	else
	{
		painter.drawImage(padding, padding, a);
		painter.drawImage(padding + a.width() + padding, padding, b);
	}
	painter.end();
	return img;
}

QImage drawText(QString message, QImage a, int x = 14, int y = 14, QColor color = Qt::black)
{
	QPainter painter;
	painter.setRenderHint(QPainter::Antialiasing);
	painter.setRenderHint(QPainter::HighQualityAntialiasing);
	painter.begin(&a);
	painter.setBrush(Qt::black);
	painter.setOpacity(0.2);
	painter.drawText(QPoint(x+1, y+1), message);
	painter.setOpacity(1.0);
	painter.setBrush(color);
	painter.drawText(QPoint(x, y), message);
	painter.end();
	return a;
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

	QElapsedTimer allTimer; allTimer.start();
	int allTime = 0;

	// Progress
	pd->setMaximum(jobsArray.size());

	for (int idx = 0; idx < jobsArray.size(); idx++ )
	{
		auto & j = jobsArray[idx];

		int searchTime = 0;

		/// Input Shapes:
		auto job = j.toObject(); if (job.isEmpty()) continue;
		auto source = job["source"].toString();
		auto target = job["target"].toString();
		auto title = job["title"].toString();

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
		QImage img;
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

			int thumbWidth = 256;

			auto cur_solution_img = stitchImages(
				renderer->render(shapeA).scaledToWidth(thumbWidth, Qt::TransformationMode::SmoothTransformation), 
				renderer->render(shapeB).scaledToWidth(thumbWidth, Qt::TransformationMode::SmoothTransformation));

			cur_solution_img = drawText(QString("cost = %1").arg(cost), cur_solution_img);

			img = stitchImages(img, cur_solution_img, true, 0);
		}

		QString msg = QString("Solution time (%1 s)").arg(double(searchTime) / 1000.0);
		int msgWidth = QFontMetrics(QFont()).width(msg) + 14;
		img = drawText(msg, img, img.width() - msgWidth, 14);

		auto output_file = QString("%1/%2.png").arg(outputPath).arg(title);
		QDir d(""); d.mkpath(QFileInfo(output_file).absolutePath());
		img.save(output_file);
	}

	allTime = allTimer.elapsed();

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
		json["resultsCount"] = 15;
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
	glClearColor(0.92f, 0.92f, 0.92f, 1);

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
	qglviewer::Vec cameraPos(-1.5,-1.75,1.0);

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