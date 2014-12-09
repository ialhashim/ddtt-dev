#include "BatchProcess.h"
#include "myglobals.h"
#include <QCoreApplication>
#include <QMatrix4x4>
#include "EncodeDecodeGeometry.h"

#include "AStarSearch.h"

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
    renderer->connect(this, SIGNAL(allJobsFinished()), SLOT(deleteLater()));

	// Progress
	pd = new QProgressDialog("Searching..", "Cancel", 0, 0);
	pd->setValue(0);
	pd->show();
	pd->connect(this, SIGNAL(jobFinished(int)), SLOT(setValue(int)));
    pd->connect(this, SIGNAL(allJobsFinished()), SLOT(deleteLater()));
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

		// Results:
		QMap <double, Energy::SearchNode*> sorted_solutions;

		/// Search solutions:
		Energy::GuidedDeformation egd;			
		
		// Load shapes
		shapeA = new Structure::ShapeGraph(source);
		shapeB = new Structure::ShapeGraph(target);

		// Set initial correspondence
		QVector<Energy::SearchNode> search_roots;
		Energy::SearchNode path(shapeA, shapeB, QStringList(), assignments);
		path.unassigned = path.unassignedList();
		search_roots << path;

		QElapsedTimer searchTimer; searchTimer.start();

		if (true)
		{
			sorted_solutions[0] = AStar::search(path);
		}
		else
		{
			// Search for all solutions
			egd.searchAll(shapeA, shapeB, search_roots);
		}

		emit(jobFinished(std::min(idx + 1, jobsArray.size() - 1)));
		QCoreApplication::processEvents();
		searchTime = searchTimer.elapsed();

		/// Rank solutions:
		if (sorted_solutions.empty())
		{
			auto all_solutions = egd.solutions();
			for (auto s : all_solutions)
			{
				double cost = roundDecimal(s->cost, 2);
				sorted_solutions[cost] = s;
			}
		}

		/// Draw top solutions:
		QImage img;
		for (int r = 0; r < resultsCount; r++)
		{
			if (r > sorted_solutions.size() - 1) continue; // less solutions than expected

			Energy::SearchNode * selected_path = NULL;
			double cost = 0;

			if (!egd.searchTrees.empty())
			{
				cost = sorted_solutions.keys().at(r);

				auto entire_path = egd.getEntirePath(sorted_solutions[cost]);
				egd.applySearchPath(entire_path);
				auto selected_path = entire_path.back();
			}
			else
			{
				selected_path = sorted_solutions.first();
			}

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

			// Show deformed
			{
				auto shapeAcopy = new Structure::ShapeGraph(*shapeA);
				for (auto n : shapeAcopy->nodes)
				{
					auto orig_mesh = n->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >().data();
					QSharedPointer<SurfaceMeshModel> new_mesh_ptr(orig_mesh->clone());
					new_mesh_ptr->updateBoundingBox();
					new_mesh_ptr->update_face_normals();
					new_mesh_ptr->update_vertex_normals();
					n->property["mesh"].setValue(new_mesh_ptr);
				}

				shapeAcopy->setAllControlPoints(shapeAcopy->animation.front());
				ShapeGeometry::encodeGeometry(shapeAcopy);
				shapeAcopy->setAllControlPoints(shapeAcopy->animation.back());
				ShapeGeometry::decodeGeometry(shapeAcopy);

				auto deformedImg = renderer->render(shapeAcopy).scaledToWidth(thumbWidth, Qt::TransformationMode::SmoothTransformation);
				deformedImg = drawText("[Deformed source]", deformedImg, 14, deformedImg.height() - 20);

				cur_solution_img = stitchImages(cur_solution_img, deformedImg);
			}

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

    auto bbox = cur_shape->bbox();

	// Setup camera
    qglviewer::Vec cameraPos(-1.5,-1.75,1.0);

	auto s = bbox.diagonal().maxCoeff();
	if (s > 1.0)
	{
		auto delta = cameraPos;
		cameraPos += (delta * s) * 0.2;
	}

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
