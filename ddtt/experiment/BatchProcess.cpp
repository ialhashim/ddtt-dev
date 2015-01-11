#include "BatchProcess.h"
#include "myglobals.h"
#include <QGuiApplication>
#include <QApplication>
#include <QGridLayout>
#include <QPushButton>
#include <QDialogButtonBox>
#include <QListWidget>
#include <QMatrix4x4>
#include "EncodeDecodeGeometry.h"

#include "AStarSearch.h"

Structure::ShapeGraph *shapeA, *shapeB;
static QVector<QColor> myrndcolors = rndColors2(100);

QDialog * dialog = nullptr;
BatchProcess * bp = nullptr;

void BatchProcess::init()
{
	bp = this;

	// Default options
	resultsCount = 6;
	outputPath = "outputPath";
	isSwapped = false;
	isSaveReport = true;
	thumbWidth = 256;

	// Clean up my self
	connect(this, SIGNAL(finished()), this, SLOT(deleteLater()));

	// Rendering	
    renderer = new RenderingWidget(1024, NULL);
	renderer->move(0, 0);
	renderer->connect(this, SIGNAL(allJobsFinished()), SLOT(deleteLater()));

	// Progress
	pd = new QProgressDialog("Searching..", "Cancel", 0, 0);
	pd->setValue(0);
	pd->connect(this, SIGNAL(jobFinished(int)), SLOT(setValue(int)));
	pd->connect(this, SIGNAL(allJobsFinished()), SLOT(deleteLater()));
	pd->connect(this, SIGNAL(setLabelText(QString)), SLOT(setLabelText(QString)));

	// Show progress
	pd->show();
	renderer->show();
}

BatchProcess::BatchProcess(QString filename) : jobfilename(filename)
{
	init();

	// Load job's data
	{
		QFile file;
		file.setFileName(filename);
		if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;
		QJsonDocument jdoc = QJsonDocument::fromJson(file.readAll());
		auto json = jdoc.object();

		resultsCount = json["resultsCount"].toInt();
		outputPath = json["outputPath"].toString();
		isSwapped = json["isSwap"].toBool();
		isSaveReport = json["isSaveReport"].toBool();
		thumbWidth = std::max(256, json["thumbWidth"].toInt());
		jobsArray = json["jobs"].toArray();

		// Clear past results in output folder
		QDir d(""); d.mkpath(outputPath);
		QDir dir(outputPath);
		for (auto filename : dir.entryList(QDir::Files))
			if (filename.endsWith(".png") || filename.endsWith(".txt"))
				dir.remove(dir.absolutePath() + "/" + filename);
	}

	// Selective jobs dialog
	if (QGuiApplication::queryKeyboardModifiers().testFlag(Qt::ShiftModifier))
	{
		dialog = new QDialog;
		dialog->setWindowTitle("Select jobs");
		dialog->setModal(true);
		QDialogButtonBox * box = new QDialogButtonBox(Qt::Horizontal);
		auto mainLayout = new QGridLayout;
		auto button = new QPushButton("OK");
		auto buttonSave = new QPushButton("Save selected..");
		button->setDefault(true);
		box->addButton(buttonSave, QDialogButtonBox::NoRole);
		box->addButton(button, QDialogButtonBox::AcceptRole);
		mainLayout->addWidget(box, 1, 0);
		dialog->setLayout(mainLayout);
		dialog->connect(box, SIGNAL(accepted()), SLOT(accept()));
		dialog->setAttribute(Qt::WA_DeleteOnClose);

		// List jobs
		QListWidget * list = new QListWidget;
		for (auto & j : jobsArray){
			auto job = j.toObject(); if (job.isEmpty()) continue;
			QListWidgetItem* item = new QListWidgetItem(job["title"].toString(), list);
			item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
			item->setCheckState(Qt::Unchecked);
			item->setData(Qt::UserRole, job);
		}
		mainLayout->addWidget(list, 0, 0);

		// Actions
		dialog->connect(box, &QDialogButtonBox::accepted, [&]{
			QJsonArray selectedJobs;
			for (int row = 0; row < list->count(); row++){
				QListWidgetItem *item = list->item(row);
				if (item->checkState() == Qt::Checked) selectedJobs << item->data(Qt::UserRole).toJsonObject();
			}
			if (!selectedJobs.isEmpty()) bp->setJobsArray(selectedJobs);
		});
		dialog->connect(box, SIGNAL(accepted()), SLOT(accept()));
		dialog->connect(buttonSave, &QPushButton::clicked, [&]{
			QJsonArray selectedJobs;
			for (int row = 0; row < list->count(); row++){
				QListWidgetItem *item = list->item(row);
				if (item->checkState() == Qt::Checked) selectedJobs << item->data(Qt::UserRole).toJsonObject();
			}
			if (!selectedJobs.isEmpty()) bp->setJobsArray(selectedJobs);
			bp->exportJobFile("currentJobs.json");
			QMessageBox::information(0, "Save", "Jobs saved.");
		});

		// Show
		dialog->exec();
	}	
}

void BatchProcess::run()
{
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

		if (isSwapped)
		{
			std::swap(source, target);
			auto splitTitle = title.split("-");
			title = splitTitle.size() > 1 ? splitTitle.back() + "-" + splitTitle.front() : title;
		}

		emit(setLabelText(QString("Corresponding: %1...").arg(title)));

		/// Initial Assignments:
		Energy::Assignments assignments;
		for (auto a : job["assignments"].toArray())
		{
			QStringList la, lb;
			for (auto part : a.toObject()["source"].toArray().toVariantList()) la << part.toString();
			for (auto part : a.toObject()["target"].toArray().toVariantList()) lb << part.toString();
			if (isSwapped) std::swap(la, lb);
			assignments.push_back(qMakePair(la, lb));
		}

		// Job report
		QVariantMap jobReport;

		// Results:
		QMap <double, Energy::SearchNode> sorted_solutions;
		QVector < QVector <Energy::SearchNode> > solution_vec;

		/// Search solutions:
		Energy::GuidedDeformation egd;			
		
		// Load shapes
		auto shapeA = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(source));
		auto shapeB = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(target));

		if (shapeA->nodes.isEmpty() || shapeB->nodes.isEmpty()) continue;

        if(job["isAnisotropy"].toBool())
        {
            QMatrix4x4 mat;
            auto bboxA = shapeA->bbox();
            auto bboxB = shapeB->bbox();
            Vector3 s = bboxA.diagonal().array() / bboxB.diagonal().array();
            mat.scale(s.x(), s.y(), s.z());
            shapeB->transform(mat, true);
        }

		// Set initial correspondence
		QVector<Energy::SearchNode> search_roots;
		Energy::SearchNode path(shapeA, shapeB, QSet<QString>(), assignments);
		path.unassigned = path.unassignedList();
		search_roots << path;

		QElapsedTimer searchTimer; searchTimer.start();

		bool isSearchAstar = true;
		unsigned int numNodesSearched = 0;

		if ( isSearchAstar )
		{
			for (auto & solution : AStar::search(path, 100, &numNodesSearched))
			{
				egd.origShapeA = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(*shapeA));
				egd.origShapeB = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(*shapeB));

				solution_vec.push_back(QVector<Energy::SearchNode>());
				for (auto & state : solution) solution_vec.back() << state;
				QVector<Energy::SearchNode*> ptrs;
				for (auto & node : solution_vec.back()) ptrs << &node;
				egd.applySearchPath(ptrs);
				sorted_solutions[solution.back().energy] = solution_vec.back().back();
			}
		}
		else
		{
			// Search for all solutions
			egd.searchAll(shapeA.data(), shapeB.data(), search_roots);

			auto all_solutions = egd.solutions();
			numNodesSearched = all_solutions.size();

			for (auto s : all_solutions)
			{
				double cost = roundDecimal(s->energy, 2);
				sorted_solutions[cost] = *s;
			}
		}

		emit(jobFinished(std::min(idx + 1, jobsArray.size() - 1)));
		QCoreApplication::processEvents();
		searchTime = searchTimer.elapsed();

		/// Draw top solutions:
		QImage img;
		for (int r = 0; r < resultsCount; r++)
		{
			if (r > sorted_solutions.size() - 1) continue; // less solutions than expected

			Energy::SearchNode * selected_path = NULL;
			double cost = sorted_solutions.keys().at(r);

			if (!isSearchAstar)
			{
				auto entire_path = egd.getEntirePath(&sorted_solutions[cost]);
				egd.applySearchPath(entire_path);
				selected_path = entire_path.back();
			}
			else
			{
				selected_path = &sorted_solutions[cost];
			}

			auto shapeA = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(*selected_path->shapeA.data()));
			auto shapeB = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(*selected_path->shapeB.data()));

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

			auto cur_solution_img = stitchImages(
				renderer->render(shapeA.data()).scaledToWidth(thumbWidth, Qt::TransformationMode::SmoothTransformation), 
				renderer->render(shapeB.data()).scaledToWidth(thumbWidth, Qt::TransformationMode::SmoothTransformation));

			// Show deformed
			{
				auto shapeAcopy = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(*shapeA));
				for (auto n : shapeAcopy->nodes)
				{
					if (!n->property.contains("mesh")) continue;

					auto orig_mesh = n->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >().data();
					if (!orig_mesh) continue;

					QSharedPointer<SurfaceMeshModel> new_mesh_ptr(orig_mesh->clone());
					new_mesh_ptr->updateBoundingBox();
					new_mesh_ptr->update_face_normals();
					new_mesh_ptr->update_vertex_normals();
					n->property["mesh"].setValue(new_mesh_ptr);
				}

				// Reset shape
				shapeAcopy->setAllControlPoints(shapeAcopy->animation.front());

				// Adjust for splitting cases
				for (auto n : shapeAcopy->nodes){
					if (n->id.contains("@")){
						QString origNode = n->id.split("@").front();
						n->setControlPoints(shapeAcopy->getNode(origNode)->controlPoints());
					}
				}

				// Compute geometry encoding
				ShapeGeometry::encodeGeometry(shapeAcopy.data());

				// Deform abstractions
				shapeAcopy->setAllControlPoints(shapeAcopy->animation.back());

				// Compute deformed underlying surface
				ShapeGeometry::decodeGeometry(shapeAcopy.data());

				auto deformedImg = renderer->render(shapeAcopy.data()).scaledToWidth(thumbWidth, Qt::TransformationMode::SmoothTransformation);
				deformedImg = drawText("[Deformed source]", deformedImg, 14, deformedImg.height() - 20);

				cur_solution_img = stitchImages(cur_solution_img, deformedImg);
			}

			// Details of solution if requested
			if ( isSaveReport )
			{
				EvaluateCorrespondence::evaluate(selected_path);
				QVariantMap details = selected_path->shapeA->property["costs"].value<QVariantMap>();
				QStringList reportItems;
				for (auto key : details.keys())
				{
					auto detail = details[key].toString();

					// Visualization
					if (key == "zzShapeCost")
						cur_solution_img = drawText(QString("[%1]").arg(detail), cur_solution_img, 12, cur_solution_img.height() - 20);

					if (key == "zzShapeCost") key += QString("(%1)").arg(r);
					reportItems += key + " : " + detail;
				}
				
				// Record found correspondence
				QString mapping;
				for (auto key : selected_path->mapping.keys())
					mapping += QString("(%1,%2) ").arg(key).arg(selected_path->mapping[key]);
				reportItems += mapping;
				jobReport[QString("solution-%1").arg(r)] = reportItems;
			}

			// Visualization
			cur_solution_img = drawText(QString("s %2: cost = %1").arg(cost).arg(r), cur_solution_img);
			img = stitchImages(img, cur_solution_img, true, 0);
		}

		QString msg = QString("Solution time (%1 s)").arg(double(searchTime) / 1000.0);
		int msgWidth = QFontMetrics(QFont()).width(msg) + 14;
		img = drawText(msg, img, img.width() - msgWidth, 14);
		img = drawText(QString("steps %1").arg(numNodesSearched), img, img.width() - msgWidth, 30);

		auto output_file = QString("%1/%2.png").arg(outputPath).arg(title);
		img.save(output_file);

		std::cout << " Saving image: " << qPrintable(output_file);

		if ( isSaveReport )
		{
			auto report_file = QString("%1/%2.txt").arg(outputPath).arg(title);
			QFile file(report_file);
			file.open(QFile::WriteOnly | QFile::Text);
			QTextStream out(&file);
			for (auto key : jobReport.keys())
			{
				out << key << ":" << "\n";
				for (auto item : jobReport[key].toStringList())
				{
					out << item << "\n";
				}
				out << "\n======================\n";
			}
		}
	}

	allTime = allTimer.elapsed();

	emit(jobFinished(jobsArray.size()));
    emit(allJobsFinished());
	emit(reportMessage(QString("Batch process time (%1 s)").arg(double(allTimer.elapsed()) / 1000), 0));
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
		json["resultsCount"] = 6;
		json["isSaveReport"] = true;
	}

	auto jobs = json["jobs"].toArray();

	auto jj = QJsonObject::fromVariantMap(job);
	jobs.push_front(jj);

	json["jobs"] = jobs;

	QJsonDocument saveDoc(json);
	QFile saveFile(filename);
	if (!saveFile.open(QIODevice::WriteOnly)) return;
	saveFile.write(saveDoc.toJson());
}

void BatchProcess::exportJobFile(QString filename)
{
	if (filename.isEmpty()) filename = jobfilename;

	QJsonObject json;

	json["outputPath"] = outputPath;
	json["resultsCount"] = resultsCount;
	json["isSaveReport"] = isSaveReport;
	json["jobs"] = jobsArray;

	QJsonDocument saveDoc(json);
	QFile saveFile(filename);
	if (!saveFile.open(QIODevice::WriteOnly)) return;
	saveFile.write(saveDoc.toJson());
}

void BatchProcess::setJobsArray(QJsonArray fromJobsArray)
{
	jobsArray = fromJobsArray;
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
		cameraPos += (delta * s) * 0.3;
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

// All automatic option
BatchProcess::BatchProcess(QString sourceFilename, QString targetFilename)
{
	init();

	QDir d(""); d.mkpath(outputPath);

	QVariantMap job;
	job["source"].setValue(sourceFilename);
	job["target"].setValue(targetFilename);
	job["title"].setValue(QFileInfo(sourceFilename).baseName() + "-" + QFileInfo(targetFilename).baseName());

	jobsArray.push_back(QJsonObject::fromVariantMap(job));
}
