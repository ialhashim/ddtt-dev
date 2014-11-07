#include "deform.h"
#include "deform-widget.h"
#include "ui_deform-widget.h"
#include "RenderObjectExt.h"

#define GEO_HEAT_USE_EIGEN_SOLVER
#include "../surfacemesh_filter_geoheat/GeoHeatHelper.h"

#include "Constraint.h"
#include "Force.h"

#ifndef DISABLE_TETGEN
#include "tetgenLib.h"
#endif

QTimer *timer = NULL;

void deform::create()
{
	if (widget) return;

	// Viewer
	worldRadius = mesh()->bbox().diagonal().norm();
	{
		//drawArea()->setAxisIsDrawn(true);
		drawArea()->camera()->setType(qglviewer::Camera::PERSPECTIVE);

		drawArea()->camera()->setUpVector(qglviewer::Vec(0, 0, 1));
		drawArea()->camera()->setPosition(qglviewer::Vec(worldRadius * 2, -worldRadius * 2, worldRadius * 1.5));
		drawArea()->camera()->lookAt(qglviewer::Vec());
		drawArea()->camera()->setSceneRadius(worldRadius);
		drawArea()->camera()->showEntireScene();
	}

	// Rendering
	{
		drawArea()->setRenderer(mesh(), "Flat Wire");
	}

	ModePluginDockWidget * dockwidget = new ModePluginDockWidget("Deform", mainWindow());
	DeformWidget * dw = new DeformWidget();
	widget = dw;

	dockwidget->setWidget(widget);
	mainWindow()->addDockWidget(Qt::RightDockWidgetArea, dockwidget);

	// UI
	connect(dw->ui->createHandle, &QPushButton::released, [&]{
		auto vid = Vertex(((DeformWidget *)widget)->ui->vertexID->value());
		this->create_handle(mesh()->vertex_coordinates()[vid], vid.idx());
	});
	connect(dw->ui->deleteHandle, &QPushButton::released, [&]{
		auto newHandles = handles;
		newHandles.clear();

		auto vid = ((DeformWidget *)widget)->ui->vertexID->value();
		for (size_t i = 0; i < handles.size(); i++)
			if (!handles[i]->element_id.contains(vid))
				newHandles << handles[i];

		if (newHandles.size() == handles.size()) newHandles.removeLast();

		handles = newHandles;

		if (handles.empty())
		{
			//delete solver; // memory leak
			this->solver = NULL;
			this->isDeformReady = false;
		}

		drawArea()->update();
	});
	connect(dw->ui->createROI, SIGNAL(released()), SLOT(create_ROI()));

	connect(dw->ui->deformButton, &QPushButton::released, [&]{
		if (handles.empty()) return;

		// Options:
		bool is_laplacian = ((DeformWidget *)widget)->ui->laplacianConstraints->isChecked();
		bool is_area = ((DeformWidget *)widget)->ui->areaConstraints->isChecked();
		bool is_volume = ((DeformWidget *)widget)->ui->volumeConstraints->isChecked();
		bool is_surface = !is_volume;
		bool is_dynamic = ((DeformWidget *)widget)->ui->gravityEnabled->isChecked();
		int num_iterations = ((DeformWidget *)widget)->ui->numSolverIterations->value();
		int dimensions = ((DeformWidget *)widget)->ui->dimensions->value();

		// Weights:
		double w1 = ((DeformWidget *)widget)->ui->weight1->value();
		double w2 = ((DeformWidget *)widget)->ui->weight2->value();
		double w3 = ((DeformWidget *)widget)->ui->weight3->value();

		/// 1) Create the solver
		solver = new ShapeOp::Solver;

		/// 2) Set the vertices
		size_t nb_points = mesh()->n_vertices();
		Eigen::Map<ShapeOp::Matrix3X> p(mesh()->vertex_coordinates().data()->data(), 3, nb_points);

		if (is_surface)
		{
			solver->setPoints(p);

			/// 3) Setup the constraints and forces
			// Triangle strain constraints
			{
				double triangle_weight = w1;

				for (auto & face : mesh()->faces())
				{
					std::vector<int> id_vector;
					for (auto & v : mesh()->vertices(face)) id_vector.push_back(v.idx());

					auto c = std::make_shared<ShapeOp::TriangleStrainConstraint>(id_vector, triangle_weight, p);
					solver->addConstraint(c);
				}
			}
			// Bending constraints
			{
				double bending_weight = 1.0;

				for (auto & edge : mesh()->edges())
				{
					if (mesh()->is_boundary(edge)) continue;

					auto h2 = mesh()->halfedge(edge, 1);
					auto h3 = mesh()->halfedge(edge, 0);

					auto v2 = mesh()->to_vertex(h2);
					auto v3 = mesh()->to_vertex(h3);
					auto v1 = mesh()->to_vertex(mesh()->next_halfedge(h3));
					auto v4 = mesh()->to_vertex(mesh()->next_halfedge(h2));

					std::vector<int> id_vector;
					id_vector.push_back(v2.idx());
					id_vector.push_back(v3.idx());
					id_vector.push_back(v1.idx());
					id_vector.push_back(v4.idx());

					auto c = std::make_shared<ShapeOp::BendingConstraint>(id_vector, bending_weight, p);
					solver->addConstraint(c);
				}
			}
			// Area constraints
			if (is_area)
			{
				double area_weight = w2;

				for (auto & face : mesh()->faces())
				{
					std::vector<int> id_vector;
					for (auto & v : mesh()->vertices(face)) id_vector.push_back(v.idx());

					auto c = std::make_shared<ShapeOp::AreaConstraint>(id_vector, area_weight, p);
					solver->addConstraint(c);
				}
			}

			if (is_laplacian)
			{
				double laplacian_weight = w3;

				for (auto & v : mesh()->vertices())
				{
					//if (mesh()->is_boundary(v)) continue;

					std::vector<int> id_vector;
					id_vector.push_back(v.idx());
					for (auto & h : mesh()->onering_hedges(v)) id_vector.push_back(mesh()->to_vertex(h).idx());

					auto c = std::make_shared<ShapeOp::UniformLaplacianConstraint>(id_vector, laplacian_weight, p, false);
					solver->addConstraint(c);
				}
			}
		}

		// Volume constraints
		if (is_volume)
		{
			TetGen tet(mesh());
			auto all_tet_points = tet.getPoints();

			ShapeOp::Matrix3X tetpnts(3, all_tet_points.size());
			auto cells = tet.getCells();

			for (size_t i = 0; i < all_tet_points.size(); i++)
				tetpnts.col(i) = all_tet_points[i];

			// DEBUG:
			if (false){
				for (auto p : all_tet_points){
					auto ps = new starlab::PointSoup;
					ps->addPoint(p);
					drawArea()->addRenderObject(ps);
				}
			}

			solver->setPoints(tetpnts);

			// Volume constraints
			double volume_weight = 1.0;
			for (auto id_vector : cells)
			{
				auto c = std::make_shared<ShapeOp::VolumeConstraint>(id_vector, volume_weight, tetpnts);
				solver->addConstraint(c);

				for (int i = 0; i < 4; i++)
				{
					for (int j = i + 1; j < 4; j++)
					{
						std::vector<int> ids;
						ids.push_back(id_vector[i]);
						ids.push_back(id_vector[j]);
						auto c = std::make_shared<ShapeOp::EdgeStrainConstraint>(ids, 1.0, tetpnts);
						solver->addConstraint(c);
					}
				}
			}
		}

		// Closeness constraints
		{
			double close_weight = 1.0;
			for (auto & handle : handles)
			{
				for (size_t i = 0; i < handle->element_id.size(); i++)
				{
					std::vector<int> id_vector;
					id_vector.push_back(handle->element_id[i]);
					auto c = std::make_shared<ShapeOp::ClosenessConstraint>(id_vector, close_weight, p);
					handle->constraint_id.push_back(solver->addConstraint(c));
				}
			}
		}

		// Forces:
		if (is_dynamic)
		{
			auto f = std::make_shared<ShapeOp::GravityForce>(Vector3(0, 0, -1));
			solver->addForces(f);
		}

		/// 4) Initalize the solver
		solver->initialize(is_dynamic, 0.01, 0.5, 1.0);

		/// 5) Optimize
		solver->solve(num_iterations, dimensions);

		/// 6) Get back the vertices
		auto final_points = solver->getPoints();
		for (size_t i = 0; i < p.cols(); i++) p.col(i) = final_points.col(i);

		this->isDeformReady = true;
		mainWindow()->setStatusBarMessage(QString("Solver ready."));
	});
}

void deform::create_handle(const Vector3 & p, size_t vid)
{
	auto handleRadius = worldRadius * 0.1;
	QSharedPointer<DeformHandle> handle(new DeformHandle(p, handleRadius, ((DeformWidget *)widget)->ui->dimensions->value() == 2));

	handle->element_orig_pos.push_back(p);
	handle->element_id.push_back(vid);

	this->connect(handle.data(), SIGNAL(manipulated()), SLOT(apply_deformation()));
	//drawArea()->setManipulatedFrame(handle.data());
	handles << handle;
	drawArea()->update();
}

void deform::create_ROI()
{
	if (last_selected < 0 || last_selected > mesh()->n_vertices()) return;

	auto vid = Vertex(last_selected);

	// Compute geodesic distance from selected point
	GeoHeatHelper h(mesh());
	ScalarVertexProperty distance = h.getUniformDistance(QSet<Vertex>() << vid);

	// Create a handle of points geodesically within threshold
	auto threshold = ((DeformWidget *)widget)->ui->roiDistance->value();

	QSharedPointer<DeformHandle> handle(new DeformHandle(mesh()->vertex_coordinates()[vid], threshold, ((DeformWidget *)widget)->ui->dimensions->value() == 2));

	for (auto v : mesh()->vertices())
	{
		double dist = distance[v];
		if (dist > threshold) continue;

		handle->element_orig_pos.push_back(mesh()->vertex_coordinates()[v]);
		handle->element_id.push_back(v.idx());
	}

	this->connect(handle.data(), SIGNAL(manipulated()), SLOT(apply_deformation()));
	handles << handle;
	last_selected = -1; // clear selection
	drawArea()->update();
}

void deform::apply_deformation()
{
	if (!solver) return;

	isSolving = true;

	// Update positional constraints
	auto gaussWeight = ((DeformWidget *)widget)->ui->gaussWeight->value();
	for (auto & handle : handles)
	{
		for (size_t i = 0; i < handle->constraint_id.size(); i++)
		{
			auto c = std::dynamic_pointer_cast <ShapeOp::ClosenessConstraint>(solver->getConstraint(handle->constraint_id[i]));
			c->setPosition(handle->transformed(handle->element_orig_pos[i], gaussWeight));
		}
	}

	size_t nb_points = mesh()->n_vertices();
	Eigen::Map<ShapeOp::Matrix3X> p(mesh()->vertex_coordinates().data()->data(), 3, nb_points);
	
	int dimensions = ((DeformWidget *)widget)->ui->dimensions->value();

	solver->solve(((DeformWidget *)widget)->ui->numSolverIterations->value(), dimensions);

	auto final_points = solver->getPoints();
	for (size_t i = 0; i < p.cols(); i++) p.col(i) = final_points.col(i);

	mesh()->update_face_normals();
	//mesh()->update_vertex_normals();

	isSolving = false;
}

void deform::decorate()
{
	double worldRadius = mesh()->bbox().diagonal().norm();
	auto handleRadius = worldRadius * 0.01;

	if (last_selected >= 0)
	{
		starlab::SphereSoup sphere;
		sphere.addSphere(mesh()->vertex_coordinates()[Vertex(last_selected)], handleRadius);
		sphere.draw();
	}

	//drawArea()->drawText(20, drawArea()->height() - 30, QString("Press space for simulation mode."), 10, Qt::gray);

	if (handles.isEmpty()) return;

	// Draw handles
	starlab::FrameSoup fs(handleRadius);
	starlab::PointSoup ps(30);
	for (auto & handle : handles)
	{
		auto handlepos = handle->position();
		auto p = Vector3(handlepos[0], handlepos[1], handlepos[2]);
		ps.addPoint(p, handle->isActive ? Qt::red : Qt::blue);
		fs.addFrame(Vector3(Vector3::UnitX()), Vector3(Vector3::UnitY()), Vector3(Vector3::UnitZ()), p);

		starlab::PointSoup fixed(6);
		if (handle->isActive){
			for (auto vid : handle->element_id){
				fixed.addPoint(mesh()->vertex_coordinates()[Vertex(vid)]);
			}
		}
		glDisable(GL_DEPTH_TEST);
		fixed.draw();
		glEnable(GL_DEPTH_TEST);
	}
	ps.draw();
	fs.draw();
}

void deform::drawWithNames()
{

}

bool deform::postSelection(const QPoint &)
{
	return true;
}

bool deform::keyPressEvent(QKeyEvent *event)
{
	// Simulation mode
	if (event->key() == Qt::Key_Space)
	{
		if (timer){
			timer->deleteLater();
			timer = NULL;
			return true;
		}

		timer = new QTimer(this);
		connect(timer, &QTimer::timeout, [&]() { 
			if (this->isSolving) return;
			this->apply_deformation(); 
			drawArea()->update(); 
		});
		timer->start(30);
		return true;
	}

	return false;
}

bool deform::mouseMoveEvent(QMouseEvent *)
{
	return false;
}

bool deform::mousePressEvent(QMouseEvent * event)
{
	if (event->buttons().testFlag(Qt::MouseButton::RightButton))
	{
		create_ROI();
	}

	if (event->modifiers() & Qt::SHIFT)
	{
		last_selected = -1;

		bool found = false;
		auto pos = drawArea()->camera()->pointUnderPixel(event->pos(), found);
		Vector3 p(pos[0], pos[1], pos[2]);

		// Look up closest handle
		if (isDeformReady)
		{
			if (!found)
			{
				for (auto & h : handles)
					h->isActive = false;
				drawArea()->setManipulatedFrame(drawArea()->camera()->frame());
				return true;
			}

			if (!handles.empty())
			{
				QMap<double, size_t> distMap;
				for (size_t i = 0; i < handles.size(); i++){
					double dist = (p - handles[i]->pos()).norm();
					distMap[dist] = i;
				}
				auto & selectedHandle = handles[distMap.values().front()];
				selectedHandle->isActive = true;
				drawArea()->setManipulatedFrame(selectedHandle.data());
				return true;
			}
		}

		Eigen::Map<ShapeOp::Matrix3X> points(mesh()->vertex_coordinates().data()->data(), 3, mesh()->n_vertices());
		if (found)
		{
			QMap<double, size_t> distMap;
			for (auto v : mesh()->vertices())
			{
				double dist = (p - points.col(v.idx())).norm();
				distMap[dist] = v.idx();
			}
			last_selected = distMap.values().front();
		}

		return true;
	}

	return false;
}
