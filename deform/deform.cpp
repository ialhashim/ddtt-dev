#include "deform.h"
#include "deform-widget.h"
#include "ui_deform-widget.h"
#include "RenderObjectExt.h"

#include "Constraint.h"

int num_iterations = 20;

void deform::create()
{
    if(widget) return;

    // Viewer
    double worldRadius = mesh()->bbox().diagonal().norm();
    {
        //drawArea()->setAxisIsDrawn(true);
        drawArea()->camera()->setType(qglviewer::Camera::PERSPECTIVE);

        drawArea()->camera()->setUpVector(qglviewer::Vec(0,0,1));
        drawArea()->camera()->setPosition(qglviewer::Vec(worldRadius * 2,-worldRadius * 2, worldRadius * 1.5));
        drawArea()->camera()->lookAt(qglviewer::Vec());
        drawArea()->camera()->setSceneRadius( worldRadius );
        drawArea()->camera()->showEntireScene();
    }

    // Rendering
    {
        drawArea()->setRenderer(mesh(), "Flat Wire");
    }

    ModePluginDockWidget * dockwidget = new ModePluginDockWidget("Deform", mainWindow());
    DeformWidget * dw = new DeformWidget();
    widget = dw;

    dockwidget->setWidget( widget );
    mainWindow()->addDockWidget(Qt::RightDockWidgetArea, dockwidget);

    // UI
    connect(dw->ui->createHandle, &QPushButton::released, [&]{
		DeformWidget * dw = (DeformWidget *)widget;
        auto vid = Vertex( dw->ui->vertexID->value() );
        Vector3 p = mesh()->vertex_coordinates()[vid];

        auto handleRadius = worldRadius * 0.2;

        QSharedPointer<DeformHandle> handle( new DeformHandle(p, vid.idx(), handleRadius) );

		this->connect(handle.data(), SIGNAL(manipulated()), SLOT(apply_deformation()));

		drawArea()->setManipulatedFrame( handle.data() );

		handles << handle;
        drawArea()->update();
    });

	connect(dw->ui->deformButton, &QPushButton::released, [&]{
		if ( handles.empty() ) return;

		// 1) Create the solver
		solver = new ShapeOp::Solver;
		
		// 2) Set the vertices
		size_t nb_points = mesh()->n_vertices();
		Eigen::Map<ShapeOp::Matrix3X> p(mesh()->vertex_coordinates().data()->data(), 3, nb_points);
		solver->setPoints(p);

		// 3) Setup the constraints and forces

		// Edge strain constraints
		if (false)
		{
			double edge_weight = 1.0;

			for (auto & edge : mesh()->edges())
			{
				auto v1 = mesh()->vertex(edge, 0);
				auto v2 = mesh()->vertex(edge, 1);

				std::vector<int> id_vector;
				id_vector.push_back(v1.idx());
				id_vector.push_back(v2.idx());

				auto c = std::make_shared<ShapeOp::EdgeStrainConstraint>(id_vector, edge_weight, p);
				solver->addConstraint(c);
			}
		}

		// Triangle strain constraints
		if (true)
		{
			double triangle_weight = 1.0;

			for (auto & face : mesh()->faces())
			{
				std::vector<int> id_vector;
				for (auto & v : mesh()->vertices(face)) id_vector.push_back(v.idx());

				auto c = std::make_shared<ShapeOp::TriangleStrainConstraint>(id_vector, triangle_weight, p);
				solver->addConstraint(c);
			}
		}

		// Laplacian constraints
		if (true)
		{
			double laplacian_weight = 1.0;

			for (auto & v : mesh()->vertices())
			{
				std::vector<int> id_vector;
				id_vector.push_back(v.idx());
				for (auto hi : mesh()->onering_hedges(v)){ 
					auto vj = mesh()->to_vertex(hi);
					id_vector.push_back(vj.idx());
				}

				auto c = std::make_shared<ShapeOp::UniformLaplacianConstraint>(id_vector, laplacian_weight, p, true);
				solver->addConstraint(c);
			}
		}

		// Closeness constraints
		{
			double close_weight = 1.0;
			for (auto & handle : handles)
			{
				std::vector<int> id_vector;
				id_vector.push_back( handle->id );
				auto c = std::make_shared<ShapeOp::ClosenessConstraint>(id_vector, close_weight, p);
				c->setPosition( handle->pos() );

				handle->constraint_id = solver->addConstraint(c);
			}
		}

		// Forces:
		//

		// 4) Initalize the solver
		solver->initialize();

		// 5) Optimize
		solver->solve( num_iterations );

		// 6) Get back the vertices
		p = solver->getPoints();

		mainWindow()->setStatusBarMessage(QString("Solver created."));
	});
}

void deform::apply_deformation()
{
	if (!solver) return;

	// Update constraints
	for (auto & handle : handles)
	{
		auto c = std::dynamic_pointer_cast < ShapeOp::ClosenessConstraint >(solver->getConstraint(handle->constraint_id));
		c->setPosition( handle->pos() );
	}

	size_t nb_points = mesh()->n_vertices();
	Eigen::Map<ShapeOp::Matrix3X> p(mesh()->vertex_coordinates().data()->data(), 3, nb_points);

	solver->solve( num_iterations );

	p = solver->getPoints();

	mesh()->update_face_normals();
	mesh()->update_vertex_normals();
}

void deform::decorate()
{
	if (handles.isEmpty()) return;
	
	double worldRadius = mesh()->bbox().diagonal().norm();

	// Draw handles
	auto handleRadius = worldRadius * 0.2;
	starlab::FrameSoup fs(handleRadius);
	starlab::PointSoup ps (20);
	for (auto & handle : handles)
	{
		auto handlepos = handle->position();
		auto p = Vector3(handlepos[0], handlepos[1], handlepos[2]);
		ps.addPoint(p);
		fs.addFrame(Vector3(Vector3::UnitX()), Vector3(Vector3::UnitY()), Vector3(Vector3::UnitZ()), p);
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

bool deform::keyPressEvent(QKeyEvent *)
{
    return false;
}

bool deform::mouseMoveEvent(QMouseEvent *)
{
    return false;
}

bool deform::mousePressEvent(QMouseEvent * event)
{
    if (event->modifiers() & Qt::SHIFT)
    {
        drawArea()->select(event->pos());
        return true;
    }

    return false;
}
