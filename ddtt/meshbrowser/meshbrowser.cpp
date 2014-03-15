#include "meshbrowser.h"
#include "ui_meshbrowser.h"
#include "StarlabDrawArea.h"
#include <QPushButton>
#include <QFileDialog>

#include "mydrawarea.h"
MeshOperation curOp = ROTATE_LEFT;

using namespace SurfaceMesh;

MeshBrowser::MeshBrowser(QWidget *parent) : QWidget(parent), ui(new Ui::MeshBrowser)
{
	// Anti-aliasing when using QGLWidget or subclasses
	QGLFormat glf = QGLFormat::defaultFormat();
	glf.setSamples(8);
	QGLFormat::setDefaultFormat(glf);

    ui->setupUi(this);

    {
        QString datasetPath = QFileDialog::getExistingDirectory(0, "Dataset");

        // Get list of files
        QStringList filters;
        filters << "*.obj" << "*.off";
        QStringList files = QDir(datasetPath).entryList(filters);

        foreach(QString mesh, files)
        {
            database << (datasetPath + "/" + mesh);
		}
	}

	// Create viewers
	{
		nU = 2;
		nV = 3;

		// Scroll bars
		int numPages = ceil(double(database.size()) / double(nU * nV)) - 1;
		ui->viewersScroll->setValue(0);
		ui->viewersScroll->setRange(0, numPages);
		connect(ui->viewersScroll, SIGNAL(valueChanged(int)), SLOT(loadMeshes()));

		// Load meshes
		loadMeshes();
	}

    // Connections
    connect(ui->rotateLeft, &QPushButton::released, [=]() {curOp = MeshOperation::ROTATE_LEFT; ui->curLabel->setText(meshOps[curOp]);});
    connect(ui->rotateRight, &QPushButton::released, [=]() {curOp = MeshOperation::ROTATE_RIGHT; ui->curLabel->setText(meshOps[curOp]);});
    connect(ui->flipNormals, &QPushButton::released, [=]() {curOp = MeshOperation::FLIP; ui->curLabel->setText(meshOps[curOp]);});
    connect(ui->invertPart, &QPushButton::released, [=]() {curOp = MeshOperation::INVERT_PART; ui->curLabel->setText(meshOps[curOp]);});
}

void MeshBrowser::loadMeshes()
{
	QGridLayout * layout = (QGridLayout *)ui->viewersFrame->layout();
	QLayoutItem* item;
	while ( ( item = layout->takeAt( 0 ) ) != NULL ){
		delete item->widget(); delete item;
	}

	int offset = ui->viewersScroll->value() * (nU * nV);

	for(int i = 0; i < nU; i++){
		for(int j = 0; j < nV; j++){

			int idx = offset + (i*nV) + j;
			if(idx + 1 > database.size()) continue;
			QString filename = database[idx];

            //Structure::Graph * g = new Structure::Graph( filename );
            //g->property["showNodes"] = false;

            SurfaceMesh::SurfaceMeshModel * m = new SurfaceMeshModel(filename, QFileInfo(filename).baseName());
            m->read( qPrintable(filename) );

            Vector3VertexProperty points = m->vertex_coordinates();

            /// Normalize, center, and move to base

            // Center to zero point
            Vector3 mean(0,0,0);
            foreach(Vertex v, m->vertices()) mean += points[v];
            mean /= m->n_vertices();
            foreach(Vertex v, m->vertices()) points[v] -= mean;

            // Scale maximum dimension to 1.0
            m->updateBoundingBox();
            Eigen::AlignedBox3d orig_bbox = m->bbox();
            Vector3 d = orig_bbox.sizes();
            double s = (d.x() > d.y())? d.x():d.y();
            s = (s>d.z())? s : d.z();
            foreach(Vertex v, m->vertices()) points[v] /= s;

            // Move to floor
            double minZ = DBL_MAX;
            foreach(Vertex v, m->vertices()) minZ = qMin(minZ, points[v].z());
            foreach(Vertex v, m->vertices()) points[v] -= Vector3(0,0,minZ);

            m->updateBoundingBox();
            m->update_face_normals();
            m->update_vertex_normals();

            MyDrawArea * viewer = new MyDrawArea(m, filename);
			layout->addWidget(viewer, i, j, 1, 1);

			viewer->setGridIsDrawn();

			Eigen::AlignedBox3d bbox = m->bbox();

            //viewer->setGridIsDrawn();
			viewer->camera()->setSceneRadius(m->bbox().sizes().norm() * 2);
			viewer->camera()->setUpVector(qglviewer::Vec(0,0,1));
			viewer->camera()->setPosition(qglviewer::Vec(2,-2,1.5));
			viewer->camera()->lookAt(qglviewer::Vec());
			viewer->camera()->showEntireScene();

            double S = bbox.sizes().norm() * 0.1;
            bbox.extend(bbox.min() + (Vector3(-1,-1,-1) * S));
            bbox.extend(bbox.max() + (Vector3(1,1,1) * S));

            viewer->camera()->setRevolveAroundPoint(qglviewer::Vec(bbox.center()));
            viewer->camera()->setSceneCenter(qglviewer::Vec(bbox.center()));
			viewer->camera()->fitBoundingBox(qglviewer::Vec(bbox.min().data()),qglviewer::Vec(bbox.max().data()));
		}
	}

	this->activateWindow();
	this->setFocus();
	this->raise();
}

MeshBrowser::~MeshBrowser()
{
    delete ui;
}
