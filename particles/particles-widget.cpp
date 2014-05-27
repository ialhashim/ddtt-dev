#include <QPushButton>
#include <QFileDialog>
#include "particles-widget.h"
#include "ui_particles-widget.h"

#include "SurfaceMeshModel.h"
using namespace SurfaceMesh;

#include "ParticleMesh.h"
#include "voxelization.h"

ParticlesWidget::ParticlesWidget(QWidget *parent) : QWidget(parent), ui(new Ui::ParticlesWidget), isReady(false)
{
    ui->setupUi(this);

	int gridsize = 128;

	connect(ui->loadShapes, &QPushButton::released, [=]{
		QStringList files = QFileDialog::getOpenFileNames(nullptr, "Open Shapes", "", "All Supported (*.obj *.off)");
		for(auto filename : files)
		{
			SurfaceMeshModel mesh;
			mesh.read( filename.toStdString() );

			std::vector<Eigen::Vector3d> mesh_points;
			VoxelContainer container = ComputeVoxelization(&mesh, gridsize, true, false);
			double unitlength = container.unitlength;

			for(auto voxel : container.data)
			{
				unsigned int x,y,z;
				mortonDecode(voxel.morton, z, y, x);
				Vector3 p = container.translation + Vector3(Vector3(x,y,z) * unitlength) + ( 0.5 * Vector3(unitlength,unitlength,unitlength) );
				
				mesh_points.push_back( p );
			}

			pmeshes.push_back( new ParticleMesh( mesh_points, 0.01 ) );
		}

		// Match particles
		for(size_t i = 0; i < pmeshes.size(); i++){
			for(size_t j = i+1; j < pmeshes.size(); j++){
				NanoKdTree * itree = pmeshes[i]->kdtree;
				NanoKdTree * jtree = pmeshes[j]->kdtree;

				for(auto & iparticle : pmeshes[i]->particles){
					iparticle.correspondence = jtree->closest( iparticle.pos );
				}

				for(auto & jparticle : pmeshes[j]->particles){
					jparticle.correspondence = itree->closest( jparticle.pos );
				}
			}
		}

		isReady = true;
	});
}

ParticlesWidget::~ParticlesWidget()
{
    delete ui;
}
