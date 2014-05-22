#include <QElapsedTimer>
#include <QFileDialog>

#include "RenderObjectExt.h"

#include "voxelize.h"
#include "voxelization.h"

void voxelize::initParameters(RichParameterSet *pars)
{
	pars->addParam(new RichInt("Resolution", 128, "Resolution"));
	pars->addParam(new RichBool("Solid", false, "Solid"));
	pars->addParam(new RichBool("Intersection", false, "Intersection"));
	pars->addParam(new RichInt("Operation", 0, "Operation"));
	pars->addParam(new RichBool("Visualize", true, "Visualize"));
}

void voxelize::applyFilter(RichParameterSet *pars)
{
	drawArea()->clear();

	drawArea()->setAxisIsDrawn();

	size_t gridsize = pars->getInt("Resolution");

	if(pars->getBool("Intersection"))
	{
		QList<Starlab::Model*> models = mainWindow()->document()->models();

		if(models.size() < 2){
			foreach( QString filename, QFileDialog::getOpenFileNames(0, "Open Meshes", 
				mainWindow()->settings()->getString("lastUsedDirectory"), "Mesh Files (*.obj)") )
			{
				SurfaceMeshModel * m = new SurfaceMeshModel(filename, QFileInfo(filename).baseName());
				m->read( filename.toStdString() );
				m->updateBoundingBox();
				m->update_face_normals();
				m->update_vertex_normals();

				document()->addModel(m);
			}

			models = mainWindow()->document()->models();
		}

		SurfaceMeshModel * modelA = dynamic_cast<SurfaceMeshModel*>(models.at(models.size() - 2));
		SurfaceMeshModel * modelB = dynamic_cast<SurfaceMeshModel*>(models.at(models.size() - 1));

		double unitlength = 0;
		vector<VoxelData> data = ComputeVoxelizationCSG(modelA, modelB, unitlength, gridsize, BooleanOperation(pars->getInt("Operation")));

		starlab::CubeSoup * cs = new starlab::CubeSoup;
		starlab::PointSoup * ps = new starlab::PointSoup(1);

		for(auto voxel : data)
		{
			unsigned int x,y,z;
			mortonDecode(voxel.morton, z, y, x);
			Vector3 p = Vector3(Vector3(x,y,z) * unitlength) + ( 0.5 * Vector3(unitlength,unitlength,unitlength) );

			if(pars->getBool("Visualize"))
			{
				if(gridsize < 64)
					cs->addCube(p, unitlength);
				else
					ps->addPoint(p, Qt::yellow);
			}
		}

		drawArea()->addRenderObject( cs );
		drawArea()->addRenderObject( ps );
	}
	else
	{
		QElapsedTimer timer; timer.start();

		double unitlength = 1.0;
		VoxelContainer voxels = ComputeVoxelization( mesh(), unitlength, gridsize, pars->getBool("Solid") );

		mainWindow()->setStatusBarMessage( QString("Time (%1 ms) / Count (%2 voxels)").arg(timer.elapsed()).arg(voxels.data.size()));

		starlab::CubeSoup * cs = new starlab::CubeSoup(1.0, false);

		std::vector<Vector3> centers;

		for(auto voxel : voxels.data)
		{
			unsigned int x,y,z;

			mortonDecode(voxel.morton, z, y, x);

			Vector3 p = voxels.translation + Vector3(Vector3(x,y,z) * unitlength) + ( 0.5 * Vector3(unitlength,unitlength,unitlength) );
			centers.push_back( p );
		}

		if(pars->getBool("Visualize"))
			for(auto p : centers) cs->addCube( p, unitlength );

		drawArea()->addRenderObject( cs );
	}
}
