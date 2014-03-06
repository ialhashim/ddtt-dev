#include <QElapsedTimer>
#include <QFileDialog>

#include "RenderObjectExt.h"

#include "voxelize.h"
#include "voxelization.h"

void voxelize::initParameters(RichParameterSet *pars)
{
	pars->addParam(new RichInt("Resolution", 32, "Resolution"));
	pars->addParam(new RichBool("Intersection", true, "Intersection"));
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
		Vector3VertexProperty points = mesh()->vertex_coordinates();
		Vector3 corner = mesh()->bbox().min();
		Vector3 delta = mesh()->bbox().center() - corner;

		for(auto v : mesh()->vertices()){
			points[v] -= corner;
		}

		QElapsedTimer timer; timer.start();

		vector<VoxelData> data = ComputeVoxelization( mesh(), gridsize );

		mainWindow()->setStatusBarMessage( QString("Time (%1 ms)").arg(timer.elapsed()) );

		Eigen::AlignedBox3d mesh_bbox = mesh()->bbox();
		double unitlength = (mesh_bbox.max()[0] - mesh_bbox.min()[0]) / (double)gridsize;

		starlab::CubeSoup * cs = new starlab::CubeSoup;

		std::vector<Vector3> centers;

		for(auto voxel : data)
		{
			unsigned int x,y,z;

			mortonDecode(voxel.morton, z, y, x);

			Vector3 p = Vector3(Vector3(x,y,z) * unitlength) + ( 0.5 * Vector3(unitlength,unitlength,unitlength) );
			p += corner;

			centers.push_back( p );
		}

		for(auto v : mesh()->vertices()){
			points[v] += corner;
		}

		if(pars->getBool("Visualize"))
			for(auto p : centers) cs->addCube( p, unitlength );

		drawArea()->addRenderObject( cs );
	}
}
