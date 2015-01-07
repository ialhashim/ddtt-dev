#pragma once
#include "SurfaceMeshModel.h"
using namespace SurfaceMesh;
#include <qgl.h>

class MeshModel
{
public:
	MeshModel();
	~MeshModel();
	void loadObj(QString filename);
	void normalization(double scale);
	void translate(double x, double y, double z);
	void draw();
	void buildDisplayList();


	GLuint m_displayListID;
	SurfaceMeshModel *m_mesh;
};

