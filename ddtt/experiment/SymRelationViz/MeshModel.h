#pragma once
#include "SurfaceMeshModel.h"
using namespace SurfaceMesh;
#include "PointCloud.h"
#include <qgl.h>

class MeshModel
{
public:
	MeshModel();
	~MeshModel();
	void loadObj(const QString &filename);
	void loadParts(const QString &pathname);
	void normalization(double scale);
	void translate(double x, double y, double z);
	void draw();
	void buildDisplayList(double ptSize=10.0);
private:
	void clear();
	void clearMesh();
	void clearParts();
public:
	GLuint m_displayListID;
	SurfaceMeshModel *m_mesh;
	QVector<PointCloud*> m_parts;
};

