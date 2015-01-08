#pragma once
#include "SurfaceMeshHelper.h"
#include <qgl.h>

class PointCloud
{
public:
	PointCloud();
	~PointCloud();

	void load(const QString &filename);
	void draw();

	int size() { return m_points.size(); };
	bool hasNormals(){ return !m_normals.empty(); };
	//void resize(int n) { m_points.resize(n); };
	Surface_mesh::Point& operator[](int i) { return m_points[i]; };
	Surface_mesh::Point& points(int i) { return m_points[i]; };
	Surface_mesh::Point& normals(int i) { return m_normals[i]; };
	//std::vector<SurfaceMesh::Point>::iterator begin() { return m_points.begin(); };
	//std::vector<SurfaceMesh::Point>::iterator end() { return m_points.end(); };

private:
	std::vector<SurfaceMesh::Point> m_points;
	std::vector<SurfaceMesh::Point> m_normals;
};

