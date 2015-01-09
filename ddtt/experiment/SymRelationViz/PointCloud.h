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
	const Eigen::AlignedBox3d& bbox(){ return m_bbox; }
	void updateBoundingBox();

private:
	std::vector<SurfaceMesh::Point> m_points;
	std::vector<SurfaceMesh::Point> m_normals;
	Eigen::AlignedBox3d m_bbox;
};

