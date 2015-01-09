#include "PointCloud.h"
#include <qfile.h>

PointCloud::PointCloud()
{
}


PointCloud::~PointCloud()
{
}

void PointCloud::load(const QString &filename)
{
	QFile file(filename);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		QString str("can't open ");
		str.append(filename);
		QMessageBox::warning(0, "Warnning", str, QMessageBox::Yes);
		return;
	}
	m_points.clear(); m_normals.clear();

	QTextStream in(&file);
	QString str;
	for (int i = 0; i < 4; ++i)
	{
		str = in.readLine();
	}
	float x,y,z,nx,ny,nz;
	while (!in.atEnd())
	{
		in >> x >> y >> z >> nx >> ny >> nz;
		m_points.push_back(Surface_mesh::Point(x, y, z));
		m_normals.push_back(Surface_mesh::Point(nx, ny, nz));
	}


	file.close();
}
void PointCloud::updateBoundingBox()
{
	m_bbox.setNull();
	for (std::vector<SurfaceMesh::Point>::iterator it = m_points.begin(); it != m_points.end(); ++it)
	{
		m_bbox.extend(*it);
	}
}