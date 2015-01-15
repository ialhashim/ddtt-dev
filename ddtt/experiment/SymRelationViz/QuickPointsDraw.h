#pragma once
#include "PointCloud.h"
#include <qgl.h>

#define glVertQt(v) glVertex3d(v.x(), v.y(), v.z())
#define glColorQt(c) glColor4d(c.redF(), c.greenF(), c.blueF(), c.alphaF())

#include <time.h>
static std::vector<double> randomColor()
{
	std::vector<double> color;

	float r = float(qMin(((rand() % 225) + 30), 255)) / 255.0f;
	float g = float(qMin(((rand() % 230) + 25), 255)) / 255.0f;
	float b = float(qMin(((rand() % 235) + 20), 255)) / 255.0f;

	color.push_back(r);
	color.push_back(g);
	color.push_back(b);
	color.push_back(1.0);

	return color;
}

static std::vector< std::vector<double> > randomColors(int count)
{
	srand(time(NULL));

	std::vector< std::vector<double> > colors(count);
	for (int i = 0; i < count; i++)
		colors[i] = randomColor();
	return colors;
}

static QColor qRandomColor()
{
	std::vector<double> c = randomColor();
	return QColor::fromRgbF(c[0], c[1], c[2], c[3]);
}

struct QuickPointsDraw
{
	static void drawPoints(PointCloud * pts, float size = 1.0, QColor c = QColor(255, 255, 255, 255))
	{
		if (!pts) return;

		glPointSize(size);
		glColorQt(c);
		glDisable(GL_LIGHTING);
		//if (pts->hasNormals())
		//{
		//	glEnable(GL_LIGHTING);
		//	glBegin(GL_POINTS);
		//	for (int i = 0; i < pts->size(); ++i)
		//	{
		//		SurfaceMesh::Point& n = pts->normals(i);
		//		glNormal3d(n.x(), n.y(), n.z());
		//		glVertQt(pts->points(i));
		//	}
		//	glEnd();
		//}
		//else
		//{
		    //glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glBegin(GL_POINTS);
			for (int i = 0; i < pts->size(); ++i)
			{
				glVertQt(pts->points(i));
			}
			glEnd();
		//}

		glEnable(GL_LIGHTING);
	}
};