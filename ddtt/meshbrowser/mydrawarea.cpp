#include "mydrawarea.h"
#include "surface_mesh/gl_wrappers.h"
#include <QMouseEvent>

QString curLabel = "";
QStringList AllLabels;
QVector<QColor> UniqueColors;
QVector<QString> labelNames; 

using namespace SurfaceMesh;

void MyDrawArea::draw()
{
	SurfaceMesh::SurfaceMeshModel * mesh = m;

	if(!mesh) return;

	if(!mesh->property("hasNormals").toBool())
	{
		mesh->update_face_normals();
		mesh->update_vertex_normals();
		mesh->setProperty("hasNormals",true);
	}

	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LIGHTING);

	

	Surface_mesh::Vertex_property<Vector3> points = mesh->vertex_property<Vector3>("v:point");
	Surface_mesh::Face_property<Vector3> fnormals = mesh->face_property<Vector3>("f:normal");

	Surface_mesh::Face_iterator fit, fend = mesh->faces_end();
	Surface_mesh::Vertex_around_face_circulator fvit, fvend;

    glColor3f(1,1,1);
	glBegin(GL_TRIANGLES);
	for (fit=mesh->faces_begin(); fit!=fend; ++fit){
		glNormal3d( fnormals[fit][0], fnormals[fit][1], fnormals[fit][2] );
		fvit = fvend = mesh->vertices(fit);
		do{ glVertex3d(points[fvit][0],points[fvit][1],points[fvit][2]); } while (++fvit != fvend);
	}
	glEnd();

    /// Render wireframe
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    Surface_mesh::Edge_iterator eit, eend=m->edges_end();
    for (eit=m->edges_begin(); eit!=eend; ++eit){
        Surface_mesh::Edge e = eit;
        Surface_mesh::Vertex v0 = m->vertex(eit,0);
        Surface_mesh::Vertex v1 = m->vertex(eit,1);
        gl::glColor(Eigen::Vector4d(0,0,0,0.5));
        gl::glVertex(points[v0]);
        gl::glVertex(points[v1]);
    }
    glEnd();
    glEnable(GL_LIGHTING);
}

void MyDrawArea::mousePressEvent(QMouseEvent * event)
{
	QGLViewer::mousePressEvent(event);

    if(!(event->buttons() & Qt::RightButton)) return;

    Surface_mesh::Vertex_property<Vector3> points = m->vertex_property<Vector3>("v:point");

    if(curOp == MeshOperation::ROTATE_LEFT){
        Eigen::AngleAxisd rot(0.5 * M_PI, Vector3::UnitZ());
        for(auto v : m->vertices())
            points[v] = rot * points[v];
    }

    if(curOp == MeshOperation::ROTATE_RIGHT){
        Eigen::AngleAxisd rot(-0.5 * M_PI, Vector3::UnitZ());
        for(auto v : m->vertices())
            points[v] = rot * points[v];
    }

    m->updateBoundingBox();
    m->update_face_normals();
    m->update_vertex_normals();
}

void MyDrawArea::postSelection(const QPoint&)
{

}

MyDrawArea::~MyDrawArea()
{
    m->write( qPrintable(filename) );
}
