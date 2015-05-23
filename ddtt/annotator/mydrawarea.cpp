#include "mydrawarea.h"
#include "surface_mesh/gl_wrappers.h"
#include <QMouseEvent>
#include <QDir>
#include <QFileInfo>
#include <stack>
#include "Octree.h"

#include "SurfaceMeshHelper.h"
using namespace SurfaceMesh;

QStringList AllLabels;
QVector<QColor> LabelColors; 
QVector<QColor> ParentColors;
QVector<QString> labelNames; 

Q_DECLARE_METATYPE(std::vector< std::vector<double> >)

#pragma warning(disable:4309)
#pragma warning(disable:4267)

MyDrawArea::MyDrawArea(QSharedPointer<Structure::Graph> shape, QString filename) : m(shape), filename(filename)
{
	isDeleted = false;
	isDrawWireframe = true;
	isDoubleLight = false;

	basename = QFileInfo(filename).baseName();

	bg = QColor(51, 51, 51);
	fg = QColor(255, 255, 255);

	// default options
	shape->property["showNodes"] = false;
	shape->property["showMeshes"] = true;
	int ni = 0;
	for (auto n : shape->nodes)
	{
		auto c = LabelColors[ni++];
		c.setAlphaF(0.3);
		shape->setColorFor(n->id, c.darker(300));
		n->vis_property["meshSolid"].setValue(true);
	}
	for (auto group : shape->groups)
	{
		auto repNode = group.front();
		auto c = shape->getNode(repNode)->vis_property["color"].value<QColor>();

		for (auto element : group)
		{
			shape->setColorFor(element, c);
		}
	}

	// Previously labeled
	for (auto n : shape->nodes)
	{
		if (n->meta.contains("label"))
		{
			auto label = n->meta["label"].toString();
			if (!labelNames.contains(label))
			{
				n->meta.remove("label");
				n->meta.remove("labelParent");
				continue;
			}

			n->property["label"].setValue(label);
			shape->setColorFor(n->id, LabelColors[labelIndices[label]]);
		}
	}
}

void MyDrawArea::draw()
{
    if(m.isNull()) return;

	if(isDoubleLight){
		glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	}

	glClearColor( bg.redF(), bg.greenF(), bg.blueF(), bg.alphaF() );
	glClear(GL_COLOR_BUFFER_BIT);

	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LIGHTING);

    m->draw(this);

	if( isDeleted )
	{
		this->startScreenCoordinatesSystem();
		glLineWidth(5);

		glColor3f(1,0,0); //red
		glBegin(GL_LINES);
		gl::glVertex(Vector3(0,0,0));
		gl::glVertex(Vector3(this->width(), this->height(),0));
		glEnd();

		glLineWidth(2);
		this->stopScreenCoordinatesSystem();
    }

    glPushAttrib(GL_ALL_ATTRIB_BITS);
    startScreenCoordinatesSystem();
    glColor3d(0,0,0);
	renderText(11, 21, basename);
	glColor3d(1,1,1);
	renderText(10, 20, basename);

	int numTaged = 0;
	for (auto n : m->nodes) if (n->property.contains("label")) numTaged++;
	bool isTagFull = numTaged == m->nodes.size();
	if (isTagFull) glColor3d(0.8, 1, 0.8); else glColor3d(1.0, 0.7, 0.7);
	QString tagString = QString("tags %1/%2").arg(numTaged).arg(m->nodes.size());
	renderText(width() - QFontMetrics(QFont()).width(tagString) - 10, 20, tagString);
    stopScreenCoordinatesSystem();
    glPopAttrib();

    glEnable( GL_MULTISAMPLE );
}

struct IsectData{
	double u,v,t;
};

bool rayTriangleIntersectionTest( const Ray &ray, Vector3 v0, Vector3 v1, Vector3 v2, HitResult & res, bool allowBack = true )
{
	res.hit = false;
	res.distance = DBL_MAX;

	double EPS = 1e-7;

	Eigen::Vector3d vertex1 = v0;
	Eigen::Vector3d vertex2 = v1;
	Eigen::Vector3d vertex3 = v2;

	// Compute vectors along two edges of the triangle.
	Eigen::Vector3d edge1 = vertex2 - vertex1;
	Eigen::Vector3d edge2 = vertex3 - vertex1;

	// Compute the determinant.
	Eigen::Vector3d directionCrossEdge2 = cross(ray.direction, edge2);

	double determinant = dot(edge1, directionCrossEdge2);

	// If the ray is parallel to the triangle plane, there is no collision.
	if (fabs(determinant) < EPS)
		return false;

	double inverseDeterminant = 1.0 / determinant;

	// Calculate the U parameter of the intersection point.
	Eigen::Vector3d distVector = ray.origin - vertex1;
	double triangleU = dot(distVector, directionCrossEdge2);
	triangleU *= inverseDeterminant;

	// Make sure it is inside the triangle.
	if (triangleU < 0 - EPS || triangleU > 1 + EPS)
		return false;

	// Calculate the V parameter of the intersection point.
	Eigen::Vector3d distanceCrossEdge1 = cross(distVector, edge1);
	double triangleV = dot(ray.direction, distanceCrossEdge1);
	triangleV *= inverseDeterminant;

	// Make sure it is inside the triangle.
	if (triangleV < 0 - EPS || triangleU + triangleV > 1 + EPS)
		return false;

	// Compute the distance along the ray to the triangle.
	double rayDistance = dot(edge2, distanceCrossEdge1);
	rayDistance *= inverseDeterminant;

	if(!allowBack){
		// Is the triangle behind the ray origin?
		if (rayDistance < 0)
			return false;
	}

	res.hit = true;
	res.distance = rayDistance;

	res.u = triangleU;
	res.v = triangleV;
	return true;
}


void MyDrawArea::mousePressEvent(QMouseEvent * event)
{
	QGLViewer::mousePressEvent(event);

    if(!(event->buttons() & Qt::RightButton)) return;

    if(curOp == ShapeOperation::ROTATE_LEFT){
        Eigen::AngleAxisd rot(0.5 * M_PI, Vector3::UnitZ());
        //for(auto v : m->vertices())
        //    points[v] = rot * points[v];
    }

    if(curOp == ShapeOperation::ROTATE_RIGHT){
        Eigen::AngleAxisd rot(-0.5 * M_PI, Vector3::UnitY());
        //for(auto v : m->vertices())
        //    points[v] = rot * points[v];
    }

	if (curOp == ShapeOperation::LABEL_PART){
		qglviewer::Vec _orig, _dir;
		camera()->convertClickToLine(event->pos(), _orig, _dir);
		Vector3 orig(_orig[0], _orig[1], _orig[2]);
		Vector3 dir(_dir[0], _dir[1], _dir[2]);
		Ray r(orig, dir);

		if (curLabel < 0) return;

		for (auto n : m->nodes)
		{
			auto mesh = m->getMesh(n->id);	
			bool isHitMesh = false;
			
			Vector3VertexProperty points = mesh->vertex_coordinates();

			for (auto f : mesh->faces()){
				std::vector<Vector3> v;
				for (auto vert : mesh->vertices(f))
					v.push_back(points[vert]);

				HitResult hit;
				if (rayTriangleIntersectionTest(r, v[0], v[1], v[2], hit))
				{
					isHitMesh = true;
					break;
				}
			}

			if (isHitMesh)
			{
				n->property["label"].setValue(curLabel);
				m->setColorFor(n->id, LabelColors[curLabelIdx]);
				break;
			}
		}

	}
}

void MyDrawArea::postSelection(const QPoint&)
{

}

MyDrawArea::~MyDrawArea()
{
	if(isDeleted) return;

	for (auto n : m->nodes)
	{
		if (!n->property.contains("label")) continue;

		QString nodeLabel = n->property["label"].toString();
		n->meta["label"].setValue(nodeLabel);
		n->meta["parentLabel"].setValue(nodeLabel.split("-").front());
	}

	m->saveToFile(filename);
}

void MyDrawArea::focusInEvent(QFocusEvent * event)
{
	QGLViewer::focusInEvent(event);

	emit( gotFocus(this) );
}
