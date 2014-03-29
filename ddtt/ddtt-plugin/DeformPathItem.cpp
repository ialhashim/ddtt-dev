#include "ShapeCorresponder.h"
#include "DeformPathItem.h"
#include <QGraphicsScene>
#include <QGraphicsView>

DeformPathItem::DeformPathItem(DeformationPath *usedPath) : path(usedPath){
	m_rect = QRectF(0,0,800,600);
}

QRectF DeformPathItem::boundingRect() const{
	return m_rect;
}

void qgluPerspective(GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar)
{
	const GLdouble ymax = zNear * tan(fovy * M_PI / 360.0);
	const GLdouble ymin = -ymax;
	const GLdouble xmin = ymin * aspect;
	const GLdouble xmax = ymax * aspect;
	glFrustum(xmin, xmax, ymin, ymax, zNear, zFar);
}

void setCamera()
{
	qgluPerspective(60, 1.0, 0.1, 40);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glTranslatef(0, 0, -3);
	glRotatef(50, 1, 0, 0);
	glRotatef(10, 0, 0, 1);
}

QRectF setCamera2(QRectF boundingRect, Structure::Graph * g)
{
	qglviewer::Camera * camera = new qglviewer::Camera;

	camera->setUpVector(qglviewer::Vec(0,0,1));
	camera->setPosition(qglviewer::Vec(2,-2,1.5));
	camera->lookAt(qglviewer::Vec());
	camera->setSceneRadius( 10 );
	camera->showEntireScene();

	if(camera->type() != qglviewer::Camera::PERSPECTIVE) camera->setType(qglviewer::Camera::PERSPECTIVE);

	QRectF r = boundingRect;

	qglviewer::Vec viewDir = camera->viewDirection();
	Eigen::AlignedBox3d graphBBox = g->bbox();
	double distance = graphBBox.sizes().maxCoeff() * 2.25;
	Vector3 center = graphBBox.center();
	Vector3 newPos = center - (distance * Vector3(viewDir[0], viewDir[1], viewDir[2]));
	camera->setRevolveAroundPoint( qglviewer::Vec(center) );
	qglviewer::Vec new_pos(newPos);
	camera->frame()->setPositionWithConstraint(new_pos);
	camera->setScreenWidthAndHeight(r.width(), r.height());
	camera->loadProjectionMatrix();
	camera->loadModelViewMatrix();

	return r;
}

void setLights()
{
	glEnable(GL_LIGHTING);
	GLfloat lightColor[] = {0.9f, 0.9f, 0.9f, 1.0f};
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor);
	glEnable(GL_LIGHT0);

	// Material
	float mat_ambient[] = {0.1745f, 0.01175f, 0.01175f, 1.0f};
	float mat_diffuse[] = {0.65f, 0.045f, 0.045f, 1.0f};
	float mat_specular[] = {0.09f, 0.09f, 0.09f, 1.0f};
	float high_shininess = 100;

	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialf(GL_FRONT, GL_SHININESS, high_shininess);
}

void DeformPathItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *, QWidget *)
{
	int viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	QGraphicsView * view = scene()->views().front();
	QPointF sceneP = mapToScene( QPointF(0,0) );
	QPoint viewP = view->mapFromScene(sceneP);

    double dy = -viewP.y();
    double dx = viewP.x();

    if( dy > viewport[3] || dy < -viewport[3] || dx > viewport[2] || dx < -viewport[2] ) return;

    QRectF m_rect_scaled = view->transform().mapRect(m_rect);
    int width = m_rect_scaled.width();
    int height = m_rect_scaled.height();
    int smallWidth = height * 0.5;
    int tinyWidth = height * (1.0 / 8.0);

    QRectF srect(width - smallWidth ,0,smallWidth,smallWidth);
    QRectF trect(width - smallWidth,smallWidth,smallWidth,smallWidth);
    QRectF inbetween(0, tinyWidth, width - smallWidth, height - tinyWidth);

	// DEBUG:
	if( false ){
		painter->drawText(QPoint(10,200), QString("IDX %3  = %1, %2, dy = %4").arg(viewP.x()).arg(viewP.y()).arg(path->idx).arg(dy));
		painter->setPen( QPen(Qt::green, 2) );	painter->drawRect( srect );
		painter->setPen( QPen(Qt::blue, 2) );	painter->drawRect( trect );
		painter->setPen( QPen(Qt::black, 2) );	painter->drawRect( inbetween );
	}

	// Draw score
	painter->drawText( QPoint(10,inbetween.top() + 20), QString("%1").arg(path->weight) );

	// Draw border
    painter->setPen(QPen(Qt::gray, 1));
    if(this->parentItem()->isSelected()) painter->setPen(QPen(Qt::blue, 5));
    painter->drawRect(boundingRect());

	// Draw 3D parts
	painter->beginNativePainting();
	setLights();
	setCamera2(srect, path->gcorr->sg);

	// Draw source
	{
		Structure::Graph * g = path->gcorr->sg;

		glViewport(dx + srect.x(), dy + viewport[3] - srect.height() - srect.top(), srect.width(), srect.height());

		for(auto n : g->nodes) 
		{
			glColor3d(0,0,0);
			if(path->scolors.contains(n->id)) glColorQt(path->scolors[n->id]);
			n->draw( false, true );
		}
	}

	// Draw target
	{
		Structure::Graph * g = path->gcorr->tg;

		glViewport(dx + trect.x(), dy + viewport[3] - trect.height() - trect.top(), trect.width(), trect.height());

		for(auto n : g->nodes) 
		{
			glColor3d(0,0,0);
			if(path->tcolors.contains(n->id)) glColorQt(path->tcolors[n->id]);
			n->draw( false, true );
		}
	}

	// Draw in between
	{ 
        Structure::Graph * g = path->scheduler->allGraphs[path->i];

		glViewport(dx + inbetween.x(), dy + viewport[3] - inbetween.height() - inbetween.top(), inbetween.width(), inbetween.height());

		for(auto n : g->nodes) 
		{
			glColor3d(0,0,0);

			QString sid = n->property["original_ID"].toString();
			if(path->scolors.contains(sid)) glColorQt(path->scolors[sid]);

			n->draw( false, true );
		}
	}

	glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);
	painter->endNativePainting();

	// Draw graph
	{
		QRectF graphRect = inbetween;
		graphRect.setHeight( inbetween.height() * 0.10 );

		QPainterPath graph;

		painter->save();
		painter->translate( graphRect.topLeft() );

		double sum = 0;

		for(int i = 0; i < path->errors.size(); i++)
		{
			double dx = (double(i) / (path->errors.size()-1));
			double dy = (sum += path->errors[i]) / path->weight;

			double x = dx * graphRect.width();
			double y = graphRect.height() - (dy * graphRect.height());

			if( i == 0 ) graph.moveTo(x,y);
			else graph.lineTo(x,y);
		}

		painter->setPen( QPen(QColor(0,128,0), 2) );
		painter->drawPath(graph);

		painter->restore();
	}
}
