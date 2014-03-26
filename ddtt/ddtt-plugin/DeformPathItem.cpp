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

void DeformPathItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
	int viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	int width = m_rect.width();
	int height = m_rect.height();
	int smallWidth = height * 0.5;
	int tinyWidth = height * (1.0 / 8.0);

	QGraphicsView * view = scene()->views().front();
	QPointF sceneP = mapToScene( QPointF(0,0) );
	QPoint viewP = view->mapFromScene(sceneP);

	double dy = -viewP.y();
	double dx = viewP.x();

	if( dy > viewport[3] || dy < -viewport[3] ) return;

	painter->drawText(QPoint(10,200), QString("IDX %3  = %1, %2, dy = %4").arg(viewP.x()).arg(viewP.y()).arg(path->idx).arg(dy));

	QRectF srect(width - smallWidth ,0,smallWidth,smallWidth);
	painter->setPen( QPen(Qt::green, 2) );
	painter->drawRect( srect );

	QRectF trect(width - smallWidth,smallWidth,smallWidth,smallWidth);
	painter->setPen( QPen(Qt::blue, 2) );
	painter->drawRect( trect );

	QRectF inbetween(0, tinyWidth, width - smallWidth, height - tinyWidth);
	painter->setPen( QPen(Qt::black, 2) );
	painter->drawRect( inbetween );

	painter->beginNativePainting();

	setCamera2(srect, path->gcorr->sg);

	// Draw source
	{
		Structure::Graph * g = path->gcorr->sg;

		glViewport(dx + srect.x(), dy + viewport[3] - srect.height() - srect.top(), srect.width(), srect.height());
		for(auto n : g->nodes) n->draw( true );
	}

	// Draw target
	{
		Structure::Graph * g = path->gcorr->tg;

		glViewport(dx + trect.x(), dy + viewport[3] - trect.height() - trect.top(), trect.width(), trect.height());
		for(auto n : g->nodes) n->draw( true );
	}

	// Draw in between
	{
		Structure::Graph * g = path->scheduler->allGraphs.front();

		glViewport(dx + inbetween.x(), dy + viewport[3] - inbetween.height() - inbetween.top(), inbetween.width(), inbetween.height());
		for(auto n : g->nodes) n->draw( true );
	}

	glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);

	painter->endNativePainting();
}
