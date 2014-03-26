#include "DeformScene.h"
#include "ShapeCorresponder.h"
#include <QGraphicsView>
#include <QHBoxLayout>
#include <QDesktopWidget>
#include <QApplication>

#include "DeformPathItem.h"
#include "DeformPathItemWidget.h"

DeformScene::DeformScene()
{
	// Anti-aliasing when using QGLWidget or subclasses
	{
		QGLFormat glf = QGLFormat::defaultFormat();
		glf.setSamples(8);
		QGLFormat::setDefaultFormat(glf);
	}

	QGLWidget * viewport = new QGLWidget;
	viewport->makeCurrent();

	QGraphicsView * gview = new QGraphicsView;

	gview->setViewport( viewport );
	gview->setViewportUpdateMode( QGraphicsView::FullViewportUpdate );
	gview->setScene( this );
	gview->setAlignment(Qt::AlignLeft | Qt::AlignTop);
	gview->setRenderHint(QPainter::HighQualityAntialiasing, true);
	gview->setRenderHint(QPainter::SmoothPixmapTransform, true);
	gview->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
	gview->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);

	// Show the widget
	QWidget * mainwidget = new QWidget;
	QHBoxLayout * layout = new QHBoxLayout;
	layout->addWidget(viewport);
	mainwidget->setLayout(layout);

	gview->resize(800,600);

	// Center to screen
	QDesktopWidget* m = QApplication::desktop();
	QRect desk_rect = m->screenGeometry(m->screenNumber(QCursor::pos()));
	int desk_x = desk_rect.width();
	int desk_y = desk_rect.height();
	int x = this->width();
	int y = this->height();

	mainwidget->move(desk_x / 2 - x / 2 + desk_rect.left(), desk_y / 2 - y / 2 + desk_rect.top());
	mainwidget->show();
	mainwidget->raise();
}

void DeformScene::drawBackground(QPainter *painter, const QRectF &rect)
{
    QGraphicsScene::drawBackground(painter,rect);
	painter->drawRect(this->sceneRect());
}

void DeformScene::drawForeground(QPainter *painter, const QRectF &rect)
{
    QGraphicsScene::drawForeground(painter, rect);
}

void DeformScene::addDeformationPath(DeformationPath * path)
{
	// Create a DeformationItem and its widget
	DeformPathItem * ditem = new DeformPathItem(path);
	DeformPathItemWidget * dwidget = new DeformPathItemWidget(path);
	
	//addItem( ditem );
	QGraphicsItemGroup * group = new QGraphicsItemGroup;
	group->addToGroup(ditem);
	group->addToGroup(dwidget);

	group->setPos(0, path->idx * 600);

	group->setAcceptHoverEvents(true);
	group->setFlags( QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable );

	addItem(group);
}
