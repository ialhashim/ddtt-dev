#pragma once

#include <QGraphicsObject>

struct DeformationPath;
typedef QMap<QString,QVariant> PropertyMap;

class DeformPathItem : public QGraphicsObject
{
	Q_OBJECT
public:
    DeformPathItem(DeformationPath * usedPath = NULL);
    PropertyMap property;

	DeformationPath * path;
	QRectF m_rect;

public:
	QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
};
