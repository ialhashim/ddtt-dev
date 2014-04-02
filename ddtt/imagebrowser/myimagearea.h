#pragma once
#include <QStringList>
#include <QColor>
#include <QWidget>
#include <QImage>

extern QString curLabel;
extern QStringList AllLabels;
extern QVector<QColor> UniqueColors;
extern QVector<QString> labelNames;

static QString imageOps[] = { "none" };
enum ImageOperation{ NONE_OP };
extern ImageOperation curOp;

class MyImageArea : public QWidget
{
	Q_OBJECT
public:
    MyImageArea(QString filename) : filename(filename) { isDeleted = false; img = QImage(filename); }
	~MyImageArea();

	void draw();

	void mousePressEvent(QMouseEvent * event);
	void focusInEvent(QFocusEvent * event);

    QImage img;
    QString filename;
	bool isDeleted;

signals:
	void gotFocus(MyImageArea *);
};

extern MyImageArea * lastSelected;

