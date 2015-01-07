#pragma once

#include <QString>
#include <QVector>
#include <QObject>
#include <QColor>
#include <QMessageBox>
#include <QPainter>

#define SAFE_DELETE( ptr ) \
	if(ptr) \
			{ \
		delete ptr; \
		ptr = NULL; \
			}
#define SAFE_DELETE_ARRAY( ptr ) \
	if(ptr) \
			{ \
		delete [] ptr; \
		ptr = NULL; \
			}

const int ColorNum = 16;

const int ColorSet[16][3] = {
		{ 130, 130, 240 }, { 255, 120, 120 }, { 46, 254, 100 }, { 250, 88, 172 },
		{ 250, 172, 88 }, { 129, 247, 216 }, { 200, 200, 50 }, { 226, 169, 143 }, { 8, 138, 41 },
		{ 1, 223, 215 }, { 11, 76, 95 }, { 190, 182, 90 },
		{ 245, 169, 242 }, { 75, 138, 8 }, { 247, 254, 46 }, { 88, 172, 250 }
};

static QColor GetColorFromSet(int i)
{
	i = i%ColorNum;
	return QColor(ColorSet[i][0], ColorSet[i][1], ColorSet[i][2], 255);
}

// order: TopRightCeil  TopRight BottomRight BottomRightCeil TopLeftCeil TopLeft BottomLeft BottomLeftCeil    
const int BoxFaceVertex[12][3] = {
	// triangular faces
		{ 0, 1, 2 }, { 2, 3, 0 },
		{ 6, 5, 4 }, { 4, 7, 6 },
		{ 1, 0, 4 }, { 4, 5, 1 },
		{ 1, 5, 6 }, { 6, 2, 1 },
		{ 2, 6, 7 }, { 7, 3, 2 },
		{ 0, 3, 7 }, { 7, 4, 0 }
};

const int BoxFaceNormalCoeff[12][2] = {
	// face normal orientation along axis
	// {axis_id (along which axis), orientation (1 means pointing to positive and -1 to negative)}
		{ 0, 1 }, { 0, 1 },
		{ 0, -1 }, { 0, -1 },
		{ 2, 1 }, { 2, 1 },
		{ 1, -1 }, { 1, -1 },
		{ 2, -1 }, { 2, -1 },
		{ 1, 1 }, { 1, 1 }
};

static void Simple_Message_Box(const QString &text)
{
	QMessageBox msg;
	msg.setText(text);
	msg.setButtonText(QMessageBox::Ok, "OK");
	msg.exec();
}

static void Simple_Message_Box(const char* text)
{
	QMessageBox msg;
	msg.setText(QString(text));
	msg.setButtonText(QMessageBox::Ok, "OK");
	msg.exec();
}

static void DrawText_QPixmap(QPixmap *p, const QString &text, const QPoint pos = QPoint(50, 50), const QColor c = QColor(0, 255, 255, 255))
{
	QPainter painter(p);
	QPen pen(c);
	painter.setPen(c);

	QFont font(QFont("Arial"));
	font.setPointSize(24);
	painter.setFont(font);
	painter.drawText(pos, text);
}

