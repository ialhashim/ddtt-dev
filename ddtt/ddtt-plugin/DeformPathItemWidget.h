#pragma once

#include <QGraphicsProxyWidget>
#include <QSlider>
#include <QLabel>

class DeformationPath;

class DeformPathItemWidget : public QGraphicsProxyWidget
{
	Q_OBJECT
public:
	DeformPathItemWidget(DeformationPath * usedPath = NULL);

	DeformationPath * path;

    QSlider * slider;
    QLabel * label;

public slots:
	void init();
    void sliderValueChanged(int val);

signals:
	void widgetCreated();
};
