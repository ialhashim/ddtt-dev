#include <QWidget>
#include <QVBoxLayout>
#include <QSlider>
#include <QLabel>
#include "DeformPathItemWidget.h"

DeformPathItemWidget::DeformPathItemWidget(DeformationPath * usedPath) : path(usedPath)
{
	// Create a widget with a slider and a progress bar
	QWidget * w = new QWidget;
	w->setMinimumSize(800 - (600 * 0.5), 75);
	w->setStyleSheet("*{ border: 1px solid red; }");

	QVBoxLayout * layout = new QVBoxLayout;

	QSlider * slider = new QSlider(Qt::Horizontal);
	slider->setMinimum(0);
	slider->setMaximum(100);
	slider->setValue(0);
	slider->setTickPosition(QSlider::TicksBothSides);

	layout->addWidget(slider);
	layout->addWidget(new QLabel("label"));

	w->setLayout(layout);

	// Connect signals

	// Set widget
	setWidget(w);
}
