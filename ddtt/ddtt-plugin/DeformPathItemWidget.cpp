#pragma warning(disable:4267)

#include <QWidget>
#include <QVBoxLayout>
#include <QSlider>
#include <QLabel>
#include <QLineEdit>
#include "DeformPathItemWidget.h"
#include "ShapeCorresponder.h"

DeformPathItemWidget::DeformPathItemWidget(DeformationPath * usedPath) : path(usedPath)
{
	// Create a widget with a slider and a progress bar
	QWidget * w = new QWidget;

	int width = 800 - (600 * 0.5);
	int height = 75;
	w->setMinimumSize(width, height);
	w->setMaximumSize(width, height);

	//w->setStyleSheet("*{ border: 1px solid red; }");

	// Set widget
	setWidget(w);

	connect(this, SIGNAL(widgetCreated()), SLOT(init()));
	emit( widgetCreated() );
}

void DeformPathItemWidget::init()
{
	QVBoxLayout * layout = new QVBoxLayout;

    slider = new QSlider(Qt::Horizontal);
	slider->setTickPosition(QSlider::TicksBothSides);
    slider->setMinimum(0);
    slider->setMaximum(path->scheduler->allGraphs.size()-1);

    QString style = "QSlider::groove:horizontal {background: blue;height: 4px;}";
    style += "QSlider::handle:horizontal {background: #ff0000;width: 20px;margin: -16px 0px -16px 0px;}";
    slider->setStyleSheet( style );

	layout->addWidget(slider);
    layout->addWidget(label = new QLabel("label"));
	//QLineEdit *numberEdit = new QLineEdit;layout->addWidget(numberEdit);

	widget()->setLayout(layout);

    // Connections
    connect(slider,SIGNAL(valueChanged(int)),this,SLOT(sliderValueChanged(int)));
}

void DeformPathItemWidget::sliderValueChanged(int val)
{
    label->setText( QString("In between: %1").arg(val) );
    path->i = val;
}
