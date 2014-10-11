#include "experiment-widget.h"
#include "ui_experiment-widget.h"

ExperimentWidget::ExperimentWidget(QWidget *parent) : QWidget(parent), ui(new Ui::ExperimentWidget)
{
    ui->setupUi(this);
}

ExperimentWidget::~ExperimentWidget()
{
    delete ui;
}
