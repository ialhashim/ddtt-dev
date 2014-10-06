#include "deform-widget.h"
#include "ui_deform-widget.h"

DeformWidget::DeformWidget(QWidget *parent) : QWidget(parent), ui(new Ui::DeformWidget)
{
    ui->setupUi(this);
}

DeformWidget::~DeformWidget()
{
    delete ui;
}
