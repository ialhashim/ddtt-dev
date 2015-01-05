#include "SymRelationViz-widget.h"
#include "ui_SymRelationViz-widget.h"

SymRelationVizWidget::SymRelationVizWidget(QWidget *parent) : QWidget(parent), ui(new Ui::SymRelationVizWidget)
{
    ui->setupUi(this);
}

SymRelationVizWidget::~SymRelationVizWidget()
{
    delete ui;
}
