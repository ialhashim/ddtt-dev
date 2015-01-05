#pragma once

#include <QWidget>

namespace Ui {
class SymRelationVizWidget;
}

class SymRelationVizWidget : public QWidget
{
    Q_OBJECT

public:
    explicit SymRelationVizWidget(QWidget *parent = 0);
    ~SymRelationVizWidget();

//private:
    Ui::SymRelationVizWidget *ui;
};
