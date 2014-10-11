#ifndef EXPERIMENTWIDGET_H
#define EXPERIMENTWIDGET_H

#include <QWidget>

namespace Ui {
class ExperimentWidget;
}

class ExperimentWidget : public QWidget
{
    Q_OBJECT

public:
    explicit ExperimentWidget(QWidget *parent = 0);
    ~ExperimentWidget();

private:
    Ui::ExperimentWidget *ui;
};

#endif // EXPERIMENTWIDGET_H
