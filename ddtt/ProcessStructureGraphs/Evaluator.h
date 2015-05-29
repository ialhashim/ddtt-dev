#pragma once

#include <QWidget>

namespace Ui {
class Evaluator;
}

class Evaluator : public QWidget
{
    Q_OBJECT

public:
    explicit Evaluator(QString datasetPath, bool isSet = false, QWidget *parent = 0);
    ~Evaluator();

    QString datasetPath;
	bool isSet;

private:
    Ui::Evaluator *ui;
};
