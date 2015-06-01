#pragma once

#include <QWidget>
#include <QThread>

namespace Ui {
class Evaluator;
}

class Evaluator : public QWidget
{
    Q_OBJECT

public:
    explicit Evaluator(QString datasetPath, bool isSet = false, bool optClustering = true, bool optGTMode = false, QWidget *parent = 0);
    ~Evaluator();

	void run();

    QString datasetPath;
	bool isSet;
	bool optClustering;
	bool optGTMode;

private:
    Ui::Evaluator *ui;
};
