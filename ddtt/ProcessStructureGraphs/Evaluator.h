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
	explicit Evaluator(QString datasetPath, bool isSet = false, QWidget *parent = 0);
	explicit Evaluator(QString datasetPath, std::vector<std::vector<std::pair<QString, QString>>> &allMaps,	std::vector<std::vector<std::pair<QString, QString>>> &allMapsLabel);
	~Evaluator();

	void run();
	void compareWithGreedyOBB(std::vector<std::vector<std::pair<QString, QString>>> &allMaps,
		std::vector<std::vector<std::pair<QString, QString>>> &allMapsLabel);

    QString datasetPath;
	bool isSet;

private:
    Ui::Evaluator *ui;
};
