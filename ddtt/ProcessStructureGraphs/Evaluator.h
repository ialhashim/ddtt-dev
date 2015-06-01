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
<<<<<<< HEAD
	explicit Evaluator(QString datasetPath, bool isSet = false, QWidget *parent = 0);
	explicit Evaluator(QString datasetPath, std::vector<std::vector<std::pair<QString, QString>>> &allMaps,	std::vector<std::vector<std::pair<QString, QString>>> &allMapsLabel);
	~Evaluator();
=======
    explicit Evaluator(QString datasetPath, bool isSet = false, bool optClustering = true, bool optGTMode = false, QWidget *parent = 0);
    ~Evaluator();
>>>>>>> 600293f0f3486770c72e24bb81fc9c31c72b23f4

	void run();
	void compareWithGreedyOBB(std::vector<std::vector<std::pair<QString, QString>>> &allMaps,
		std::vector<std::vector<std::pair<QString, QString>>> &allMapsLabel);

    QString datasetPath;
	bool isSet;
	bool optClustering;
	bool optGTMode;

private:
    Ui::Evaluator *ui;
};
