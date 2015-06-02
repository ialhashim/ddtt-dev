#pragma once

#include <QWidget>
#include <QThread>
#include <QVariant>

namespace Ui {
class Evaluator;
}

class Evaluator : public QWidget
{
    Q_OBJECT

public:
    explicit Evaluator(QString datasetPath, bool isSet = false, bool optClustering = true, 
		bool optGTMode = false, QVariantMap otherOptions = QVariantMap(), QWidget *parent = 0);
    ~Evaluator();

	Evaluator::Evaluator(QString datasetPath, std::vector<std::vector<std::pair<QString, QString>>> &allMaps, 
        std::vector<std::vector<std::pair<QString, QString>>> &allMapsLabel, bool isSet = false) : datasetPath(datasetPath){
        compareWithGreedyOBB(allMaps, allMapsLabel, isSet);
	}

	void run();

    void compareWithGreedyOBB(std::vector<std::vector<std::pair<QString, QString>>> &allMaps,
                              std::vector<std::vector<std::pair<QString, QString>>> &allMapsLabel, bool isSet = false);

    QString datasetPath;
	bool isSet;
	bool optClustering;
	bool optGTMode;

	QVariantMap otherOptions;

private:
    Ui::Evaluator *ui;
};
