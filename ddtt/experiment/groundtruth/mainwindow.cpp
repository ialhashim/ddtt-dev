#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDir>
#include <QFile>
#include <QDebug>
#include <QMessageBox>

struct Record{
    QString shapeA, shapeB;
    QString p,q;
    QString user;
    Record(){}
    Record(QString shapeA, QString shapeB, QString p, QString q, QString user):shapeA(shapeA), shapeB(shapeB),p(p),q(q),user(user){}

    bool operator < (const Record& rhs) const {
        return (p+"-"+q+"-"+shapeA+"-"+shapeB) < (rhs.p+"-"+rhs.q+"-"+rhs.shapeA+"-"+rhs.shapeB);
    }

    bool operator ==(const Record& rhs) const {
        return p == rhs.p && q == rhs.q && shapeA == rhs.shapeA && shapeB == rhs.shapeB && user == rhs.user;
    }
};

void computePerformance(double threshold = 0.8)
{
	QVector<Record> records;

	QString datasetPath = "C:/Users/Ibraheem/Dropbox/DDTT-dataset/evaluation/userstudy";

	QDir datasetDir(datasetPath);
	QStringList subdirs = datasetDir.entryList(QDir::Dirs | QDir::NoSymLinks | QDir::NoDotAndDotDot);

	QStringList users;

	for (QString subdir : subdirs)
	{
		QDir d(datasetPath + "/" + subdir);

		auto matchFiles = d.entryList(QStringList() << "*.match", QDir::Files);
		if (matchFiles.isEmpty()) continue;

		users << subdir;

		for (QString matchFile : matchFiles)
		{
			// Get matches from file
			QFile ff(d.absolutePath() + "/" + matchFile);
			ff.open(QFile::ReadOnly | QFile::Text);
			QTextStream in(&ff);
			auto matches = in.readAll().split("\n", QString::SkipEmptyParts);

			// Shape names
			QStringList shapes = matchFile.replace(".match", "").split("-", QString::SkipEmptyParts);

			bool isSwap = shapes.front() > shapes.back();

			// Add to vector
			QVector< QPair<QString, QString> > pairMatches;
			for (QString matchLine : matches)
			{
				QStringList match = matchLine.split(" ", QString::SkipEmptyParts);

				QString p = match.front();
				QString q = match.back();

				if (isSwap) std::swap(p, q);

				pairMatches << qMakePair(p, q);
			}

			//qDebug() << pairMatches;

			// Add record
			if (isSwap) std::reverse(shapes.begin(), shapes.end());
			for (auto a : pairMatches) records << Record(shapes.front(), shapes.back(), a.first, a.second, subdir);
		}
	}

	// Debug records:
	qSort(records);
	for (Record & record : records)
	{
		//qDebug() << record.p << " " << record.q << " " << record.shapeA << " " << record.shapeB << " " << record.user;
	}

	//qDebug() << "Num records =" << records.size();

	// ShapePairs/part_p/part_q/users
	QMap<QString, QMap<QString, QMap<QString, QStringList> > > part_matches;

	// Extract set of good pair-wise part correspondences
	for (Record r : records)
	{
		part_matches[r.shapeA + "|" + r.shapeB][r.p][r.q] << r.user;
	}

	// Passed threshold
	int num_users = users.size();
	
	QVector<Record> C;

	for (QString shapePair : part_matches.keys())
	{
		auto cur_shape_pairs = shapePair.split("|");

		for (auto p : part_matches[shapePair].keys())
		{
			auto part_matching = part_matches[shapePair][p];

			for (auto q : part_matching.keys())
			{
				auto match = part_matching[q];

				int num_assigned = match.size();

				double v = double(num_assigned) / num_users;

				if (v >= threshold)
				{
					C << Record(cur_shape_pairs.front(), cur_shape_pairs.back(), p, q, "from_users");
				}
			}
		}
	}

	//qDebug() << "Num records passed =" << C.size();

	// Ground truth
	// ShapePair/part_p/part_q
	QMap <QString, QMap<QString, QString> > ground_truth;
	for (Record r : C)
	{
		ground_truth[r.shapeA + "|" + r.shapeB][r.p] = r.q;
	}

	// Test set:
	QString testPath = "C:/Users/Ibraheem/Dropbox/DDTT-dataset/evaluation/userstudy/geotopo";
	QDir testDir(testPath);

	int wins = 0;
	int loss = 0;

	auto matchFiles = testDir.entryList(QStringList() << "*.match", QDir::Files);
	for (auto matches_filename : matchFiles)
	{
		// Get matches from file
		QFile ff(testDir.absolutePath() + "/" + matches_filename);
		ff.open(QFile::ReadOnly | QFile::Text);
		QTextStream in(&ff);
		auto matches = in.readAll().split("\n", QString::SkipEmptyParts);

		QStringList shapes = matches_filename.replace(".match", "").split("-", QString::SkipEmptyParts);
		bool isSwap = shapes.front() > shapes.back();

		// Shape names sorted
		if (isSwap) std::reverse(shapes.begin(), shapes.end());
		QString shapeA = shapes.front();
		QString shapeB = shapes.back();

		// Get pair matches
		for (QString matchLine : matches)
		{
			QStringList match = matchLine.split(" ", QString::SkipEmptyParts);

			QString p = match.front();
			QString q = match.back();

			if (isSwap) std::swap(p, q);

			if (ground_truth[shapeA + "|" + shapeB].contains(p))
			{
				QString truth = ground_truth[shapeA + "|" + shapeB][p];

				if (truth == q)
					wins++;
				else
					loss++;
			}
		}
	}

	double performance = (double(wins) / (wins + loss));

	QString result = QString("Performance = %1, agreement %6%, selected = %4, wins = %2, loss = %3, all pairs = %5")
		.arg(performance).arg(wins).arg(loss).arg(C.size()).arg(records.size()).arg(threshold * 100);

	QMessageBox::information(0, "Performance", result);

	qDebug() << result;
}

MainWindow::MainWindow(QWidget *parent) :
QMainWindow(parent),
ui(new Ui::MainWindow)
{
	ui->setupUi(this);

	computePerformance(0.25);
	computePerformance(0.50);
	computePerformance(0.75);
	computePerformance(0.80);
	computePerformance(0.90);
	computePerformance(1.00);
}

MainWindow::~MainWindow()
{
    delete ui;
}
