#include "mainwindow.h"
#include <QApplication>
#include <QFileDialog>
#include <QCommandLineParser>
#include "BatchProcess.h"

QVariantMap options;

enum CommandLineParseResult{
    CommandLineOk,
    CommandLineError,
    CommandLineVersionRequested,
    CommandLineHelpRequested
};

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    //w.show();

    QCommandLineParser parser;
    parser.addHelpOption();
    parser.addOptions({
			{ "nogui", QString("Using command line arguments, do not show GUI.") },
            { { "p", "path" }, QString("Path for job files or shape files."), QString("path") },
            { { "g", "align" }, QString("Input is not aligned. Find lowest cost alignment.") },
            { { "o", "roundtrip" }, QString("Compute least cost from source to target, and target to source.") },
            { { "s", "sourceShape" }, QString("Path for source shape file."), QString("source") },
            { { "t", "targetShape" }, QString("Path for target shape file."), QString("target") },
            { { "k", "k" }, QString("(k) parameter for DP search."), QString("k") },
            { { "a", "auto" }, QString("Automatically try to find initial correspondence. Not used for job files.") },
			{ { "j", "job" }, QString("Job file to load."), QString("job") },
	});

    if (!parser.parse(QCoreApplication::arguments())) {
        QString errorText = parser.errorText();
        std::cout << qPrintable(errorText);
        return CommandLineError;
    }
	else
		parser.process(a);

	QString path = parser.value("p");
	QDir::setCurrent(path);

    QString jobs_filename;
    if(parser.isSet("nogui") || parser.isSet("auto") || parser.isSet("sourceShape"))
    {
        if (parser.isSet("auto") || parser.isSet("sourceShape") || parser.isSet("targetShape"))
		{
			QString sourceShape = parser.value("sourceShape");
			QString targetShape = parser.value("targetShape");

            if(parser.isSet("g")) options["align"].setValue(true);
            if(parser.isSet("o")) options["roundtrip"].setValue(true);
            if(parser.isSet("k")) options["k"].setValue(parser.value("k").toInt());

            if(options["roundtrip"].toBool() || options["align"].toBool())
            {
				QTimer::singleShot(1000, [sourceShape, targetShape] {
					int numJobs = 0;

					QVector< QVector<QVariantMap> > reports;

					options["isManyTypesJobs"].setValue(QString("isManyTypesJobs"));

					int numIter = 1;

					// Command line now only supports two tests.. GUI has more options
					if (options["align"].toBool()) numIter = 2;

					for (int iter = 0; iter < numIter; iter++)
					{
						auto bp = new BatchProcess(sourceShape, targetShape, options);
						//QObject::connect(bp, SIGNAL(allJobsFinished()), &w, SLOT(close()));
						QObject::connect(bp, SIGNAL(finished()), bp, SLOT(deleteLater()));

						bp->jobUID = numJobs++;
						bp->start();
						while (bp->isRunning())
						{
							qApp->processEvents();
						}
                        qApp->processEvents();

						reports << bp->jobReports;

						if (options["roundtrip"].toBool())
						{
                            auto bp2 = new BatchProcess(targetShape, sourceShape, options);

							bp2->jobUID = numJobs++;
							bp2->start();
                            while (bp2->isRunning())
                            {
                                qApp->processEvents();
                            }
                            qApp->processEvents();

							reports << bp2->jobReports;
						}

						options["isFlip"].setValue(true);
					}

					// Look at reports
					double minEnergy = 1.0;
					int totalTime = 0;
					QVariantMap minJob;
					for (auto reportVec : reports)
					{
						for (auto report : reportVec)
						{
							totalTime += report["search_time"].toInt();

							double c = report["min_cost"].toDouble();
							if (c < minEnergy){
								minEnergy = c;
								minJob = report;
							}
						}
					}

					std::cout << "\nJobs computed: " << numJobs << "\n";
					std::cout << minEnergy << " - " << qPrintable(minJob["img_file"].toString());

					// Aggregate results
					{
						QFile file("log.txt");
						if (file.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Append))
						{
							QTextStream out(&file);
							out << sourceShape << "," << targetShape << "," << minJob["min_cost"].toDouble() << "\n";
						}
					}
					
				});
            }
            else
            {
                auto bp = new BatchProcess(sourceShape, targetShape, options);
                QObject::connect(bp, SIGNAL(allJobsFinished()), &w, SLOT(close()));
                QObject::connect(bp, SIGNAL(finished()), bp, SLOT(deleteLater()));
                bp->start();
            }

			return a.exec();
		}
		else
		{
			jobs_filename = parser.value("job");
			if (jobs_filename.isEmpty()){
				std::cout << "Please provide a job file.";
				return CommandLineError;
			}
		}

		std::cout << "Not enough arguments.";
		return CommandLineError;
    }
    else
    {
		jobs_filename = QFileDialog::getOpenFileName(&w, "Load Jobs", "", "Jobs File (*.json)");
    }

    QTimer::singleShot(0, [&] {
        BatchProcess * bp = new BatchProcess(jobs_filename);
        QObject::connect(bp, SIGNAL(allJobsFinished()), &w, SLOT(close()));
		QObject::connect(bp, SIGNAL(finished()), bp, SLOT(deleteLater()));
        bp->start();
    });

    return a.exec();
}
