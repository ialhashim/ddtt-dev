#include "mainwindow.h"
#include <QApplication>
#include <QFileDialog>
#include <QCommandLineParser>
#include "BatchProcess.h"

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
    parser.addOptions({
			{ "nogui", QString("Using command line arguments, do not show GUI.") },
			{ { "p", "path" }, QString("Path for job files or shape files."), QString("path") },
			{ { "a", "auto" }, QString("Automatically try to find initial correspondence. Not used for job files.") },
			{ { "s", "sourceShape" }, QString("Path for job files or shape files."), QString("source") },
			{ { "t", "targetShape" }, QString("Path for job files or shape files."), QString("target") },
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
    if(parser.isSet("nogui") || parser.isSet("auto"))
    {
		if (parser.isSet("auto") && parser.isSet("sourceShape") && parser.isSet("targetShape"))
		{
			QString sourceShape = parser.value("sourceShape");
			QString targetShape = parser.value("targetShape");

			BatchProcess * bp = new BatchProcess(sourceShape, targetShape);
			QObject::connect(bp, SIGNAL(allJobsFinished()), &w, SLOT(close()));
			QObject::connect(bp, SIGNAL(finished()), bp, SLOT(deleteLater()));
			bp->start();

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
