#include "mainwindow.h"
#include <QApplication>
#include <QFileDialog>
#include <QCommandLineParser>
#include <QProcess>
#include "BatchProcess.h"

enum CommandLineParseResult{
    CommandLineOk,
    CommandLineError,
    CommandLineVersionRequested,
    CommandLineHelpRequested
};

typedef QMap<QString, PropertyMap> DatasetMap;
DatasetMap shapesInDataset(QString datasetPath)
{
	DatasetMap dataset;

	QDir datasetDir(datasetPath);
	QStringList subdirs = datasetDir.entryList(QDir::Dirs | QDir::NoSymLinks | QDir::NoDotAndDotDot);

	for(auto subdir : subdirs)
	{
		// Special folders
		if (subdir == "corr") continue;

		QDir d(datasetPath + "/" + subdir);

		// Check if no graph is in this folder
		if (d.entryList(QStringList() << "*.xml", QDir::Files).isEmpty()) continue;

		auto xml_files = d.entryList(QStringList() << "*.xml", QDir::Files);

		dataset[subdir]["Name"] = subdir;
		dataset[subdir]["graphFile"] = d.absolutePath() + "/" + (xml_files.empty() ? "" : xml_files.front());
		dataset[subdir]["thumbFile"] = d.absolutePath() + "/" + d.entryList(QStringList() << "*.png", QDir::Files).join("");
		dataset[subdir]["objFile"] = d.absolutePath() + "/" + d.entryList(QStringList() << "*.obj", QDir::Files).join("");
	}

	return dataset;
}

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
			{ { "f", "folder" }, QString("Folder for a shape dataset."), QString("folder") },
			{ { "q", "quiet" }, QString("Skip visualization.") },
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

    // Process shape sets
    if(parser.isSet("folder"))
    {
		QElapsedTimer timer; timer.start();

        QString dataset = parser.value("folder");
		QDir d("");

		// 1) get sorted set of pairs to compare
		QVector< QPair<PropertyMap, PropertyMap> > shapePairs;
		auto folders = shapesInDataset(dataset);
		auto folderKeys = folders.keys();
		for (int i = 0; i < folderKeys.size(); i++)
			for (int j = i + 1; j < folderKeys.size(); j++)
				shapePairs << qMakePair(folders[folderKeys.at(i)], folders[folderKeys.at(j)]);
		int shapePairsCount = shapePairs.size();
		int curShapePair = 0;
		
		// 2) perform correspondence for shape pairs
		{
			// Remove previous log file
			d.remove("log.txt");

			// go over pairs
			for (auto shapePair : shapePairs)
			{
				QProcess p;
				QString source = shapePair.first["graphFile"].toString();
				QString target = shapePair.second["graphFile"].toString();

				auto sg = QFileInfo(source).baseName();
				auto tg = QFileInfo(target).baseName();

				std::cout << QString("Now: %1 / %2").arg(sg).arg(tg).leftJustified(35, ' ', true).toStdString();

				QStringList pargs;

				// input
				pargs << "-s" << source << "-t" << target;	

				// search options
				if (parser.isSet("o")) pargs << "-o";		
				if (parser.isSet("k")) pargs << "-k" << parser.value("k");
				if (parser.isSet("q")) pargs << "-q";

				p.start(a.applicationFilePath(), pargs);
				p.waitForFinished(-1);

				auto percent = (double(curShapePair++) / shapePairsCount) * 100.0;
				std::cout << QString("[%1 %] - %2 - ").arg((int)percent).arg(shapePairsCount-curShapePair).toStdString();
				int secPerPercent = (timer.elapsed() / 1000) / (percent + 1);
				int timeMinLeft = (secPerPercent * (100 - percent)) / 60;
				std::cout << QString("Time (%1s=%2m) ETA(%3m)\n").arg(timer.elapsed() / 1000).arg(timer.elapsed() / 60000).arg(timeMinLeft).toStdString();
			}
		}

		// 3) prepare results folder
		QString job_name = QString("job_%1").arg(QDateTime::currentDateTime().toString("dd_MM_yyyy_hh_mm_ss"));
		d.mkpath(job_name);

		// 4) collect results
		{
			// output sorted set of shapes
			{
				QFile file(d.absolutePath() + "/" + job_name + "/" + "shapes.txt");
				if (file.open(QIODevice::WriteOnly | QIODevice::Text)){
					QTextStream out(&file);
					for (int i = 0; i < folderKeys.size(); i++)
						out << QString("%1 %2\n").arg(i).arg(folderKeys.at(i));
				}
			}

			// read log and record the minimum costs and correspondence results
			{
				QFile ff(d.absolutePath() + "/" + "log.txt");
				ff.open(QFile::ReadOnly | QFile::Text);
				QTextStream in(&ff);
				auto datasetLogLines = in.readAll().split("\n", QString::SkipEmptyParts);

				// Hopefully everything went well
				assert(datasetLogLines.size() == shapePairs.size());
				if (datasetLogLines.size() != shapePairs.size()) return -1;

				// Record final results
				QJsonArray results;
				int idx = 0;
				for (int i = 0; i < folderKeys.size(); i++){
					for (int j = i + 1; j < folderKeys.size(); j++){
						// Read matching pair info
						auto log_line = datasetLogLines[idx].split(",", QString::SkipEmptyParts);
						QJsonObject matching;
						matching["source"] = log_line.at(0);
						matching["target"] = log_line.at(1);
						matching["cost"] = log_line.at(2);

						// Get correspondence data
						QString correspondenceFile = log_line.at(3);
						bool isSwapped = (log_line.at(4).toInt() % 2) == 0;
						QJsonArray correspondence;
						{
							// Open correspondence file
							QFile cf(correspondenceFile);
							cf.open(QFile::ReadOnly | QFile::Text);
							QTextStream cfin(&cf);
							auto corrLines = cfin.readAll().split("\n", QString::SkipEmptyParts);

							// Read correspondence file (swapping if needed)
							for (auto line : corrLines)
							{
								auto matched_pair = line.split(" ", QString::SkipEmptyParts);
								QString sid = matched_pair.at(0);
								QString tid = matched_pair.at(1);
								if (isSwapped) std::swap(sid, tid);

								QJsonArray part_pair;
								part_pair.push_back(sid);
								part_pair.push_back(tid);
								correspondence.push_back(part_pair);
							}
						}
						matching["correspondence"] = correspondence;

						// indexing
						matching["i"] = i;
						matching["j"] = j;

						// Record result
						results.push_back(matching);
						idx++;
					}
				}

				// Write all results in JSON format
				{
					QJsonDocument saveDoc(results);
					QFile saveFile(d.absolutePath() + "/" + job_name + "/" + "dataset.corr");
					saveFile.open(QIODevice::WriteOnly);
					saveFile.write(saveDoc.toJson());
				}
			}
		}

		return folders.size();
    }

    QString jobs_filename;
    if(parser.isSet("nogui") || parser.isSet("auto") || parser.isSet("sourceShape"))
    {
        if (parser.isSet("auto") || parser.isSet("sourceShape") || parser.isSet("targetShape"))
		{
			QString sourceShape = parser.value("sourceShape");
			QString targetShape = parser.value("targetShape");
			QVariantMap options;

            if(parser.isSet("g")) options["align"].setValue(true);
            if(parser.isSet("o")) options["roundtrip"].setValue(true);
            if(parser.isSet("k")) options["k"].setValue(parser.value("k").toInt());
			if(parser.isSet("q")) options["isQuietMode"].setValue(true);

            if(options["roundtrip"].toBool() || options["align"].toBool())
			{
				options["isManyTypesJobs"].setValue(true);
				options["isOutputMatching"].setValue(true);

				QTimer::singleShot(1, [sourceShape, targetShape, options]
				{
					int numJobs = 0;

					auto cur_options = options;

					QVector< QVector<QVariantMap> > reports;

					int numIter = 1;

					// Command line now only supports two tests.. GUI has more options
					if (cur_options["align"].toBool()) numIter = 2;

					for (int iter = 0; iter < numIter; iter++)
					{
						auto bp = QSharedPointer<BatchProcess>(new BatchProcess(sourceShape, targetShape, cur_options));

						bp->jobUID = numJobs++;
						bp->run();

						reports << bp->jobReports;

						if (cur_options["roundtrip"].toBool())
						{
							auto bp2 = QSharedPointer<BatchProcess>(new BatchProcess(targetShape, sourceShape, cur_options));

							bp2->jobUID = numJobs++;
							bp2->run();

							reports << bp2->jobReports;
						}

						cur_options["isFlip"].setValue(true);
					}

					// Look at reports
					double minEnergy = 1.0;
					int totalTime = 0;
					QVariantMap minJob;
					for (auto & reportVec : reports)
					{
						for (auto & report : reportVec)
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
							out << sourceShape << "," << targetShape << "," << minJob["min_cost"].toDouble() << "," << minJob["match_file"].toString() << "," << minJob["is_swapped"].toInt() << "\n";
						}
					}
					
				}
				);
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
