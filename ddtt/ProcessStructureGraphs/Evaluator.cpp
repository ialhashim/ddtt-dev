#include "Evaluator.h"
#include "ui_Evaluator.h"

#include "globals.h"
#include "StructureGraph.h"

#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QFileDialog>

// QString exeCorresponder = "C:/Development/ddtt/ddtt/experiment/build-standalone-Qt_5_4-Release/release/geotopCorrespond.exe";
QString exeCorresponder = "C:/Users/Yixin/Desktop/StarlabPackageRelease/geotopCorrespond.exe";
//QString experiment = "C:/Temp/___Data/two_chairs";

struct MultiStrings{
	MultiStrings(QVector < QPair<QString, QString> > pairings = QVector < QPair<QString, QString> >()){
		for (auto p : pairings){
			groups[p.first] << p.first << p.second;
			groups[p.second] << p.first << p.second;
		}
	}
	QMap < QString, QSet<QString> > groups;
	QString representative(QString label){
		auto sorted = groups[label].toList();
		qSort(sorted);
		return sorted.front();
	}
};

struct MatchingRecord{
	QString S, T, sid, tid;
	MatchingRecord(QString sid="sid", QString tid="tid", QString label_s = "label_s",
		QString label_t = "label_t") : sid(sid), tid(tid), S(label_s), T(label_t) {
		if (this->tid == "") this->tid = sid;
		if (this->sid == "") this->sid = tid;
	}
	enum MyEnum{CORRECT_REC,INVALID_Rec,OTHER_REC} state;
};
typedef QVector < MatchingRecord > MatchingRecords;

struct LabelOracle{
	LabelOracle(){}
	QMap < QString, QSet<QString> > mapping;
	void push( QString first, QString second )
	{
		// An item can be itself
		mapping[first] << first << second;

		// Sibling can be item
		mapping[second] << second;
	}

	void build(){
		QVector < QPair<QString, QString> > combinations;
		for (auto l : mapping.keys()) for (auto k : mapping[l]) combinations << qMakePair(l, k);
		gt.possible = MultiStrings(combinations);
	}

	struct GroundTruth{
		QMap < QString, int > truth;
		MultiStrings possible;

		struct PrecisionRecall{
			double precision, recall, G, M, R;
			PrecisionRecall(double p = 0, double r = 0, double G = 0, double M = 0, double R = 0) 
				: precision(p), recall(r), G(G), M(M), R(R)
			{}
		};

		PrecisionRecall compute( MatchingRecords records ){
			auto correct_matches = truth;

			// Count of ground truth
			int G = 0;
			for (auto count : truth) G += count;

			// Count of correct matches
			int R = 0;
			for (auto record : records)
			{
				bool isExactMatch = false, isAcceptableMatch = false;

				// in case of broken data files
				if (record.T.trimmed() == "") record.T = record.S;
				if (record.S.trimmed() == "") record.S = record.T;

				auto coarse_s = possible.representative(record.S);
				auto coarse_t = possible.representative(record.T);

				isExactMatch = (record.S == record.T);

				// Only if one side is coarse do we go up a level
				bool isSourceCoarse = possible.representative(record.S) == record.S;
				bool isTargetCoarse = possible.representative(record.T) == record.T;
				if (isSourceCoarse || isTargetCoarse){
					if (!isExactMatch)
						isAcceptableMatch = (coarse_s == record.T || coarse_t == record.S);
				}

				if ( isExactMatch || isAcceptableMatch )
				{
					R++;
				}
			}

			// Count of returned matches
			int M = records.size();

			assert(M);
			assert(G);

			double p = double(R) / M;
			double r = double(R) / G;

			return PrecisionRecall(p, r, G, M, R);
		}
	};

	GroundTruth gt;

	void makeGroundTruth(QStringList source, QStringList target)
	{
		/// Remove labels with no equivalent whatsoever:
		auto removeIrrelevant = [](QStringList me, QStringList other, MultiStrings possible){
			QStringList result = me;

			for (auto label : me)
			{
				bool isFound = false;

				for (auto t : other){
					isFound = possible.representative(t) == possible.representative(label);
					if (isFound) break;
				}

				if (isFound) continue;

				result.removeAll(label);
			}
			
			return result;
		};

		// (b) find the only labels that should be considered
		source = removeIrrelevant(source, target, gt.possible);
		target = removeIrrelevant(target, source, gt.possible);

		// (c) each label will have the maximum number of appearance between two graphs
		QMap<QString, int> source_labels_counter, target_labels_counter;
		
		for (auto l : source) source_labels_counter[gt.possible.representative(l)]++;
		for (auto l : target) target_labels_counter[gt.possible.representative(l)]++;

		auto all_labels = source_labels_counter.keys().toSet() + target_labels_counter.keys().toSet();

		for (auto l : all_labels) gt.truth[l] = std::max(source_labels_counter[l], target_labels_counter[l]);
	}

	QVector<GroundTruth::PrecisionRecall> pr_results;
};

Evaluator::Evaluator(QString datasetPath, bool isSet, QWidget *parent) : QWidget(parent), ui(new Ui::Evaluator), datasetPath(datasetPath), isSet(isSet)
{
    ui->setupUi(this);

	run();
}
Evaluator::Evaluator(QString datasetPath, std::vector<std::vector<std::pair<QString, QString>>> &allMaps, std::vector<std::vector<std::pair<QString, QString>>> &allMapsLabel) : datasetPath(datasetPath)
{
	compareWithGreedyOBB(allMaps,allMapsLabel);
}

void Evaluator::run()
{
	//if (true) this->datasetPath = datasetPath = experiment;

	QDir dir(datasetPath);

	QString outputPath = datasetPath;
	QString resultsFile = outputPath + "/" + dir.dirName() + "_corr.json";

	QElapsedTimer ours_timer; ours_timer.start();

	if (!QFileInfo(exeCorresponder).exists())
	{
		exeCorresponder = QFileDialog::getOpenFileName(0, "geoCorresponder", "", "*.exe");
	}

	QString cmd = QString("%1 -o -q -k 2 -f %2 -z %3").arg(exeCorresponder).arg(datasetPath).arg(datasetPath);

	// Check first for cached results
	if (!QFileInfo(resultsFile).exists())
		system(qPrintable(cmd));

	auto all_pair_wise_time = ours_timer.elapsed();

	// Now process results
	if (!isSet)
	{
		/// Looking at pair-wise comparisons:

		// Open labels JSON file
		QMap<QString, QStringList> lables;
		LabelOracle oracle;
		{
			auto labelsFilename = datasetPath + "/labels.json";
			QFile file;
			file.setFileName(labelsFilename);
			if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;
			QJsonDocument jdoc = QJsonDocument::fromJson(file.readAll());
			auto json = jdoc.object();

			// Get lables
			auto labelsArray = json["labels"].toArray();
			for (auto l : labelsArray)
			{
				auto label = l.toObject();
				auto parent = label["parent"].toString();

				lables[parent] << label["title"].toString();
			}

			// Get cross-labels
			auto crossLabelsArray = json["cross-labels"].toArray();

			// regular nodes
			for (auto k : lables.keys()) for (auto l : lables[k]) oracle.push(l, l);

			// cross-labeled
			for (auto l : crossLabelsArray)
			{
				auto crosslabel = l.toObject();
				oracle.push(crosslabel["first"].toString(), crosslabel["second"].toString());
			}

			oracle.build();
		}

		// Open results JSON file
		QFile file;
		file.setFileName(resultsFile);
		if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;
		QJsonDocument jdoc = QJsonDocument::fromJson(file.readAll());
		auto corrArray = jdoc.toVariant().value<QVariantList>();

		QElapsedTimer stats_timer; stats_timer.start();

		for (auto c : corrArray)
		{
			auto obj = c.toMap();
			auto i = obj["i"].toInt();
			auto j = obj["j"].toInt();
			auto cost = obj["cost"].toDouble();
			auto corr = obj["correspondence"].value<QVariantList>();

			// Program might have crashed
			if (corr.isEmpty()) continue;

			if (obj["source"].toString() == obj["target"].toString()) continue;

			// Load graphs
			Structure::Graph source(obj["source"].toString()), target(obj["target"].toString());

			// Collect all labels from both shapes
			QStringList sourceLabels, targetLabels;
			for (auto n : source.nodes) if (n->meta.contains("label")) sourceLabels << n->meta["label"].toString();
			for (auto n : target.nodes) if (n->meta.contains("label")) targetLabels << n->meta["label"].toString();

			// Build expected ground truth
			oracle.makeGroundTruth(sourceLabels, targetLabels);

			MatchingRecords M;

			for (auto match : corr)
			{
				auto matching = match.value<QVariantList>();
				auto sid = matching.front().toString();
				auto tid = matching.back().toString();

				if (sid.isEmpty() || tid.isEmpty()) continue;

				auto sn = source.getNode(sid);
				auto tn = target.getNode(tid);

				M << MatchingRecord(sn->id, tn->id, sn->meta["label"].toString(), tn->meta["label"].toString());
			}

			oracle.pr_results << oracle.gt.compute(M);
		}

		// General stats
		int numG = 0, numM = 0, numR = 0;

		// Sum results
		double P = 0, R = 0;
		for (auto pr : oracle.pr_results)
		{
			P += pr.precision;
			R += pr.recall;

			numG += pr.G;
			numM += pr.M;
			numR += pr.R;
		}

		// Get average
		int N = oracle.pr_results.size();
		P /= N;
		R /= N;

		QString report = QString("[%5] Avg. P = %1, R = %2, Pair-wise time (%3 ms) - post (%4 ms)").arg(P).arg(R)
			.arg(all_pair_wise_time).arg(stats_timer.elapsed()).arg(dir.dirName());

		report += QString("\nG_count %1 / M_count %2 / R_count %3").arg(numG).arg(numM).arg(numR);

		// Save log of P/R measures
		{
			QFile logfile(datasetPath + "/log.txt");
			logfile.open(QIODevice::WriteOnly | QIODevice::Text);
			QTextStream out(&logfile);
			out << report + "\n";
		}

		debugBox(report);
	}
	else
	{

	}
}
void Evaluator::compareWithGreedyOBB( std::vector<std::vector<std::pair<QString, QString>>> &allMaps, 
	std::vector<std::vector<std::pair<QString, QString>>> &allMapsLabel)
{
	QDir dir(datasetPath);

	QString outputPath = datasetPath;
	// Now process results
	{
		/// Looking at pair-wise comparisons:
		std::cout << "path:" << datasetPath.toStdString() << std::endl;
		// Open labels JSON file
		QMap<QString, QStringList> lables;
		LabelOracle oracle;
		{
			auto labelsFilename = datasetPath + "/labels.json";
			QFile file;
			file.setFileName(labelsFilename);
			if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;
			QJsonDocument jdoc = QJsonDocument::fromJson(file.readAll());
			auto json = jdoc.object();

			// Get lables
			auto labelsArray = json["labels"].toArray();
			for (auto l : labelsArray)
			{
				auto label = l.toObject();
				auto parent = label["parent"].toString();

				lables[parent] << label["title"].toString();
			}

			// Get cross-labels
			auto crossLabelsArray = json["cross-labels"].toArray();

			// regular nodes
			for (auto k : lables.keys()) for (auto l : lables[k]) oracle.push(l, l);

			// cross-labeled
			for (auto l : crossLabelsArray)
			{
				auto crosslabel = l.toObject();
				oracle.push(crosslabel["first"].toString(), crosslabel["second"].toString());
			}

			oracle.build();
		}

		QElapsedTimer stats_timer; stats_timer.start();

		for (int pi = 0; pi < allMaps.size(); pi++)
		{
			auto corr = allMaps[pi];
			auto corrl = allMapsLabel[pi];


			QStringList sourceLabels, targetLabels;
			for (auto n : corrl) {
				sourceLabels << n.first;
				targetLabels << n.second;
			}

			// Build expected ground truth
			oracle.makeGroundTruth(sourceLabels, targetLabels);

			MatchingRecords M;

			for (int ci = 0; ci < corr.size(); ci++)
			{
				auto sn = corr[ci].first;
				auto tn = corr[ci].second;
				auto snl = corrl[ci].first;
				auto tnl = corrl[ci].second;

				M << MatchingRecord(sn, tn, snl, tnl);
			}

			oracle.pr_results << oracle.gt.compute(M);
		}

		// General stats
		int numG = 0, numM = 0, numR = 0;

		// Sum results
		double P = 0, R = 0;
		for (auto pr : oracle.pr_results)
		{
			P += pr.precision;
			R += pr.recall;

			numG += pr.G;
			numM += pr.M;
			numR += pr.R;
		}

		// Get average
		int N = oracle.pr_results.size();
		P /= N;
		R /= N;
		QString report = QString("[%3] Avg. P = %1, R = %2").arg(P).arg(R).arg(dir.dirName());

		report += QString("\nG_count %1 / M_count %2 / R_count %3").arg(numG).arg(numM).arg(numR);

		// Save log of P/R measures
		{
			QFile logfile(datasetPath + "/log.txt");
			logfile.open(QIODevice::WriteOnly | QIODevice::Text);
			QTextStream out(&logfile);
			out << report + "\n";
		}

		debugBox(report);
	}
}

Evaluator::~Evaluator()
{
    delete ui;
}
