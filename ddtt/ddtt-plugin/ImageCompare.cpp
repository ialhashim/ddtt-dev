#include <omp.h>

#include "ImageCompare.h"
#include <QDir>
#include <QFileInfo>
#include <QDebug>

#include "fft.h"
#include "QC.h"

ImageCompare::ImageCompare()
{
}

void ImageCompare::loadKnowledge(QString folderPath, QString datasetName)
{
	QDir d(folderPath);

	if(!d.exists()) {
		qDebug() << "Warning: Cannot find folder " << folderPath;
		return;
	}

	if(!datasetName.length()) datasetName = d.dirName();

	DataSet dataset;
	dataset.name = datasetName;
	dataset.path = folderPath;
	
	QStringList imageFiles = d.entryList( QStringList() << "*.png" );
	dataset.data.resize( imageFiles.size() );

	// Go over image files in the folder
	#pragma omp parallel for
	for(int i = 0; i < imageFiles.size(); i++)
	{
		QString file = imageFiles[i];
		QString fullFilename = d.absolutePath() + "/" + file;
		
		ImageCompare::Instance inst;
		inst.index = dataset.data.size();
		inst.filename = fullFilename;
		inst.id = QFileInfo(inst.filename).baseName();

		// Check for pre computed contour file
		QString contourFilename = inst.id + ".txt";
		if( d.exists(contourFilename) )
		{
			// Load from disk
			QFile contourFile( d.absolutePath() + "/" + contourFilename );
			QTextStream in( &contourFile );
			if( contourFile.open(QIODevice::ReadOnly | QIODevice::Text) )
			{
				while( !in.atEnd() )
				{
					QStringList data = in.readLine().split(" ", QString::SkipEmptyParts);
					if(data.isEmpty() || data.size() < 2) continue;

					inst.contour.push_back( std::make_pair(data[0].toDouble(), data[1].toDouble()) );
				}
			}
		}
		else
		{
			// Compute contour..
		}

		// Check for signatures
		QString sigFilename = inst.id + ".sig";
		if( d.exists(sigFilename) )
		{
			// Load from disk
			QFile sigFile( d.absolutePath() + "/" + sigFilename );
			QTextStream in( &sigFile );
			if( sigFile.open(QIODevice::ReadOnly | QIODevice::Text) )
			{
				QStringList signature = in.readLine().split(" ", QString::SkipEmptyParts);
				for(auto s : signature) inst.signature.push_back( s.toDouble() );
			}
		}
		else
		{
			// Compute signature..
			inst.signature = generateSignature(inst, true);
		}
		
		dataset.data[i] = inst;
	}

	datasets[datasetName] = dataset;
}

std::vector<double> ImageCompare::generateSignature( Instance inst, bool isSaveToFile )
{		
	int sig_length = 128;

	if(!inst.contour.size()) return std::vector<double>();

	// Find centroid
	std::pair<double,double> c(0,0);
	for(auto p : inst.contour) { c.first += p.first ; c.second += p.second; }
	c.first /= inst.contour.size(); c.second /= inst.contour.size();

	// Centroid Distance signature
	std::vector<double> dists;
	for(auto p : inst.contour) 
	{
		double dist = std::sqrt( pow(p.first - c.first,2) + pow(p.second - c.second,2) );
		dists.push_back( dist );
	}

	// FFT
	std::vector<double> real = dists, imag(real.size(), 0);
	Fft::transform(real, imag);

	bool isNormalizeFD = true;
	if( isNormalizeFD )
	{
		// Translation invariance
		real[0] = 0;
		imag[1] = 0;

		// Scale invariance
		double si = std::abs( std::complex<double>(real[1], imag[1]) );
		for(size_t i = 0; i < real.size(); i++)
		{
			real[i] /= si;
			imag[i] /= si;
		}

		// Rotation and changes in starting point
		for(size_t i = 0; i < real.size(); i++)
		{
			real[i] = std::abs( std::complex<double>(real[i], imag[i]) );
			imag[i] = 0;
		}
	}

	// Check against signature length
	if( inst.contour.size() < sig_length)
	{
		size_t padding = sig_length - inst.contour.size();
		size_t h = inst.contour.size() / 2;
		real.insert(real.begin() + h, padding, 0);
		imag.insert(imag.begin() + h, padding, 0);
	}

	// Flip
	Fft::fftshift(real, imag);

	// Shorten to assigned size 
	Fft::shrink(sig_length, real, imag);

	if( isSaveToFile )
	{
		// Save to disk
		QString sigFilename = inst.id + ".sig";
		QFile sigFile( QFileInfo(inst.filename).absolutePath() + "/" + sigFilename );
		QTextStream out( &sigFile );

		if( sigFile.open(QIODevice::WriteOnly | QIODevice::Text) )
			for(auto s : real) out << QString::number(s) << " ";
	}

	return real;
}

double ImageCompare::distance(const Instance &instanceA, const Instance &instanceB)
{
	/* Quadratic-Chi Histogram Distance */
	//double dist = QC::distance( j.signature, instance.signature );

	/* L-2 */
	double dist = QC::L2distance( instanceA.signature, instanceB.signature );	

	return dist;
}

QVector<ImageCompare::Instance> ImageCompare::kNearest(const ImageCompare::Instance &instance, int k, bool isReversed)
{
	QVector<ImageCompare::Instance> result;

	typedef QPair<double, int> ScoreInstance;
    QVector< ScoreInstance > candidates;

	for(auto key : datasets.keys())
	{
		DataSet & dataset = datasets[key];

		for(auto & j : dataset.data)
		{
			candidates << ScoreInstance( distance(j, instance), j.index );
		}
	
		// Sort...
		std::sort( candidates.begin(), candidates.end(), [](const ScoreInstance& a, const ScoreInstance& b){ return a.first < b.first; } );

		// To get 'k' furthest, useful for debugging
		if(isReversed) std::reverse(candidates.begin(), candidates.end());

		// Return first 'k'
		for(int i = 0; i < k; i++) result.push_back( datasets[key].data.at( candidates[i].second ) );
	}

    return result;
}

ImageCompare::Instance ImageCompare::getInstance( QString datasetName, int idx )
{
	return datasets[datasetName].data.at(idx);
}

size_t ImageCompare::datasetSize(QString datasetName)
{
	return datasets.value(datasetName).data.size();
}

QVector< QVector<ImageCompare::Instance> > ImageCompare::duplicateSets( double threshold, QString datasetName )
{
	if(!datasetName.length() || !datasets.keys().contains(datasetName))	datasetName = datasets.keys().front();
	DataSet & dataset = datasets[datasetName];

	int N = dataset.data.size();

	typedef QPair<double, int> ScoreInstance;
	QVector< QVector< ScoreInstance > > candidates(N);

	for(int i = 0; i < N; i++)
	{
		Instance & instance = dataset.data[i];

		for(int j = 0; j < N; j++)
		{
			if(i==j) continue;

			candidates[i].push_back( qMakePair(distance(instance, dataset.data.at(j)), j) );
		}

		// Sort
		std::sort( candidates[i].begin(), candidates[i].end(), 
			[](const ScoreInstance& a, const ScoreInstance& b){ return a.first < b.first; } );
	}

	QVector< QVector<ImageCompare::Instance> > result;

	for(int i = 0; i < N; i++)
	{
		std::vector<int> similar;
		for(auto p : candidates[i]) if(p.first < threshold) similar.push_back( p.second );

		QVector<ImageCompare::Instance> mysimilar;
		for(auto s : similar) mysimilar.push_back(dataset.data.at(s));

		if( !mysimilar.empty() ) 
		{
			mysimilar.push_front( dataset.data.at(i) );
			result.push_back( mysimilar );
		}
	}

	return result;
}

void ImageCompare::removeDuplicateSets( double threshold, QString datasetName )
{
	if(!datasetName.length() || !datasets.keys().contains(datasetName))	datasetName = datasets.keys().front();
	QString datasetPath = datasets[datasetName].path;

	QStringList filenamesToRemove;

	for(auto set : duplicateSets(threshold, datasetName)){
		set.removeFirst();
		for(auto instance : set) filenamesToRemove << instance.filename;	
	}
	if( filenamesToRemove.empty() ) return;

	if(QMessageBox::question(0, "Delete duplicates", QString("Are you sure you want to delete (%1) files?").arg(
		filenamesToRemove.size()), QMessageBox::Yes|QMessageBox::No) == QMessageBox::Yes){
			for(auto filename : filenamesToRemove)
				QFile::remove( filename );
	}

	// Reload
	loadKnowledge(datasetPath, datasetName);
}

void ImageCompare::showInstances(QVector<ImageCompare::Instance> instances)
{
	QVector<QImage> imgs;
	QStringList labels;

	for(auto i : instances){
		imgs.push_back( i.image() );
		labels.push_back( QString("%1 : %2").arg(i.index).arg(i.id) );
	}

	showImages(imgs, labels);
}
