#include "CorrespondenceSearch.h"
#include "DeformEnergy.h"

#include <QProgressDialog>
QProgressDialog * pd = NULL;

CorrespondenceSearch::CorrespondenceSearch(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB, const Paths & paths) 
	: paths(paths), shapeA(shapeA), shapeB(shapeB)
{
	property["pathsCount"].setValue(paths.size());

    pd = new QProgressDialog( "Searching..", "Cancel", 0, paths.size() );
    pd->show();
    pd->setValue(0);

    this->connect( this, SIGNAL(pathComputed()), SLOT(increaseProgress()) );
}

void CorrespondenceSearch::run()
{
    QElapsedTimer allTimer; allTimer.start();

	pathScores.clear();
    pathScores.resize( paths.size(), DBL_MAX );

    bool abort = false;

	// Make copy for thread safety
	QVector<Structure::ShapeGraph*> shapeAs(omp_get_max_threads());
	QVector<Structure::ShapeGraph*> shapeBs(omp_get_max_threads());
	for (size_t i = 0; i < shapeAs.size(); i++){
		shapeAs[i] = new Structure::ShapeGraph(*shapeA);
		shapeBs[i] = new Structure::ShapeGraph(*shapeB);
	}

    // Evaluate correspondences
	#pragma omp parallel for
	for (int pi = 0; pi < paths.size(); pi++)
    {
		auto & shapeA_copy = shapeAs[omp_get_thread_num()];
		auto & shapeB_copy = shapeBs[omp_get_thread_num()];

        #pragma omp flush (abort)
        if(!abort)
        {
            if(pd->wasCanceled()){
                abort = true;
                #pragma omp flush (abort)
            }

            double score = 0;
            auto & path = paths[pi];

            // Evaluate correspondence:
            {
				score = DeformEnergy(shapeA_copy, shapeB_copy, path.first, path.second, false).error;
            }

            pathScores[pi] = score;

            emit( pathComputed() );
        }
    }

    /// Assign best correspondence:
    {
        auto bestPath = std::min_element(pathScores.begin(), pathScores.end()) - pathScores.begin();
        bestCorrespondence = paths[bestPath];
    }

    // Timing
    property["allSearchTime"].setValue( (int)allTimer.elapsed() );

    emit( done() );
}

void CorrespondenceSearch::increaseProgress()
{
    if(!pd->wasCanceled()) pd->setValue( pd->value() + 1 );
}
