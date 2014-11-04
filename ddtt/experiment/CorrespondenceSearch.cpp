#include "CorrespondenceSearch.h"
#include "DeformEnergy2.h"
#include "DeformEnergy.h"

#include <QProgressDialog>
QProgressDialog * pd = NULL;

CorrespondenceSearch::CorrespondenceSearch(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB, const Paths & paths, bool isUseAlternative)
	: paths(paths), shapeA(shapeA), shapeB(shapeB), isUseAlternative(isUseAlternative)
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

	pathDetails.clear();
	pathDetails.resize(paths.size());

    bool abort = false;

    // Evaluate correspondences
	#pragma omp parallel for
	for (int pi = 0; pi < paths.size(); pi++)
    {
		Structure::ShapeGraph shapeA_copy(*shapeA);
		Structure::ShapeGraph shapeB_copy(*shapeB);

        #pragma omp flush (abort)
        if(!abort)
        {
			// Check for global abort
            if(pd->wasCanceled()){
                abort = true;
                #pragma omp flush (abort)
            }

			auto & path = paths[pi];

            // Evaluate correspondence:
			if (!isUseAlternative)
			{
				DeformEnergy2 de(&shapeA_copy, &shapeB_copy, path.first, path.second, false);
				pathScores[pi] = de.total_energy;
				pathDetails[pi] = de.energyTerms;
            }
			else
			{
				DeformEnergy de(&shapeA_copy, &shapeB_copy, path.first, path.second, false);
				pathScores[pi] = de.total_energy;
				pathDetails[pi] = de.energyTerms;
			}

			// Report progress
            emit( pathComputed() );
        }
    }

    /// Assign best correspondence:
    {
		bestCorrespondence = std::min_element(pathScores.begin(), pathScores.end()) - pathScores.begin();
    }

    // Timing
    property["allSearchTime"].setValue( (int)allTimer.elapsed() );

    emit( done() );
}

void CorrespondenceSearch::increaseProgress()
{
    if(!pd->wasCanceled()) pd->setValue( pd->value() + 1 );
}
