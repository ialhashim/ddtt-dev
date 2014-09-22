#include "CorrespondenceSearch.h"
#include "PartCorresponder.h"
#include "myglobals.h"

#include <QProgressDialog>
QProgressDialog * pd = NULL;

CorrespondenceSearch::CorrespondenceSearch(CorrespondenceGenerator *generator) :
	generator( generator ), sA( generator->sA ), sB( generator->sB )
{
	int pathsCount = (int)generator->computedAssignments.size();
	property["pathsCount"].setValue( pathsCount );

    pd = new QProgressDialog( "Searching..", "Cancel", 0, pathsCount );
	pd->show();
	pd->setValue(0);

    this->connect( this, SIGNAL(pathComputed()), SLOT(increaseProgress()) );
}

void CorrespondenceSearch::run()
{
	QElapsedTimer allTimer; allTimer.start();

    auto & paths = generator->computedAssignments;

	QVector<ParticleMesh*> input; input << sA << sB;
	QVector<Particles> inputParticles; inputParticles << sA->particles << sB->particles;
	auto boxA = sA->bbox();

	std::vector<double> pathScores( paths.size(), 0 );
  
	bool abort = false;

    // Evaluate correspondences
    #pragma omp parallel for
    for(int pi = 0; pi < (int)paths.size(); pi++)
    {
        #pragma omp flush (abort)
        if(!abort)
        {
            if(pd->wasCanceled()){
                abort = true;
                #pragma omp flush (abort)
            }
        }

		double score = 0;
		auto & path = paths[pi];

		// Compute dense correspondence
		QVector<Particles> particles = inputParticles;
		for(auto & pairing : path) PartCorresponder::correspondSegments( pairing, input, particles );

		// Evaluate correspondence:
		{
			int numSamples = 4;
			int start = 1; // Skipping initial configuration

			for(int sample = start; sample < numSamples; sample++) 
			{
				double t = double(sample) / (numSamples-1);

				for(auto & p : particles[0])
				{
					// one-to-none
					if( p.correspondence > particles[1].size() ) 
					{
						score += (p.pos - boxA.center()).norm();
						continue;
					}

					auto blended = AlphaBlend(t, p.pos, particles[1][p.correspondence].pos);

					score += (blended - p.pos).norm();
				}
			}
		}

		pathScores[pi] = score;

		emit( pathComputed() );
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

void CorrespondenceSearch::applyCorrespondence(const PartsCorrespondence & correspondence)
{	
	QVector<ParticleMesh*> input; input << sA << sB;
	QVector<Particles> particles; 
	particles << sA->particles << sB->particles;

	for( auto & pairing : correspondence )
	{
		PartCorresponder::correspondSegments( pairing, input, particles );

		// DEBUG:
		if( true )
		{
			auto ls = new starlab::LineSegments(3);
			auto boxA = sA->segmentBoundingBox( input[0]->property["segments"].value<Segments>()[pairing.first] );
			auto boxB = sB->segmentBoundingBox( input[1]->property["segments"].value<Segments>()[pairing.second] );
			ls->addLine( Vector3( boxA.center() ), Vector3(boxB.center() + Vector3(1,0,0)), Qt::black );
			debug << ls;
		}
	}

	// Assign back
	sA->particles = particles.front();
	sB->particles = particles.back();
}
