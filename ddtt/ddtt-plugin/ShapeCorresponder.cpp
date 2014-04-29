#include "ShapeCorresponder.h"
#include "PathsGenerator.h"

// Used for Evaluation
#include "SynthesisManager.h"
#include "Octree.h"
#include "SoftwareRenderer.h"
#include "ImageCompare.h"
#include "MarchingSquares.h"
#include "ContourMappingDistance.h"

ImageCompare im;
//QString chairsDatasetFolder = "C:/Temp/_imageSearch/test_data";
QString chairsDatasetFolder = "C:/Temp/_imageSearch/all_images_apcluster_data";

ShapeCorresponder::ShapeCorresponder(Structure::Graph * g1, Structure::Graph * g2, QString knowledge) : source(g1), target(g2)
{
    QElapsedTimer prepareTimer, computeTimer, evaluateTimer;
    prepareTimer.start();

	// Load knowledge
	im.loadKnowledge( chairsDatasetFolder, "chairs" );

    // Generate possible paths
    int k = 2;
    paths = PathsGenerator( source, target, k ).paths;

    // Timing & Stats
    property["pathsCount"].setValue( (int)paths.size() );
	property["prepareTime"].setValue( (int)prepareTimer.elapsed() );
	computeTimer.start();

	// Subsample paths	
	bool isSubsample = false;
    if( isSubsample ){
		int MaxNumPaths = 2;
        std::vector<bool> mask = subsampleMask( MaxNumPaths, paths.size() );
		std::vector<DeformationPath> subsampled;
		for(auto & path : paths){
            if(mask[path.i]) subsampled.push_back(path);
		}
		paths = subsampled; 
	}

	/// Cache Octrees
    QVector<Structure::Graph*> graphs; graphs << source << target;
    for(auto g : graphs){
		for(auto n : g->nodes){
			SurfaceMesh::Model * model = n->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >().data();
			if(!model) continue;
			Octree * octree = model->property("octree").value<Octree*>();
			if( !octree ){
				octree = new Octree(model, 40);
				QVariant oct; oct.setValue(octree);
				model->setProperty("octree", oct);
			}
		}
	}

	bool isVisualizeProcess = true;

	/// Find best correspondence
	int numSamples = 5;
	int k_neighbours = 4;

	DeformationPath bestPath;
	{
		beginFastNURBS();

		// Evaluate deformations
		#pragma omp parallel for
		for(int pi = 0; pi < (int)paths.size(); pi++)
		{
			auto & path = paths[pi];
			path.i = pi;

			// Prepare blending
			path.scheduler = QSharedPointer<Scheduler>( new Scheduler );
			path.blender = QSharedPointer<TopoBlender>( new TopoBlender( path.gcorr, path.scheduler.data() ) );
			path.synthman = QSharedPointer<SynthesisManager>( new SynthesisManager(path.gcorr, path.scheduler.data(), path.blender.data()) );

			// Deform
			path.scheduler->timeStep = 0.1;
			path.scheduler->executeAll();

			// Evaluate
			bool isEvaluate = true;
			if( isEvaluate )
			{
				// Prepare
				path.weight = 0.0;
				path.synthman->makeProxies(60, 20);

				// DEBUG:
				QImage pathImage;
				int pathImageWidth = 1280;
				if(isVisualizeProcess) pathImage = QImage( pathImageWidth, pathImageWidth * 0.6, QImage::Format_ARGB32_Premultiplied );

				// Sample the computed path
				for(int s = 0; s < numSamples; s++)
				{
					double t = double(s)/(numSamples-1);
					t = ((1.0 - 0.6) / 2.0) + (t * 0.6); // middle 60%

					int idx = t * (path.scheduler->allGraphs.size()-1);
					Structure::Graph * g = path.scheduler->allGraphs[idx];

					std::vector<SimplePolygon> polys = path.synthman->drawWithProxies(g);

					if( polys.size() )
					{
						QVector< QVector<Vector3> > tris;

						for(auto p : polys){
							if(p.vertices.size() < 4)
								tris.push_back( QVector<Vector3>() << p.vertices[0] << p.vertices[1] << p.vertices[2] );
							else{
								tris.push_back( QVector<Vector3>() << p.vertices[0] << p.vertices[1] << p.vertices[2] );
								tris.push_back( QVector<Vector3>() << p.vertices[2] << p.vertices[3] << p.vertices[0] );
							}
						}

						// Setup camera
						Eigen::AlignedBox3d graphBBox = g->bbox();
						double distance = graphBBox.sizes().maxCoeff() * 4.5;
						Vector3 direction (-2,-2,0.8);
						direction.normalize();
						Vector3 target = graphBBox.center();
						Vector3 eye = (direction * distance) + target;
						Vector3 up(0,0,1);
						Eigen::MatrixXd camera = SoftwareRenderer::CreateViewMatrix(eye, target, up);

						// Draw geometry into a buffer
						Eigen::MatrixXd buffer = SoftwareRenderer::render(tris, 128, 128, camera);

						// Find outer most contour using marching squares
						std::vector< std::pair<double,double> > contour;
						for(auto p : MarchingSquares::march(buffer, 1.0))
							contour.push_back( std::make_pair(p.x(), p.y()) );

						// Compute signature
						ImageCompare::Instance sampleInstance( contour );
						sampleInstance.index = s;
						sampleInstance.id = QString("path-%1").arg(pi);

						// Compare with respect to knowledge
						ImageCompare::InstanceMatches neighbours = im.kNearest(sampleInstance, k_neighbours);

						for( auto against : neighbours )
						{
							/* Signature distance */
							//minDistance = std::min(minDistance, against.first );

							/* Contour Mapping Measure */
							std::vector<Eigen::Vector2d> cA, cB;
							for(auto p : sampleInstance.contour) cA.push_back(Eigen::Vector2d(p.first, p.second));
							for(auto p : against.second.contour) cB.push_back(Eigen::Vector2d(p.first, p.second));

							path.weight += optMapMae(cA,cB);
						}

						/// DEBUG:
						if( isVisualizeProcess )
						{
							sampleInstance.cachedImage = SoftwareRenderer::matrixToImage( buffer, false );
							sampleInstance.cachedImage = ImageCompare::visualizeInstance(sampleInstance, QString::number(path.weight));

							QPainter painter(&pathImage);

							int w = pathImageWidth / numSamples;
							painter.drawImage(s * w, 0, sampleInstance.cachedImage);
							painter.drawRect(sampleInstance.cachedImage.rect().translated(QPoint(s*w,0)));

							int top = sampleInstance.cachedImage.height(), left = s * w, padding = 5;
							int x = left, y = top + padding;

							for(int ni = 0; ni < (int)neighbours.size(); ni++)
							{
								QImage img = ImageCompare::visualizeInstance(neighbours[ni].second, QString::number(neighbours[ni].first));
								img = img.scaledToWidth(w * 0.47, Qt::SmoothTransformation);
								painter.drawImage(x, y, img);

								x += img.width(); 
								if(x > (w + left) - img.width()) { x = left; y += img.height() + padding; }
							}
						}
					}
				}

				if( isVisualizeProcess && !pathImage.isNull() )
				{
					QString filename = QString("path-%1.png").arg(pi);
					pathImage.save( filename );
				}
			}

			// Clean up
			{
				path.scheduler.clear();
				path.blender.clear();
				path.synthman.clear();
				path.errors.clear();
			}
		}

		endFastNURBS();
	}

	// Timing
	property["computeTime"].setValue( (int)computeTimer.elapsed() );

	if( !paths.size() ) return;

	// Find best
	{
		// Best = lowest error
		std::sort(paths.begin(), paths.end(), DeformationPathCompare);
		bestPath = paths.front();

		// set indices
		int j = 0; for( auto & p : paths )	p.idx = j++;

		// Color with best correspondence
		source->setColorAll(Qt::lightGray);
		target->setColorAll(Qt::lightGray);

		for( auto p : bestPath.pairs )
		{
			QColor c = starlab::qRandomColor2();
			for( auto sid : p.first ) source->setColorFor(sid, c);
			for( auto tid : p.second ) target->setColorFor(tid, c);
		}
	}
}
