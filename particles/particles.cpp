#include "particles.h"
#include "particles-widget.h"
#include "ui_particles-widget.h"
#include <omp.h>
#include <numeric>
#include <cmath>

#include "SurfaceMeshModel.h"
#include "SurfaceMeshHelper.h"
using namespace SurfaceMesh;

QTimer * timer = NULL;

QStringList files;

#include "ParticleMesh.h"
#include "myglobals.h"
#include "Bounds.h"

#include "StructureAnalysis.h"

#include "Raytracing.h"

#include "spherelib.h"
#include "SphericalHarmonic.h"
std::vector<Spherelib::Sphere*> spheres;

#include "BasicTable.h"

#include "kmeans.h"

#include "ParticleCorresponder.h"
#include "ParticleDeformer.h"

void particles::processShapes()
{
	if(!widget) return;

	ParticlesWidget * pw = (ParticlesWidget *) widget;
	pw->isReady = false;

	if(pw->pmeshes.empty()) return;

	drawArea()->clear();

	// Check grid size changes, redo voxelization
	if( pw->pmeshes.front()->grid.gridsize != pw->ui->gridsize->value() )
		reVoxelize();

	// Show voxelized meshes
	//for(auto & s : pw->pmeshes) document()->addModel(s->surface_mesh->clone());

	if(pw->pmeshes.size() > 1)
	{
		// Experiment
		blending();
		pw->isReady = true;
		return;
	}

	QElapsedTimer allTimer; allTimer.start();

	qApp->processEvents();

	// Uniform sampling on the sphere
	Spherelib::Sphere sphere( pw->ui->sphereResolution->value() );
	std::vector< Eigen::Vector3d > sampledRayDirections  = sphere.rays();
	std::vector< size_t > antiRays = sphere.antiRays();
	size_t perSampleRaysCount = sampledRayDirections.size();

	// Rotation invariant descriptor
	SphericalHarmonic<Vector3,float> sh( std::max(1,pw->ui->bands->value()) );
	std::vector< SHSample<Vector3,float> > sh_samples;
	sh.SH_setup_spherical( sampledRayDirections, sh_samples );

	if( true )
	{
		for(auto & s : pw->pmeshes)
		{
			mainWindow()->setStatusBarMessage(QString("Working with (%1) particles.").arg(s->particles.size()));
			s->debug.clear();

			QElapsedTimer timer; timer.start();

			// Keep record of sampled directions
			s->usedDirections = sampledRayDirections;
			s->antiRays = antiRays;
			s->sphereResolutionUsed = pw->ui->sphereResolution->value();

			// Hit results are saved as vector of distances
			std::vector< std::vector<float> > & descriptor = (s->desc = std::vector< std::vector<float> >( 
				s->particles.size(), std::vector<float>( perSampleRaysCount ) ));

			int rayCount = 0;

			// Accelerated raytracing
			raytracing::Raytracing<Eigen::Vector3d> rt( s->surface_mesh );

			// Shooting rays from inside the volume
			{
				// Medial points record
				std::vector<Eigen::Vector3f> ma_point(s->particles.size(), Eigen::Vector3f(0,0,0));
				std::vector<bool> ma_point_active(s->particles.size(), false);
				std::vector<float> ma_point_rad(s->particles.size(), 0);

				// Accelerated raytracing
				//raytracing::Raytracing<Eigen::Vector3d> rt( s->surface_mesh );
				// ^^^ Moved up ^^

				// Shoot rays around all particles
				#pragma omp parallel for
				for(int pi = 0; pi < (int)s->particles.size(); pi++)
				{
					auto & p = s->particles[pi];

					int r = 0;
					for(auto d : sampledRayDirections)
					{
						raytracing::RayHit hit = rt.hit( p.pos, d );
						descriptor[pi][r++] = hit.distance;
						rayCount++;
					}

					p.avgDiameter = 0;
				}

				// Smooth ray response
				int smoothRaysIter = pw->ui->fnSmoothIters->value();
				if(smoothRaysIter > 0)
				{
					for(int pi = 0; pi < (int)s->particles.size(); pi++)
					{
						auto & p = s->particles[pi];

						sphere.setValues( descriptor[pi] );
						sphere.smoothValues( smoothRaysIter );
						descriptor[pi] = sphere.values();
					}
				}

				// Report
				mainWindow()->setStatusBarMessage( QString("Ray tracing: rays (%1) / time (%2 ms)").arg( rayCount ).arg( timer.elapsed() ) );
				timer.restart();

				// Average of neighbors
				int avgNeighIters = pw->ui->avgNeighIters->value();
				if( avgNeighIters )
				{
					size_t gridsize = s->grid.gridsize;

					for(int it = 0; it < avgNeighIters; it++)
					{
						auto smoothDesc = descriptor;
						{
							#pragma omp parallel for
							for(int pi = 0; pi < (int)s->particles.size(); pi++)
							{
								const auto & p = s->particles[pi];

								std::vector<float> sum( descriptor[p.id].size(), 0 );
								int count = 0;

								unsigned int x,y,z;
								mortonDecode(p.morton, x, y, z);

								for(int u = -1; u <= 1; u++){
									for(int v = -1; v <= 1; v++){
										for(int w = -1; w <= 1; w++)
										{
											Eigen::Vector3i c(x + u, y + v, z + w);

											// Skip outside grid
											if(c.x() < 0 || c.y() < 0 || c.z() < 0) continue;
											if(c.x() > gridsize-1 || c.y() > gridsize-1 || c.z() > gridsize-1) continue;

											uint64_t mcode = mortonEncode_LUT(c.x(),c.y(),c.z());
											if(p.morton == mcode) continue;
											if(!s->grid.occupied[ mcode ]) continue;

											size_t nj = s->mortonToParticleID[ mcode ];

											std::vector<float> nei = descriptor[nj];

											// Add this neighbor to sum
											for(size_t j = 0; j < sum.size(); j++) sum[j] += nei[j];
											count++;
										}
									}
								}

								// divide over count
								for(size_t j = 0; j < sum.size(); j++) sum[j] /= count;

								smoothDesc[p.id] = sum;
							}
						}
						descriptor = smoothDesc;
					}
				}

				// Use medial points
				if( pw->ui->useMedialPoints->isChecked() )
				{
					// Compute medial points
					{
						// Points just outside the volume
						NanoKdTree surface_kdtree;
						for(auto p : s->grid.pointsOutside()) surface_kdtree.addPoint(p.cast<double>());
						surface_kdtree.build();

						// Find maximal balls inside shape
						#pragma omp parallel for
						for(int pi = 0; pi < (int)s->particles.size(); pi++)
						{
							auto & p = s->particles[pi];
							std::vector<float> & desc = descriptor[pi];

							Vector3 maPoint = p.pos;
							
							double sumDiameter = 0;
							int visitCount = 0;
							std::vector< bool > isVisisted( sampledRayDirections.size(), false );

							typedef std::pair<double,Vector3> RadiusCenter;
							std::vector<RadiusCenter> medial_balls;

							double diameter, radius, d2;
							size_t ret_index;
							bool isFullyInside;

							for(size_t idx = 0; idx < sampledRayDirections.size(); idx++)
							{
								if(isVisisted[antiRays[idx]]) continue;
								isVisisted[antiRays[idx]] = true;
								visitCount++;

								Vector3 start( p.pos + sampledRayDirections[idx] * desc[idx] );
								Vector3 end( p.pos + sampledRayDirections[antiRays[idx]] * desc[antiRays[idx]] );
								Vector3 midpoint = (start + end) * 0.5;

								diameter = desc[idx] + desc[antiRays[idx]];
								radius = diameter / 2;

								sumDiameter += diameter;

								surface_kdtree.tree->knnSearch(&midpoint[0], 1, &ret_index, &d2);
								
								isFullyInside = d2 > (radius * radius);

								if( isFullyInside )
									medial_balls.push_back(RadiusCenter(radius,midpoint));
							}

							// Select medial point
							if( medial_balls.size() )
							{
								// Sort balls by radius
								std::sort(medial_balls.begin(),medial_balls.end(), 
									[](RadiusCenter a, RadiusCenter b){ return a.first < b.first; });
								auto selected_medial = medial_balls.front();
								
								ma_point[pi] = selected_medial.second.cast<float>();
								ma_point_rad[pi] = selected_medial.first;
								ma_point_active[pi] = true;

								// Record average diameter around particle
								p.avgDiameter = sumDiameter / visitCount;
							}
						}
					}

					NanoKdTree all_medial_kdtree;
					for(size_t i = 0; i < ma_point_active.size(); i++)
						if( ma_point_active[i] )
							all_medial_kdtree.addPoint( ma_point[i].cast<double>() );
					all_medial_kdtree.build();

					// Filter out not so medial points
					if( pw->ui->filterMedialPoints->isChecked() )
					{
						#pragma omp parallel for
						for(int pi = 0; pi < (int)s->particles.size(); pi++)
						{
							float ma_rad = ma_point_rad[pi];

							// Either thin area or minor feature
							bool isPossibleThin = ma_rad < s->grid.unitlength * 2;

							if( ma_point_active[pi] && isPossibleThin )
							{
								float max_rad = ma_point_rad[pi];

								// Look around for access to wider regions
								for(auto pj : s->neighbourhood( s->particles[pi] ))
									max_rad = std::max(max_rad, ma_point_rad[pj]);

								// A neighbor has wider access than me => I'm not so medial
								double scale = 1.2;
								if( max_rad > ma_point_rad[pi] * scale ){
									ma_point_active[pi] = false;
								}
							}
						}
					}

					// KD-tree of medial points (with map)
					NanoKdTree medial_kdtree;
					std::map<size_t,size_t> ma_particle;
					for(size_t i = 0; i < ma_point_active.size(); i++){
						if( ma_point_active[i] ){
							ma_particle[ma_particle.size()] = i;
							medial_kdtree.addPoint( ma_point[i].cast<double>() );
						}
					}
					medial_kdtree.build();

					// Report
					mainWindow()->setStatusBarMessage( QString("Medial particles (%1 ms)").arg( timer.elapsed() ) );
					timer.restart();

					// Only use descriptors from medial points
					if( true )
					{
						// Compute medial descriptors
						#pragma omp parallel for
						for(int pi = 0; pi < (int)s->particles.size(); pi++)
						{
							if( ma_point_active[pi] ){
								Vector3 maPoint = ma_point[pi].cast<double>();
								int r = 0;
								for(auto d : sampledRayDirections){
									raytracing::RayHit hit = rt.hit( maPoint, d );
									descriptor[pi][r++] = hit.distance;
								}

								// Compute flatness of my neighbourhood
								double search_rad = s->grid.unitlength * 4;
								KDResults matches;
								all_medial_kdtree.ball_search(maPoint, search_rad, matches);

								Eigen::MatrixXf x(3, matches.size());
								for(size_t i = 0; i < matches.size(); i++)
									x.col(i) = all_medial_kdtree.cloud.pts[matches[i].first].cast<float>();

								Eigen::MatrixXf avg = x.rowwise().mean();
								x = x - avg.replicate(1,x.cols());
								Eigen::MatrixXf sigma = x * x.transpose() * (1.0 / x.cols());
								Eigen::JacobiSVD<Eigen::MatrixXf> svd(sigma, Eigen::ComputeFullU | Eigen::ComputeFullV);
								Eigen::Vector3f values = svd.singularValues().normalized();
								Eigen::Matrix3f V = svd.matrixV();

								auto normal = V.col(0).cross( V.col(1) ).normalized();
								float ratio = (values[1] / values[0]);
								
								s->particles[pi].axis = normal.cast<double>();
								s->particles[pi].flat = ratio;

								ma_point_rad[pi] = ratio;

								s->particles[pi].medialPos = maPoint;
							}
						}

						#pragma omp parallel for
						for(int pi = 0; pi < (int)s->particles.size(); pi++)
						{
							if( !ma_point_active[pi] )
							{
								size_t pj = ma_particle[medial_kdtree.closest( s->particles[pi].pos )];
								descriptor[pi] = descriptor[pj];

								s->particles[pi].avgDiameter = s->particles[pj].avgDiameter;
								//s->particles[pi].measure = s->particles[pj].measure;
								s->particles[pi].axis = s->particles[pj].axis;
								s->particles[pi].flat = s->particles[pj].flat;

								s->particles[pi].medialID = pj;
								s->particles[pi].isMedial = false;
								s->particles[pi].medialPos = s->particles[pj].medialPos;
							}
							else
							{
								s->particles[pi].medialID = pi;
								s->particles[pi].isMedial = true;
							}

							auto maxelement = std::max_element(descriptor[pi].begin(),descriptor[pi].end());
							s->particles[pi].direction = sampledRayDirections[ maxelement - descriptor[pi].begin() ];
						}
					}

					// Normalize diameters
					if( true )
					{
						Bounds<float> b;
						for(auto & p : s->particles) b.extend(p.avgDiameter);
						for(auto & p : s->particles) p.avgDiameter = b.normalized(p.avgDiameter);
					}

					// Report
					mainWindow()->setStatusBarMessage( QString("Projection to Medial particles (%1 ms)").arg( timer.elapsed() ) );
					timer.restart();

					// [DEBUG] medial points
					if( pw->ui->showMedial->isChecked() )
					{
						starlab::PointSoup * ps = new starlab::PointSoup;
						Boundsf interval;

						for(size_t i = 0; i < ma_point_active.size(); i++)
							if(ma_point_active[i])
								interval.extend(ma_point_rad[i]);

						for(size_t i = 0; i < ma_point_active.size(); i++) {
							if(ma_point_active[i]){
								float val = interval.normalized( ma_point_rad[i] );
								ps->addPoint(Vector3(ma_point[i].cast<double>()), starlab::qtJetColor(val));
							}
						}

						drawArea()->addRenderObject(ps);

						return;
					}

					// [DEBUG] show projected skeleton
					if( pw->ui->projectSkeleton->isChecked() )
					{
						SurfaceMeshModel * m = s->surface_mesh->clone();
						auto points = m->vertex_coordinates();
						for(auto v : m->vertices())
							points[v] = medial_kdtree.cloud.pts[ medial_kdtree.closest(points[v]) ];

						document()->addModel( m );
						drawArea()->setRenderer(m, "Flat Wire");
					}

					// [DEBUG] show main 'direction' of particles
					if( false )
					{
						starlab::LineSegments * vs = new starlab::LineSegments(2);
						for(auto & particle : s->particles){
							Vector3 vec = particle.direction;
							Eigen::Vector4d color(abs(vec[0]), abs(vec[1]), abs(vec[2]), particle.alpha);
							color[0] *= color[0];color[1] *= color[1];color[2] *= color[2];
							vs->addLine(particle.pos,  Vector3(particle.pos + vec * 0.01), QColor::fromRgbF(color[0],color[1],color[2],1));
						}
						drawArea()->addRenderObject(vs);
						return;
					}
				}
			}

			// Rotation invariant descriptor
			if( pw->ui->useRotationInv->isChecked() && pw->ui->bands->value() > 0 )
			{
				s->sig = std::vector< std::vector<float> >(s->particles.size(), std::vector<float>( pw->ui->bands->value() ));

				#pragma omp parallel for
				for(int pi = 0; pi < (int)s->particles.size(); pi++){
					auto & p = s->particles[pi];
					auto desc = s->desc[pi];
					//desc = Bounds<float>::from( desc ).normalize( desc ); // Normalize descriptor
					std::vector<float> coeff;
					sh.SH_project_function(desc, sh_samples, coeff);
					s->sig[p.id] = sh.SH_signature(coeff);
				}
			}

			// [DEBUG] visualize distance to ground
			if( pw->ui->showDistGround->isChecked() )
			{
				for(auto & particle : s->particles)
					drawArea()->drawPoint(particle.pos, 5, starlab::qtJetColor(particle.measure));
				return;
			}

			// [DEBUG] visualize average diameter
			if( pw->ui->showDiameter->isChecked() )
			{
				for(auto & particle : s->particles)
					drawArea()->drawPoint(particle.pos, 5, starlab::qtJetColor(particle.avgDiameter));
				return;
			}

			// Standardize features column-wise to zero mean and unit variance.
			{
				//size_t n = s->particles.size();
				//size_t p = s->desc.front().size();
				//auto M = toEigenMatrix<float>(s->desc);
				//M.normalize();
				//s->desc = fromEigenMatrix<float>(M);
			}
			
		}
	}

	// Report
	mainWindow()->setStatusBarMessage( QString("Descriptors ready (%1 ms)").arg( allTimer.elapsed() ) );

	// Descriptor options
	for(auto & s : pw->pmeshes)
	{
		#pragma omp parallel for
		for(int i = 0; i < (int)s->particles.size(); i++)
		{
			auto & p = s->particles[i];

			std::vector<float> new_desc;

			if(pw->ui->useDescriptor->isChecked() ) 	new_desc = s->desc[p.id];
			if(pw->ui->useRotationInv->isChecked())		new_desc = s->sig[p.id];
			if(pw->ui->useGroundDist->isChecked() )		new_desc.push_back(p.measure);
			if(pw->ui->useDiameter->isChecked()   ) 	new_desc.push_back(p.avgDiameter);
			if(pw->ui->useFlat->isChecked()		  )		new_desc.push_back(p.flat);
			if(pw->ui->useHeight->isChecked()     )		new_desc.push_back(p.pos.z());

			if(new_desc.empty()) new_desc.push_back(p.pos.z()); // simply height..

			// Numerical check
			for(auto & d : new_desc) if(isnan(d) || !isfinite(d)) d = 0;

			s->desc[p.id] = new_desc;
		}

		//showTable(s->desc, std::min(size_t(100),s->particles.size()));
	}
	
	// Pre-segmented shapes case
	if( pw->ui->isSegmentedMesh->isChecked() ) 
	{
		pw->isReady = true;
		drawArea()->update();
		emit( shapesProcessed() );
		mainWindow()->setStatusBarMessage( QString("All time (%1 ms)").arg( allTimer.elapsed() ) );
		return;
	}

	QElapsedTimer curTimer; curTimer.start();

	// k-means clustering
	if( !pw->ui->performSturctureAnalysis->isChecked() )
	{
		int K = pw->ui->kclusters->value();

		for(auto & s : pw->pmeshes)
		{
			std::vector<size_t> seeds;

			if(pw->ui->useGroundDistSeed->isChecked())	seeds = s->specialSeeding(ParticleMesh::GROUND, K);
			if(pw->ui->useDescSeed->isChecked())		seeds = s->specialSeeding(ParticleMesh::DESCRIPTOR, K);

			s->cluster(K, seeds, pw->ui->l1norm->isChecked(), pw->ui->showSeeds->isChecked());
		}
	}

	// Report
	mainWindow()->setStatusBarMessage( QString("Clustering (%1 ms)").arg( curTimer.elapsed() ) );
	curTimer.restart();

	// Merge smaller clusters with larger ones
	for(int i = 0; i < pw->ui->simplifyCluster->value(); i++)
	{
		for(auto & s : pw->pmeshes)
			s->shrinkSmallerClusters();

		// Report
		mainWindow()->setStatusBarMessage( QString("Merging clusters (%1 ms)").arg( curTimer.elapsed() ) );
		curTimer.restart();
	}

	// Analysis
	if( pw->ui->performSturctureAnalysis->isChecked() )
	{	
		mainWindow()->setStatusBarMessage( QString("Unsupervised clustering..") );

		for(auto & s : pw->pmeshes)
		{		
			s->property["showHulls"].setValue( pw->ui->showHulls->isChecked() );
			s->property["isMerge"].setValue( pw->ui->isMerge->isChecked() );

			StructureAnalysis sa( s );

			for(auto d : sa.debug) drawArea()->addRenderObject(d);
		}

		// Report
		mainWindow()->setStatusBarMessage( QString("Structure analysis (%1 ms)").arg( curTimer.elapsed() ) );
		curTimer.restart();
	}

	// Report
	mainWindow()->setStatusBarMessage( QString("All time (%1 ms)").arg( allTimer.elapsed() ) );

	pw->isReady = true;

	drawArea()->update();

	emit( shapesProcessed() );
}

void particles::blending()
{
	if(!widget) return;
	ParticlesWidget * pw = (ParticlesWidget *) widget;
	pw->isReady = false;

	if(pw->pmeshes.size() < 2) return;
	ParticleCorresponder pc(pw->pmeshes.front(), pw->pmeshes.back());
	ParticleDeformer pd(pw->pmeshes.front(), pw->pmeshes.back());
}

void particles::reVoxelize()
{	
	if(!widget) return;

	ParticlesWidget * pw = (ParticlesWidget *) widget;
	pw->pmeshes.clear();
	for(auto filename : files){
		SurfaceMeshModel fromMesh( filename, QFileInfo(filename).baseName() );
		fromMesh.read( filename.toStdString() );
		// Rotate if requested
		double angleDeg = pw->ui->rotateAngle->value();
		if( angleDeg ){
			Eigen::AngleAxisd rot( deg_to_rad(angleDeg), Vector3(0,0,1));
			for(auto v : fromMesh.vertices()){
				auto & p = fromMesh.vertex_coordinates()[v];
				p = rot * p;
			}
		}
		pw->pmeshes.push_back( new ParticleMesh( &fromMesh, pw->ui->gridsize->value() ) );
	}
}

void particles::create()
{
	if( widget ) return;

	// Viewer
	{
		//drawArea()->setAxisIsDrawn(true);
		drawArea()->camera()->setType(qglviewer::Camera::PERSPECTIVE);

		double worldRadius = 1;
		drawArea()->camera()->setUpVector(qglviewer::Vec(0,0,1));
		drawArea()->camera()->setPosition(qglviewer::Vec(2,-2,1.5));
		drawArea()->camera()->lookAt(qglviewer::Vec());
		drawArea()->camera()->setSceneRadius( worldRadius );
		drawArea()->camera()->showEntireScene();
	}

	ModePluginDockWidget * dockwidget = new ModePluginDockWidget("Particles", mainWindow());
	ParticlesWidget * pw = new ParticlesWidget();
	widget = pw;

	dockwidget->setWidget( widget );
	mainWindow()->addDockWidget(Qt::RightDockWidgetArea, dockwidget);

	// Collect experiment figures
	{
		QWidget * experimentViewer = new QWidget;
		QVBoxLayout * expViewerLayout = new QVBoxLayout;

		QWidget * thumbsWidget = new QWidget;
		QVBoxLayout * thumbsLayout = new QVBoxLayout;

		QScrollArea * scrollarea = new QScrollArea;
		scrollarea->setBackgroundRole(QPalette::Dark);
		scrollarea->setWidgetResizable(true);

		int thumbWidth = 400;

		connect(this, &particles::shapesProcessed, [=](){
			QLabel * thumb = new QLabel;
			thumb->setPixmap(QPixmap::fromImage(drawArea()->grabFrameBuffer().scaledToWidth(thumbWidth)));
			thumbsWidget->layout()->addWidget(thumb);

			// Scroll
			QTimer *timer = new QTimer(this);
			timer->setSingleShot(true);
			connect(timer, &QTimer::timeout, [=]() {
				auto vbar = scrollarea->verticalScrollBar();
				vbar->setValue(vbar->maximum() - thumb->height());
				timer->deleteLater();
			} );
			timer->start(100);
		});

		thumbsWidget->setLayout( thumbsLayout );
		scrollarea->setWidget( thumbsWidget );
		expViewerLayout->addWidget(scrollarea);
		experimentViewer->setLayout(expViewerLayout);
		experimentViewer->setMinimumSize(thumbWidth,thumbWidth);
		experimentViewer->move( mainWindow()->geometry().topLeft() - QPoint(thumbWidth,0) );
		experimentViewer->show();
	}

	// General Tests
	connect(pw->ui->testButton, &QPushButton::released, [=]
	{
		auto s = pw->pmeshes.front();
		SegmentGraph neiGraph;
		auto segments = s->segmentToComponents( s->toGraph(), neiGraph );

		int id = 0;

		drawArea()->clear();

		for(auto & seg : segments)
		{
			// Collect and map current particles
			std::vector< std::vector<float> > descs;
			QMap<size_t,size_t> pmap;

			for(auto v : seg.vertices){
				pmap[v] = pmap.size();
				descs.push_back( s->desc[v] );
			}

			// Perform binary split
			int K = 2;
			clustering::kmeans< std::vector< std::vector<float> >, clustering::lpnorm< std::vector<float> > > km( descs, K );
			km._centers.clear();
			auto seeds = s->specialSeeding(ParticleMesh::DESCRIPTOR, K, seg.vertices);

			for(auto pid : seeds) km._centers.push_back( s->desc[pid] );
			km.run();

			// Assign found clusters
			for(auto v : seg.vertices)
				s->particles[v].segment = id + km.cluster( pmap[v] );

			// show seeds
			{
				auto color = ParticleMesh::rndcolors.at( (double(rand()) / RAND_MAX) * (ParticleMesh::rndcolors.size()-1) );
				starlab::PointSoup * ps = new starlab::PointSoup(20);
				if(!seeds.empty()) for(auto pid : seeds) ps->addPoint(s->particles[pid].pos, color);
				else for(auto pid : km.initindices) ps->addPoint(s->particles[pid].pos, color);
				drawArea()->addRenderObject(ps);
			} 

			id += 2;
		}

		drawArea()->update();
	});

	// Load and process shapes:
	connect(pw->ui->loadShapes, &QPushButton::released, [=]{
		files = QFileDialog::getOpenFileNames(mainWindow(), "Open Shapes", "", "All Supported (*.obj *.off)");

		QElapsedTimer timer; timer.start();

		for(auto filename : files)
		{
			if(pw->ui->isSegmentedMesh->isChecked())
			{
				pw->pmeshes.push_back( new ParticleMesh( filename, pw->ui->gridsize->value() ) );
			}
			else
			{
				SurfaceMeshModel fromMesh( filename, QFileInfo(filename).baseName() );
				fromMesh.read( filename.toStdString() );

				// Rotate if requested
				double angleDeg = pw->ui->rotateAngle->value();
				if( angleDeg ){
					Eigen::AngleAxisd rot( deg_to_rad(angleDeg), Vector3(0,0,1));
					for(auto v : fromMesh.vertices()){
						auto & p = fromMesh.vertex_coordinates()[v];
						p = rot * p;
					}
				}

				pw->pmeshes.push_back( new ParticleMesh( &fromMesh, pw->ui->gridsize->value() ) );
			}
		}

		mainWindow()->setStatusBarMessage(QString("Shapes loaded and voxelized (%1 ms)").arg(timer.elapsed()));

		emit( pw->shapesLoaded() );
	});

	connect(pw->ui->clearMeshes, &QPushButton::released, [=]{
		pw->pmeshes.clear();
		drawArea()->update();
	});

	connect(pw->ui->saveMeshes, &QPushButton::released, [=]{
		for(auto pmesh : pw->pmeshes){
			QString filename = QFileDialog::getSaveFileName(mainWindow(),"Save mesh", "", "Particle meshes (*.pmesh)");
			if(filename.length() < 1) continue;
			QFile file( filename );
			file.open(QIODevice::WriteOnly);
			QDataStream out(&file);
			out << *pmesh;
		}
	});

	connect(pw->ui->loadMeshes, &QPushButton::released, [=]{
		pw->pmeshes.clear();
		drawArea()->clear();

		for(auto filename : QFileDialog::getOpenFileNames(mainWindow(), "Load mesh", "", "Particle meshes (*.pmesh)"))
		{
			QFile file( filename );
			file.open(QIODevice::ReadOnly);
			QDataStream in(&file);

			auto pmesh = new ParticleMesh();
			in >> *pmesh;
			pw->pmeshes.push_back(pmesh);
		}

		if(pw->pmeshes.empty()) return;

		if(pw->pmeshes.size() > 1) 
		{
			QElapsedTimer processingTimer; processingTimer.start();

			ParticleCorresponder pc(pw->pmeshes.front(), pw->pmeshes.back());
			for(auto d : pc.debug) drawArea()->addRenderObject(d);

			ParticleDeformer pd(pw->pmeshes.front(), pw->pmeshes.back());
			for(auto d : pd.debug) drawArea()->addRenderObject(d);

			mainWindow()->setStatusBarMessage(QString("Processed (%1 ms)").arg(processingTimer.elapsed()));
		}

		pw->isReady = true;
		drawArea()->update();
	});

	// Post-processing
	connect(pw, SIGNAL(shapesLoaded()), SLOT(processShapes()));
	connect(pw->ui->processShapesButton, SIGNAL(clicked()), SLOT(processShapes()));

#ifdef QT_DEBUG
	pw->ui->gridsize->setValue(16);
#endif
}

void particles::decorate()
{
	ParticlesWidget * pwidget = (ParticlesWidget*) widget;
	if(!pwidget || !pwidget->isReady || pwidget->pmeshes.size() < 1) return;

	for(auto & sphere : spheres)
	{
		sphere->draw();
	}

	bool isDebugVisualization = !pwidget->ui->isShowBlend->isChecked();
	if(pwidget->pmeshes.size() > 1 && pwidget->pmeshes.front()->property["debug"].toBool())
		isDebugVisualization = true;

	if( pwidget->pmeshes.size() < 2 || isDebugVisualization ) 
	{
		glPushMatrix();

		// Evaluation
		for( auto s : pwidget->pmeshes )
		{
			s->drawParticles( drawArea()->camera() );
			s->drawDebug( *drawArea() );

			glTranslated(1, 0, 0);
		}

		glPopMatrix();

		return;
	}

	// Experimental
	static bool isForward = true;
	static double alpha = 0;
	if(isForward) alpha += 0.01; else alpha -= 0.01;
	if(alpha > 1.0){alpha = 1.0;isForward = false;}
	if(alpha < 0.0){alpha = 0.0;isForward = true;}

	// Prepare scene once
	if(timer == NULL){
		Eigen::AlignedBox3d largeBox;
		for(auto pmesh : pwidget->pmeshes) largeBox.extend(pmesh->bbox());

		drawArea()->setSceneRadius( largeBox.sizes().norm() * 2 );
		drawArea()->showEntireScene();

		timer = new QTimer;
		connect(timer, &QTimer::timeout, [=]() { drawArea()->update(); });
		timer->start(20);
	}

	//for(auto pmesh : pwidget->pmeshes) pmesh->drawParticles();
	for(auto pmesh : pwidget->pmeshes) pmesh->drawDebug( *drawArea() );

	// Collect points
	std::vector<Eigen::Vector3f> mixedPoints;

	ParticleMesh * imesh = pwidget->pmeshes.front();
	ParticleMesh * jmesh = pwidget->pmeshes.back();

	for(auto & particle : imesh->particles){
		//mixedPoints.push_back( imesh->realPos( AlphaBlend(alpha, particle.relativePos, 
		//	jmesh->particles[particle.correspondence].relativePos) ).cast<float>() );

		//mixedPoints.push_back(AlphaBlend(alpha, particle.pos, jmesh->particles[particle.correspondence].pos).cast<float>());
	}

	/*if( !pwidget->ui->isOneSided->isChecked() )
	{
		for(auto & particle : jmesh->particles){
			//mixedPoints.push_back( imesh->realPos( AlphaBlend((1.0 - alpha), particle.relativePos, 
			//	imesh->particles[particle.correspondence].relativePos) ).cast<float>() );

			mixedPoints.push_back(AlphaBlend((1.0 - alpha), particle.pos, imesh->particles[particle.correspondence].pos).cast<float>());
		}
	}*/

	// Test meshing
	/*if( false )
	{
		SurfaceMeshModel * m = imesh->meshPoints( mixedPoints );
		m->update_face_normals();
		glBegin(GL_QUADS);
		for(auto f : m->faces()){
			glNormal3dv(m->face_normals()[f].data());
			for(auto v : m->vertices(f))
				glVertex3dv( m->vertex_coordinates()[v].data() );
		}
		glEnd();
		delete m;
		return;
	}*/

	// Setup shader
	static QMap<QString, int> location;
	static QOpenGLShaderProgram * program = NULL;
	if( !program )
	{
		program = new QOpenGLShaderProgram;
		program->addShaderFromSourceCode(QOpenGLShader::Vertex, 
			"attribute highp vec4 vertex;\n"
			"uniform highp mat4 matrix;\n"
			"void main(void)\n"
			"{\n"
			"   gl_Position = vertex * matrix;\n"
			"}");
		program->addShaderFromSourceCode(QOpenGLShader::Fragment, 
			"#version 330\n"
			"out vec4 vFragColor;\n"
			"uniform mediump vec4 color;\n"
			"uniform mediump vec3 lightDir;\n"
			"void main(void)\n"
			"{\n" 
			"vec3 N;\n"
			"N.xy = gl_PointCoord* 2.0 - vec2(1.0);    \n"
			"float mag = dot(N.xy, N.xy);\n"
			"if (mag > 1.0) discard;   // kill pixels outside circle\n"
			"N.z = sqrt(1.0-mag);\n"
			"float diffuse = max(0.0, dot(lightDir, N));\n"
			"vFragColor = color * diffuse;\n"
			"}"
			);
		program->link();
		program->bind();

		location["vertex"] = program->attributeLocation("vertex");
		location["matrix"] = program->uniformLocation("matrix");
		location["color"] = program->uniformLocation("color");
		location["lightDir"] = program->uniformLocation("lightDir");
	}
	else
	{
		GLdouble m[16];
		drawArea()->camera()->getModelViewProjectionMatrix(m);
		QMatrix4x4 pmvMatrix(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7],m[8], m[9], m[10], m[11], m[12], m[13], m[14], m[15]);

		QColor color1(255,0,0,255);
		QColor color2(255,128,0,255);

		QColor color = QColor::fromHsl(AlphaBlend(alpha,color1.hue(),color2.hue()),color1.saturation(),color1.lightness(), color1.alpha());

		QVector3D lightdir(0,0,1);

		program->bind();
		program->setAttributeArray(location["vertex"], mixedPoints.front().data(), 3);
		program->enableAttributeArray(location["vertex"]);
		program->setUniformValue(location["matrix"], pmvMatrix);
		program->setUniformValue(location["color"], color);
		program->setUniformValue(location["lightDir"], lightdir);
	}

	glPointSize(10);
	glEnable(GL_POINT_SPRITE);

	glDrawArrays(GL_POINTS, 0, (int)mixedPoints.size());

	if( program )
	{
		program->disableAttributeArray( location["vertex"] );
		program->release();
	}
}

bool particles::keyPressEvent(QKeyEvent*e)
{
	if(e->key() == Qt::Key_E)
	{
		ParticlesWidget * pwidget = (ParticlesWidget*) widget;
		if(!pwidget || !pwidget->isReady || pwidget->pmeshes.size() < 1) return false;

		auto & pmesh = pwidget->pmeshes.front();

		qglviewer::Vec cen = drawArea()->camera()->revolveAroundPoint();
		size_t pi = pmesh->closestParticles(Vector3 (cen[0],cen[1],cen[2])).front().second;

		starlab::PointSoup * ps = new starlab::PointSoup;

		std::vector<double> diff(pmesh->particles.size());
		Bounds<float> b;

		Spherelib::Sphere sphere( pwidget->ui->sphereResolution->value() );
		auto rotatedIndices = sphere.rotated(6);

		for(auto & p : pmesh->particles)
		{
			size_t pj = p.id;

			auto di = Eigen::Map<Eigen::VectorXf>(&pmesh->desc[pi][0], pmesh->desc[pi].size());
			double d = DBL_MAX;

			for( auto indices : rotatedIndices )
			{
				std::vector<float> rotated_desc;
				for(auto idx : indices) rotated_desc.push_back(pmesh->desc[pj][idx]);

				auto dj = Eigen::Map<Eigen::VectorXf>(&rotated_desc[0], rotated_desc.size());
				auto dist = (di-dj).lpNorm<1>();

				if(dist < d)
					d = dist;
			}

			diff[p.id] = d;

			b.extend(d);
		}

		diff = b.normalize(diff);

		for(auto & p : pmesh->particles)
		{
			ps->addPoint(p.pos, starlab::qtJetColor(diff[p.id]));
		}

		drawArea()->clear();
		drawArea()->addRenderObject(ps);
		
		return true;
	}

	if(e->key() == Qt::Key_R)
	{
		ParticlesWidget * pw = (ParticlesWidget *) widget;

		int steps = 40;
		double theta = 360 / steps;
		for(int i = 0; i < steps; i++)
		{
			pw->ui->rotateAngle->setValue( theta * i );
			reVoxelize();
			processShapes();
			qApp->processEvents();
			drawArea()->update();
			drawArea()->grabFrameBuffer().save(QString("rotate_%1.png").arg(i));
		}
		return true;
	}

	if(e->key() == Qt::Key_Space)
	{
		spheres.clear();

		ParticlesWidget * pwidget = (ParticlesWidget*) widget;
		if(!pwidget || !pwidget->isReady || pwidget->pmeshes.size() < 1) return false;

		auto & pmesh = pwidget->pmeshes.front();

		pmesh->distort();
	}

	if(e->key() == Qt::Key_S)
	{
		ParticlesWidget * pwidget = (ParticlesWidget*) widget;
		if(!pwidget || !pwidget->isReady || pwidget->pmeshes.size() < 1) return false;

		auto & pmesh = pwidget->pmeshes.front();

		qglviewer::Vec cen = drawArea()->camera()->revolveAroundPoint();
		Vector3 q(cen[0],cen[1],cen[2]);

		size_t pi = pmesh->closestParticles(q).front().second;

		Spherelib::Sphere * sphere = new Spherelib::Sphere( pwidget->ui->sphereResolution->value(), pmesh->particles[pi].pos, 0.01 );

		// Reconstruct when using rotation invariant
		if( pwidget->ui->bands->value() > 0 )
		{
			// Uniform sampling on the sphere
			Spherelib::Sphere ss( pwidget->ui->sphereResolution->value() );
			std::vector< Eigen::Vector3d > sampledRayDirections  = ss.rays();
			std::vector< size_t > antiRays = ss.antiRays();

			// Rotation invariant descriptor
			SphericalHarmonic<Vector3,float> sh( std::max(1,pwidget->ui->bands->value()) );
			std::vector< SHSample<Vector3,float> > sh_samples;
			sh.SH_setup_spherical( sampledRayDirections, sh_samples );
			sh.SH_project_function( pmesh->desc[ pi ], sh_samples, pmesh->desc[pi] );

			sphere->setValues( sh.SH_reconstruct(sampledRayDirections, pmesh->desc[pi]) );
		}
		else
			sphere->setValues( pmesh->desc[ pi ] );

		sphere->normalizeValues();

		mainWindow()->setStatusBarMessage( QString("Particle [%1] with maximum [%2] and minimum [%3] measure [%4] flat [%5]").arg(pi).arg(*std::max_element(
			pmesh->desc[ pi ].begin(),pmesh->desc[ pi ].end())).arg(*std::min_element(pmesh->desc[ pi ].begin(),
			pmesh->desc[ pi ].end())).arg(pmesh->particles[pi].measure).arg(pmesh->particles[pi].flat) );

		spheres.clear();
		spheres.push_back( sphere );
	}

	drawArea()->update();
	return true;
}
