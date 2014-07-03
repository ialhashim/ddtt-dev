#include "particles.h"
#include "particles-widget.h"
#include "ui_particles-widget.h"
#include <omp.h>
#include <numeric>
#include <cmath>

#include "SurfaceMeshModel.h"
#include "SurfaceMeshHelper.h"
using namespace SurfaceMesh;

#include <QOpenGLShaderProgram>
#include <QWidget>
#include <QScrollArea>
#include <QScrollBar>
#include <QFileDialog>
#include <QElapsedTimer>

#define AlphaBlend(alpha, start, end) ( ((1-alpha) * start) + (alpha * end) )
QTimer * timer = NULL;

QStringList files;

#include "ParticleMesh.h"
#include "myglobals.h"

#include "Raytracing.h"

#include "kmeans.h"

#include "spherelib.h"
#include "SphericalHarmonic.h"
std::vector<Spherelib::Sphere*> spheres;

#ifdef WIN32
namespace std{ int isfinite(double x) {return _finite(x);} }
#endif

void particles::processShapes()
{
	if(!widget) return;

	ParticlesWidget * pw = (ParticlesWidget *) widget;
	pw->isReady = false;

	if(pw->pmeshes.empty()) return;

	drawArea()->clear();

	// Check grid size changes, redo voxelization
	if( pw->pmeshes.front()->grid.gridsize != pw->ui->gridsize->value() ){
		pw->pmeshes.clear();
		for(auto filename : files){
			SurfaceMeshModel fromMesh( filename, QFileInfo(filename).baseName() );
			fromMesh.read( filename.toStdString() );
			pw->pmeshes.push_back( new ParticleMesh( &fromMesh, pw->ui->gridsize->value() ) );
		}
	}

	// Show voxelized meshes
	//for(auto & s : pw->pmeshes) document()->addModel(s->surface_mesh->clone());

	QElapsedTimer allTimer; allTimer.start();

	mainWindow()->setStatusBarMessage("Shapes loaded, now processing..");
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
			s->debug.clear();

			QElapsedTimer timer; timer.start();

			// Hit results are saved as vector of distances
			std::vector< std::vector<float> > & descriptor = (s->desc = std::vector< std::vector<float> >( 
				s->particles.size(), std::vector<float>( perSampleRaysCount ) ));

			int rayCount = 0;

			// Shooting rays from inside the volume
			{
				// Medial points record
				std::vector<Eigen::Vector3f> ma_point(s->particles.size(), Eigen::Vector3f(0,0,0));
				std::vector<bool> ma_point_active(s->particles.size(), false);
				std::vector<float> ma_point_rad(s->particles.size(), 0);

				// Accelerated raytracing
				raytracing::Raytracing<Eigen::Vector3d> rt( s->surface_mesh );

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
							double maxUsedRadius = -DBL_MAX;
							std::vector<bool> isVisisted( sampledRayDirections.size(), false );

							double sumDiameter = 0;
							int visitCount = 0;

							double diameter, radius, out_dist, d2;
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

								surface_kdtree.tree->knnSearch(&midpoint[0], 1, &ret_index, &out_dist);
								d2 = out_dist;

								isFullyInside = d2 > (radius * radius);

								if( isFullyInside && radius > maxUsedRadius )
								{
									maxUsedRadius = radius;
									maPoint = midpoint;
								}
							}

							// Found a medial point
							if( maxUsedRadius > 0 )
							{
								ma_point[pi] = maPoint.cast<float>();
								ma_point_rad[pi] = maxUsedRadius;
								ma_point_active[pi] = true;
							}

							// Record average diameter around particle
							p.avgDiameter = sumDiameter / visitCount;
						}

						// Normalize average diameter
						if( true )
						{
							double minDiameter = DBL_MAX, maxDiameter = -DBL_MAX;
							for(auto & p : s->particles) { 
								minDiameter = std::min(minDiameter, p.avgDiameter); 
								maxDiameter = std::max(maxDiameter, p.avgDiameter);
							}

							for(auto & p : s->particles) 
								p.avgDiameter = (p.avgDiameter-minDiameter) / (maxDiameter-minDiameter);
						}
					}

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
								for(auto pj : s->neighbourhood( s->particles[pi], 2 ))
									max_rad = std::max(max_rad, ma_point_rad[pj]);

								// A neighbor has wider access than me => I'm not so medial
								if( max_rad > ma_point_rad[pi] )
									ma_point_active[pi] = false;
							}
						}
					}

					// Contract medial points to their neighborhood center
					for( int i = 0; i < pw->ui->contractMedial->value(); i++ )
					{
						NanoKdTree kdtree;
						for(int pi = 0; pi < (int)s->particles.size(); pi++)
							if( ma_point_active[pi] ) kdtree.addPoint( ma_point[pi].cast<double>() );
						kdtree.build();

						// Get a set of neighbors for each medial point
						std::vector< std::vector<Vector3> > nei(s->particles.size());
						for(int pi = 0; pi < (int)s->particles.size(); pi++){
							if( !ma_point_active[pi] ) continue;

							KDResults matches;
							kdtree.ball_search( ma_point[pi].cast<double>(), s->grid.unitlength * 1.5, matches );
							for(auto m : matches) 
								if((kdtree.cloud.pts[m.first]-ma_point[pi].cast<double>()).norm() > 1e-6) 
									nei[pi].push_back(kdtree.cloud.pts[m.first]);
						}

						// Move to center of neighborhood
						for(int pi = 0; pi < (int)s->particles.size(); pi++){
							if( !ma_point_active[pi] ) continue;
							Vector3 sumv(0,0,0);
							for(auto v : nei[pi]) sumv += v;
							ma_point[pi] = sumv.cast<float>() / nei[pi].size();
						}
					}

					// Report
					mainWindow()->setStatusBarMessage( QString("Medial particles (%1 ms)").arg( timer.elapsed() ) );
					timer.restart();

					// [DEBUG] medial points
					if( pw->ui->showMedial->isChecked() )
					{
						starlab::PointSoup * ps = new starlab::PointSoup;
						float maxVal = -1;
						for(size_t i = 0; i < ma_point_active.size(); i++) maxVal = std::max(maxVal, ma_point_rad[i]);
						for(size_t i = 0; i < ma_point_active.size(); i++) {
							if(ma_point_active[i])
								ps->addPoint(Vector3(ma_point[i].cast<double>()), starlab::qtJetColor(ma_point_rad[i], 0, maxVal));
						}
						drawArea()->deleteAllRenderObjects();
						drawArea()->addRenderObject(ps);
					}

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
									raytracing::RayHit hit = rt.hit( maPoint, (d + Vector3(1,1,1) * 1e-6).normalized() ); // weird bug..
									descriptor[pi][r++] = hit.distance;
								}
							}
						}

						// Assign descriptor of all particles to their nearest medial point
						NanoKdTree medial_kdtree;
						std::map<size_t,size_t> ma_particle;
						for(size_t i = 0; i < ma_point_active.size(); i++){
							if( ma_point_active[i] ){
								ma_particle[ma_particle.size()] = i;
								medial_kdtree.addPoint( ma_point[i].cast<double>() );
							}
						}
						medial_kdtree.build();

						#pragma omp parallel for
						for(int pi = 0; pi < (int)s->particles.size(); pi++)
						{
							if( !ma_point_active[pi] )
							{
								size_t pj = ma_particle[medial_kdtree.closest( s->particles[pi].pos )];
								descriptor[pi] = descriptor[pj];
							}

							auto maxelement = std::max_element(descriptor[pi].begin(),descriptor[pi].end());
							s->particles[pi].direction = sampledRayDirections[ maxelement - descriptor[pi].begin() ];
						}
					}

					// Report
					mainWindow()->setStatusBarMessage( QString("Projection to Medial particles (%1 ms)").arg( timer.elapsed() ) );
					timer.restart();

					// [DEBUG] show projected skeleton
					if( pw->ui->projectSkeleton->isChecked() )
					{
						NanoKdTree kdtree;
						for(int pi = 0; pi < (int)s->particles.size(); pi++)
							if( ma_point_active[pi] ) kdtree.addPoint( ma_point[pi].cast<double>() );
						kdtree.build();

						SurfaceMeshModel * m = s->surface_mesh->clone();
						auto points = m->vertex_coordinates();
						for(auto v : m->vertices())
							points[v] = kdtree.cloud.pts[ kdtree.closest(points[v]) ];

						document()->addModel( m );
						drawArea()->setRenderer(m, "Flat Wire");
					}

					// [DEBUG] show main 'direction' of particles
					if( false )
					{
						starlab::LineSegments * vs = new starlab::LineSegments(2);
						for(auto & particle : s->particles){
							Eigen::Vector4d color(abs(particle.direction[0]), abs(particle.direction[1]), abs(particle.direction[2]), particle.alpha);
							color[0] *= color[0];color[1] *= color[1];color[2] *= color[2];
							vs->addLine(particle.pos,  Vector3(particle.pos + particle.direction * 0.004), QColor::fromRgbF(color[0],color[1],color[2],1));
						}
						s->debug.push_back(vs);
					}
				}
			}

			// [DEBUG] visualize distance to ground
			if( pw->ui->showDistGround->isChecked() )
			{
				for(auto & particle : s->particles)
					drawArea()->drawPoint(particle.pos, 5, starlab::qtJetColor(particle.measure));
				return;
			}

			// Normalize descriptor
			if( pw->ui->normalizeDesc->isChecked() )
			{
				float mn = FLT_MAX, mx = -FLT_MIN;
				for(auto & p : s->particles) { 
					auto & desc = s->desc[p.id];
					mn = std::min(mn, *std::min_element(desc.begin(),desc.end())); 
					mx = std::max(mx, *std::max_element(desc.begin(),desc.end()));
				}

				for(auto & p : s->particles){
					auto & desc = s->desc[p.id];
					for(auto & d : desc) d = (d-mn) / (mx-mn);
				}
			}

			// Normalize descriptor per particle
			if( pw->ui->normalizeSphereFn->isChecked() )
			{
				#pragma omp parallel for
				for(int pi = 0; pi < (int)s->particles.size(); pi++)
				{
					auto & p = s->particles[pi];
					std::vector<float> & desc = descriptor[p.id];
					double max_desc = *std::max_element(desc.begin(),desc.end());
					double min_desc = *std::min_element(desc.begin(),desc.end());
					for(auto & d : desc) d = (d-min_desc) / (max_desc-min_desc);
				}
			}

			// Rotation invariant descriptor
			if( pw->ui->bands->value() > 0 )
			{
				s->sig = std::vector< std::vector<float> >(s->particles.size(), std::vector<float>( pw->ui->bands->value() ));

				#pragma omp parallel for
				for(int pi = 0; pi < (int)s->particles.size(); pi++)
				{
					auto & p = s->particles[pi];
					std::vector<float> & desc = descriptor[p.id];

					std::vector<float> coeff;
					sh.SH_project_function(desc, sh_samples, coeff);

					s->sig[p.id] = sh.SH_signature(coeff);

					if( pw->ui->useRotationInv->isChecked() )
						desc = s->sig[p.id];
				}
			}
		}
	}

	// Special features
	if( true )
	{
		// Different descriptor options
		for(auto & s : pw->pmeshes)
		{
			for(auto & p : s->particles)
			{
				std::vector<float> new_desc;

				if(pw->ui->useDescriptor->isChecked()) 
					new_desc = s->desc[p.id];

				if(pw->ui->useGroundDist->isChecked())
					new_desc.push_back(p.measure);

				if(pw->ui->useDiameter->isChecked()) 
					new_desc.push_back(p.avgDiameter);

				s->desc[p.id] = new_desc;
			}
		}
	}

	// k-means clustering
	if( true )
	{
		int K = pw->ui->kclusters->value();
		int numIterations = 300;
		double minchangesfraction = 0.005;
		typedef clustering::l2norm_squared< std::vector<float> > dist_fn;

		bool isClusterShapesTogether = false;

		if( isClusterShapesTogether )
		{
			// Collect descriptors
			std::vector< std::vector<float> > allDesc;
			for(auto & s : pw->pmeshes){
				if( pw->ui->bands->value() > 0 )
					allDesc.insert(allDesc.end(), s->sig.begin(), s->sig.end());
				else
					allDesc.insert(allDesc.end(), s->desc.begin(), s->desc.end());
			}

			// Cluster
			clustering::kmeans< std::vector< std::vector<float> >, dist_fn > km(allDesc, K);
			km.run(numIterations, minchangesfraction);

			// Assign classes
			int pi = 0;
			for(auto & s : pw->pmeshes)
				for(auto & p : s->particles)
					p.segment = (int) km.clusters()[pi++];
		}
		else
		{
			for(auto & s : pw->pmeshes)
			{	
				if(!s->particles.size()) continue;

				// Cluster
				clustering::kmeans< std::vector< std::vector<float> >, dist_fn > km(s->desc, K);

				// Special seeding when using ground distance
				if( pw->ui->useGroundDistSeed->isChecked() )
				{
					// Sort seeds based on ground
					typedef std::pair<size_t,double> seedpair;
					std::vector< seedpair > pg;
					for(auto & p : s->particles) pg.push_back( std::make_pair(p.id, p.measure) );
					std::sort(pg.begin(),pg.end(),[](seedpair a, seedpair b){ return a.second < b.second; });

					km._centers.clear();
					for(int i = 0; i < K; i++){
						size_t pid = pg[(double(i)/K) * (pg.size()-1)].first;
						km._centers.push_back( s->desc[pid] );
					}
				}

				km.run(numIterations, minchangesfraction);

				// Assign classes
				for(auto & p : s->particles)
					p.segment = (int) km.clusters()[p.id];
			}
		}
	}

	// Merge smaller clusters with larger ones
	for(int i = 0; i < pw->ui->simplifyCluster->value(); i++)
	{
		for(auto & s : pw->pmeshes)
		{	
			for(auto & p : s->particles)
			{
				std::map<size_t,size_t> clusterHistogram;
				for( auto pj : s->neighbourhood(p, 2) )
					clusterHistogram[ s->particles[pj].segment ]++;

				std::vector< std::pair<size_t,size_t> > ch( clusterHistogram.begin(), clusterHistogram.end() );
				std::sort(ch.begin(),ch.end(),[](std::pair<size_t,size_t> a, std::pair<size_t,size_t> b){ return a.second < b.second; });

				// majority rule
				p.segment = ch.back().first;
			}
		}
	}

	// Find symmetric parts
	if( pw->ui->showSymmetry->isChecked() )
	{		
		for(auto & s : pw->pmeshes)
		{		
			// Extract groups of similar segments
			typedef GenericGraphs::Graph<uint,double> GraphType;
			std::map<int, std::vector<GraphType*> > groupedSegments;
			std::map<size_t, GraphType*> allSegs;
			GraphType rawNeiGraph, neiGraph;

			for(auto seg : s->segmentToComponents( rawNeiGraph ))
			{
				int si = s->particles[seg.FirstVertex()].segment;

				auto c_seg = new GraphType( seg );
				groupedSegments[ si ].push_back( c_seg );
				c_seg->sid = si;

				// Track segment
				allSegs[ c_seg->uid ] = c_seg;

				// Segment centroid
				{
					Vector3 sum(0,0,0);
					for(auto v : c_seg->vertices) sum += s->particles[v].pos;
					Vector3 centroid = sum / c_seg->vertices.size();

					c_seg->setProperty<Vector3>("centroid", centroid);
				}
			}

			// Function to map vertex ID to segment pointer
			auto mapVertToSegment = [=]( uint v ){ 
				GraphType * g = NULL;
				for(auto uid_segment : allSegs) {
					if( uid_segment.second->IsHasVertex(v) ) {
						g = uid_segment.second;
						break;
					}
				}
				return g; 
			};

			// Record connectivity relations between segments
			for(auto e : rawNeiGraph.GetEdgesSet())
			{
				auto s1 = mapVertToSegment( e.index );
				auto s2 = mapVertToSegment( e.target );

				neiGraph.AddEdge( s1->uid, s2->uid, 1.0 );
			}

			double threshold = s->grid.unitlength * 3; // two voxels

			starlab::PlaneSoup * ps = new starlab::PlaneSoup( threshold * 2 );
			int pcount = 0;

			for(auto & segGroup : groupedSegments)
			{
				auto & group = segGroup.second;
				for(size_t i = 0; i < group.size(); i++)
				{
					if(group[i]->vertices.size() < 5) continue;

					for(size_t j = i+1; j < group.size(); j++){
						Vector3 centroid_i = *group[i]->getProperty<Vector3>("centroid");
						Vector3 centroid_j = *group[j]->getProperty<Vector3>("centroid");

						//if( std::abs(centroid_i.z() - centroid_j.z()) > threshold ) continue;

						ps->addPlane( Vector3((centroid_i+centroid_j) * 0.5), (centroid_i-centroid_j).normalized() );
						pcount++;
					}
				}
			}

			mainWindow()->setStatusBarMessage(QString("segment count (%1) planes count (%2)").arg( allSegs.size() ).arg( pcount ));

			if( pw->ui->showSymmetry->isChecked() )
			{
				drawArea()->addRenderObject(ps);
			}
		}
	}

	/// [ Correspondence ] match particles
	if( true && pw->pmeshes.size() > 1 )
	{
		typedef clustering::l2norm_squared< std::vector<float> > dist_fn;

		for(size_t i = 0; i < pw->pmeshes.size(); i++){
			for(size_t j = i+1; j < pw->pmeshes.size(); j++){
				NanoKdTree * itree = pw->pmeshes[i]->relativeKdtree;
				NanoKdTree * jtree = pw->pmeshes[j]->relativeKdtree;

				for(auto & iparticle : pw->pmeshes[i]->particles)
				{
					if( false )
					{
						iparticle.correspondence = jtree->closest( iparticle.relativePos );
					}
					else
					{
						// Experiment
						KDResults matches;
						jtree->k_closest( iparticle.relativePos, 20, matches );

						QMap<double, size_t> measures;
						for(auto p : matches) 
						{
							double weight = p.second;
							double desc_dist = dist_fn()( pw->pmeshes[i]->desc[iparticle.id], pw->pmeshes[j]->desc[p.first] );
							measures[ weight * desc_dist ] = p.first;
						}
						iparticle.correspondence = measures[ measures.keys().front() ];
					}
				}

				for(auto & jparticle : pw->pmeshes[j]->particles)
				{
					if( true )
					{
						jparticle.correspondence = itree->closest( jparticle.relativePos );
					}
					else
					{
						// Experiment
						KDResults matches;
						itree->k_closest( jparticle.relativePos, 20, matches );

						QMap<double, size_t> measures;
						for(auto p : matches) 
						{
							double weight = p.second;
							double desc_dist = dist_fn()( pw->pmeshes[j]->desc[jparticle.id], pw->pmeshes[i]->desc[p.first] );
							measures[ weight * desc_dist ] = p.first;
						}
						jparticle.correspondence = measures[ measures.keys().front() ];
					}
				}
			}
		}
	}

	// Report
	mainWindow()->setStatusBarMessage( QString("All time (%1 ms)").arg( allTimer.elapsed() ) );

	pw->isReady = true;

	drawArea()->update();

	emit( shapesProcessed() );
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
	connect(pw->ui->testButton, &QPushButton::released, [=]{
		//for(auto p : sphere_fibonacci_points( 100 ))drawArea()->drawPoint(p, 5);
		//drawArea()->update();
	});

	// Load and process shapes:
	connect(pw->ui->loadShapes, &QPushButton::released, [=]{
		files = QFileDialog::getOpenFileNames(nullptr, "Open Shapes", "", "All Supported (*.obj *.off)");

		QElapsedTimer timer; timer.start();

		for(auto filename : files)
		{
			SurfaceMeshModel fromMesh( filename, QFileInfo(filename).baseName() );
			fromMesh.read( filename.toStdString() );

			pw->pmeshes.push_back( new ParticleMesh( &fromMesh, pw->ui->gridsize->value() ) );
		}

		mainWindow()->setStatusBarMessage(QString("Shapes loaded and voxelized (%1 ms)").arg(timer.elapsed()));

		emit( pw->shapesLoaded() );
	});

	// Post-processing
	connect(pw, SIGNAL(shapesLoaded()), SLOT(processShapes()));
	connect(pw->ui->processShapesButton, SIGNAL(clicked()), SLOT(processShapes()));

#ifdef QT_DEBUG
	pw->ui->gridsize->setValue(32);
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

	if(pwidget->pmeshes.size() < 2) 
	{
		glPushMatrix();

		// Evaluation
		for( auto s : pwidget->pmeshes )
		{
			s->drawParticles( drawArea()->camera() );
			s->drawDebug( *drawArea() );

			glTranslated(s->bbox().sizes().x() * 1.1, 0, 0);
		}

		glPopMatrix();

		return;
	}

	// Experimental
	static bool isForward = true;
	static double alpha = 0;
	if(isForward) alpha += 0.005; else alpha -= 0.005;
	if(alpha > 1.0){alpha = 1.0;isForward = false;}
	if(alpha < 0.0){alpha = 0.0;isForward = true;}

	// Prepare scene once
	Eigen::AlignedBox3d largeBox;
	for(auto pmesh : pwidget->pmeshes) largeBox.extend(pmesh->bbox());
	if(timer == NULL){
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
	for(size_t i = 0; i < pwidget->pmeshes.size(); i++)
	{
		for(size_t j = i+1; j < pwidget->pmeshes.size(); j++)
		{
			ParticleMesh * imesh = pwidget->pmeshes[i];
			ParticleMesh * jmesh = pwidget->pmeshes[j];

			for(auto particle : imesh->particles){
				Vector3 p = AlphaBlend(alpha, particle.pos, jmesh->particles[particle.correspondence].pos);
				mixedPoints.push_back( p.cast<float>() );
			}

			for(auto particle : jmesh->particles){
				Vector3 p = AlphaBlend((1.0-alpha), particle.pos, imesh->particles[particle.correspondence].pos);
				mixedPoints.push_back( p.cast<float>() );
			}
		}
	}

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

		QColor color(255,0,0,255);
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
	if(e->key() == Qt::Key_Space)
	{
		spheres.clear();
	}

	if(e->key() == Qt::Key_S)
	{
		ParticlesWidget * pwidget = (ParticlesWidget*) widget;
		if(!pwidget || !pwidget->isReady || pwidget->pmeshes.size() < 1) return false;

		auto & pmesh = pwidget->pmeshes.front();

		qglviewer::Vec cen = drawArea()->camera()->revolveAroundPoint();
		Vector3 q(cen[0],cen[1],cen[2]);

		size_t pi = 0;

		double minDist = DBL_MAX;
		for(auto & p : pmesh->particles){
			double dist = (p.pos-q).norm();
			if(dist < minDist){
				minDist = dist;
				pi = p.id;
			}
		}

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

		mainWindow()->setStatusBarMessage( QString("Particle [%1] with maximum [%2] and minimum [%3]").arg(pi).arg(*std::max_element(
			pmesh->desc[ pi ].begin(),pmesh->desc[ pi ].end())).arg(*std::min_element(pmesh->desc[ pi ].begin(),pmesh->desc[ pi ].end())) );

		spheres.clear();
		spheres.push_back( sphere );
	}

	drawArea()->update();
	return true;
}
