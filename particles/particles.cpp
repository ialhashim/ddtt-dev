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
#include "Bounds.h"

#include "Raytracing.h"

#include "kmeans.h"

#include "spherelib.h"
#include "SphericalHarmonic.h"
std::vector<Spherelib::Sphere*> spheres;

#include "BasicTable.h"

#ifdef WIN32
namespace std{ int isfinite(double x) {return _finite(x);} }
#endif

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

						if(!hit.isHit) hit.distance = s->grid.unitlength;

						descriptor[pi][r++] = hit.distance;
						rayCount++;
					}

					p.avgDiameter = 0;
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
								double scale = 1.2;
								if( max_rad > ma_point_rad[pi] * scale )
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

								s->particles[pi].avgDiameter = s->particles[pj].avgDiameter;
								//s->particles[pi].measure = s->particles[pj].measure;
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

					// Normalize diameters
					if( true )
					{
						Bounds<float> b;
						for(auto & p : s->particles) b.extend(p.avgDiameter);
						for(auto & p : s->particles) p.avgDiameter = b.normalized(p.avgDiameter);
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

			// Rotation invariant descriptor
			if( pw->ui->bands->value() > 0 )
			{
				s->sig = std::vector< std::vector<float> >(s->particles.size(), std::vector<float>( pw->ui->bands->value() ));

				#pragma omp parallel for
				for(int pi = 0; pi < (int)s->particles.size(); pi++)
				{
					auto & p = s->particles[pi];

					auto desc = s->desc[pi];
					// Normalize descriptor
					//desc = Bounds<float>::from( desc ).normalize( desc );

					std::vector<float> coeff;
					sh.SH_project_function(desc, sh_samples, coeff);

					s->sig[p.id] = sh.SH_signature(coeff);
				}
			}

			// Axis perpendicular to floor flow
			int numIterations = pw->ui->experimentIter->value();
			if( numIterations > 0 )
			{
				#pragma omp parallel for
				for(int pi = 0; pi < (int)s->particles.size(); pi++)
				{
					QMap<double, size_t> fnei;
					for(auto pj : s->neighbourhood( s->particles[pi], 1 ))
						fnei[ s->particles[pj].measure ] = pj;

					auto fromParticle = fnei[fnei.firstKey()];
					auto toParticle = fnei[fnei.lastKey()];

					Vector3 a = s->particles[ fromParticle ].pos;
					Vector3 b = s->particles[ toParticle ].pos;
					Vector3 axis = ( b - a ).normalized();

					if(fromParticle == toParticle)
						axis = Vector3(0,0,1);
					
					s->particles[pi].axis = axis;
				}

				// Average direction
				for(int itr = 0; itr < numIterations; itr++)
				{
					std::vector<Vector3> smoothedDirections(s->particles.size());

					#pragma omp parallel for
					for(int pi = 0; pi < (int)s->particles.size(); pi++)
					{
						auto nei = s->neighbourhood( s->particles[pi], 1 );
						Vector3 sumAxis(0,0,0);
						for(auto pj : nei) sumAxis += s->particles[pj].axis;
						sumAxis /= nei.size();

						smoothedDirections[pi] = sumAxis.normalized();
					}

					for(auto & p : s->particles)
						s->particles[p.id].axis = smoothedDirections[p.id];
				}

				// Average
				{
					#pragma omp parallel for
					for(int pi = 0; pi < (int)s->particles.size(); pi++)
					{
						auto nei = s->neighbourhood( s->particles[pi], 1 );
						Vector3 sumAxis(0,0,0);
						for(auto pj : nei) sumAxis += s->particles[pj].axis;
						sumAxis /= nei.size();

						s->particles[pi].avgDiameter = sumAxis.norm();
					}
				}

				{
					Bounds<float> b;
					for(auto & p : s->particles) b.extend(p.avgDiameter);
					for(auto & p : s->particles) p.avgDiameter = b.normalized(p.avgDiameter);
				}

				// [DEBUG] visualize axis experiment
				if( pw->ui->showExperiment->isChecked() )
				{
					drawArea()->clear();
					starlab::LineSegments * vs = new starlab::LineSegments(2);
					for(auto & particle : s->particles){
						Vector3 d = particle.axis;
						//double value = abs(dot(d, Vector3(0,0,1)));
						double value = particle.avgDiameter;
						vs->addLine(particle.pos,  Vector3(particle.pos + d * 0.01), starlab::qtJetColor(value));
					}
					drawArea()->addRenderObject(vs);

					return;
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
		}
	}

	// Special features
	if( true )
	{
		// Different descriptor options
		for(auto & s : pw->pmeshes)
		{
			#pragma omp parallel for
			for(int pi = 0; pi < (int)s->particles.size(); pi++)
			{
				std::vector<float> new_desc;

				auto & p = s->particles[pi];

				if(pw->ui->useDescriptor->isChecked() ) 	new_desc = s->desc[p.id];
				if(pw->ui->useRotationInv->isChecked())		new_desc = s->sig[p.id];
				if(pw->ui->useDiameter->isChecked()   ) 	new_desc.push_back(p.avgDiameter);
				if(pw->ui->useGroundDist->isChecked() )		new_desc.push_back(p.measure);
				if(pw->ui->useHeight->isChecked()     )		new_desc.push_back(p.pos.z());

				if(pw->ui->experimentIter->value())
				{
					double vertical = abs(dot(p.axis, Vector3(0,0,1)));
					new_desc.push_back(vertical);
				}

				if(new_desc.empty()) new_desc.push_back(p.pos.z()); // simply height..

				s->desc[p.id] = new_desc;
			}

			showTable(s->desc, std::min(size_t(1000),s->particles.size()));
		}
	}

	// k-means clustering
	if( true )
	{
		int K = pw->ui->kclusters->value();
		int numIterations = 300;
		double minchangesfraction = 0.005;
		typedef clustering::l2norm_squared< std::vector<float> > dist_fn;

		/*bool isClusterShapesTogether = false;
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
		else*/
		{
			for(auto & s : pw->pmeshes)
			{	
				if(!s->particles.size()) continue;

				// Cluster
				clustering::kmeans< std::vector< std::vector<float> >, dist_fn > km(s->desc, K);

				// Special seeding when using ground distance
				if( pw->ui->useGroundDistSeed->isChecked() )
				{
					std::set<size_t> seeds;
					for(int i = 0; i < K; i++)
					{
						double t = double(i) / (K-1);
						int idx = t * (s->pathFromFloor.size()-1);
						seeds.insert( s->pathFromFloor[idx] );
					}

					// Add seeds
					km._centers.clear();
					for(auto pid : seeds) 
					{
						km._centers.push_back( s->desc[pid] );

						// DEBUG
						if( pw->ui->showSeeds->isChecked() )
							drawArea()->drawPoint(s->particles[pid].pos, 20, Qt::black);
					}
				}

				km.run(numIterations, minchangesfraction);

				// Assign classes
				#pragma omp parallel for
				for(int pi = 0; pi < (int)s->particles.size(); pi++){
					auto & p = s->particles[pi];
					p.segment = (int) km.clusters()[p.id];
				}
			}
		}
	}

	// Merge smaller clusters with larger ones
	for(int i = 0; i < pw->ui->simplifyCluster->value(); i++)
	{
		for(auto & s : pw->pmeshes)
		{	
			std::vector<int> newSegmentAssignment(s->particles.size());

			for(auto & p : s->particles)
			{
				std::map<size_t,size_t> clusterHistogram;
				for( auto pj : s->neighbourhood(p, 2) )
					clusterHistogram[ s->particles[pj].segment ]++;

				std::vector< std::pair<size_t,size_t> > ch( clusterHistogram.begin(), clusterHistogram.end() );
				std::sort(ch.begin(),ch.end(),[](std::pair<size_t,size_t> a, std::pair<size_t,size_t> b){ return a.second < b.second; });

				// majority rule
				newSegmentAssignment[p.id] = (int)ch.back().first;
				//p.segment = ch.back().first;
			}

			for(auto & p : s->particles) p.segment = newSegmentAssignment[p.id];
		}
	}

	// Find symmetric parts
	if( pw->ui->showSymmetry->isChecked() )
	{	
		struct Plane{ 
			Vector3 pos,n; double weight; 
			Plane(Vector3 pos, Vector3 n, double weight):pos(pos),n(n),weight(weight){}

			static std::vector<Plane> mergePlanes(const std::vector<Plane> & groupPlanes, double similiairy_threshold)
			{
				// Compare planes
				std::vector<Plane> mergedPlanes;
				QMap<size_t,int> count;
				QSet<size_t> used;
				int k = 0;

				for(size_t i = 0; i < groupPlanes.size(); i++){
					if(used.contains(i)) continue;

					auto plane_i = groupPlanes[i];
					mergedPlanes.push_back( plane_i );

					count[k]++;

					for(size_t j = i+1; j < groupPlanes.size(); j++){
						if(used.contains(j)) continue;
						auto plane_j = groupPlanes[j];

						double dot = plane_i.n.dot(plane_j.n);
						double similiairy = abs(dot);
						if( similiairy > similiairy_threshold )
						{
							used.insert(j); // mark as used
							mergedPlanes[k].pos += plane_j.pos;
							mergedPlanes[k].n += (dot > 0) ? plane_j.n : -plane_j.n;
							mergedPlanes[k].weight = std::max(mergedPlanes[k].weight, plane_j.weight);
							count[k]++;
						}
					}

					k++;
				}

				for(size_t i = 0; i < mergedPlanes.size(); i++){
					auto & plane = mergedPlanes[i];
					plane.pos /= count[i];
					plane.n /= count[i];

					plane.n.normalize();
				}

				return mergedPlanes;
			}
		};

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

				// Segment properties
				{
					Vector3 sum(0,0,0);
					Eigen::AlignedBox3d bbox;
					Bounds<double> interval;
					
					for(auto v : c_seg->vertices) 
					{
						Vector3 p = s->particles[v].pos;
						sum += p;
						bbox.extend( p );
						interval.extend( s->particles[v].measure );
					}

					Vector3 centroid = sum / c_seg->vertices.size();
					c_seg->setProperty< Vector3 >("centroid", centroid);
					c_seg->setProperty< Eigen::AlignedBox3d >("bbox", bbox);
					c_seg->setProperty< Bounds<double> >("interval", interval);
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

				if(s1 != s2) neiGraph.AddEdge( (uint)s1->uid, (uint)s2->uid, 1.0 );
			}

			//double threshold = s->grid.unitlength * 3; // two voxels

			int pcount = 0;

			std::vector<Plane> allPlanes;

			for(auto & segGroup : groupedSegments)
			{
				auto & group = segGroup.second;

				std::vector<Plane> groupPlanes;
				
				// Compare pair of parts in the same segment
				for(size_t i = 0; i < group.size(); i++){
					for(size_t j = i+1; j < group.size(); j++){
						Vector3 centroid_i = *group[i]->getProperty<Vector3>("centroid");
						Vector3 centroid_j = *group[j]->getProperty<Vector3>("centroid");
						//if( std::abs(centroid_i.z() - centroid_j.z()) > threshold ) continue;
						Vector3 planeCenter = (centroid_i+centroid_j) * 0.5;
						Vector3 planeNormal = (centroid_i-centroid_j).normalized();

						Eigen::AlignedBox3d bbox_i = *group[i]->getProperty<Eigen::AlignedBox3d>("bbox");
						Eigen::AlignedBox3d bbox_j = *group[j]->getProperty<Eigen::AlignedBox3d>("bbox");
						double bbox_diff = (bbox_i.sizes() - bbox_j.sizes()).norm();

						Bounds<double> interval_i = *group[i]->getProperty< Bounds<double> >("interval");
						Bounds<double> interval_j = *group[j]->getProperty< Bounds<double> >("interval");
						
						// by volume
						double minVol = std::min(interval_i.count, interval_j.count);
						double maxVol = std::max(interval_i.count, interval_j.count);
						double ratio = minVol / maxVol;
						if(ratio < 0.4) continue;

						groupPlanes.push_back( Plane(planeCenter, planeNormal, ratio) );
						
						//ps->addPlane( Vector3(planeCenter, planeNormal);
						//drawArea()->drawSegment(centroid_i,centroid_j);
						pcount++;
					}
				}

				auto mergedPlanes = Plane::mergePlanes(groupPlanes, 0.9);
				for (auto & plane : mergedPlanes) allPlanes.push_back(plane);
			}

			auto rndColors = rndColors2(pcount);

			auto mergedPlanes = Plane::mergePlanes(allPlanes, 0.9);

			for(size_t i = 0; i < mergedPlanes.size(); i++)
			{
				auto & plane = mergedPlanes[i];
				auto color = QColor::fromRgbF(rndColors[i].redF(),rndColors[i].greenF(),rndColors[i].blueF());
				starlab::PlaneSoup * ps = new starlab::PlaneSoup(0.1 * plane.weight, true, color);
				ps->addPlane( plane.pos, plane.n );

				drawArea()->addRenderObject(ps);
			}

			mainWindow()->setStatusBarMessage(QString("segment count (%1) planes count (%2)").arg( allSegs.size() ).arg( pcount ));
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
