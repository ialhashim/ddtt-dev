#include "particles.h"
#include "particles-widget.h"
#include "ui_particles-widget.h"
#include "interfaces/ModePluginDockWidget.h"
#include <omp.h>
#include <numeric>

#include "SurfaceMeshModel.h"
#include "SurfaceMeshHelper.h"
using namespace SurfaceMesh;

#include <QOpenGLShaderProgram>
#include <QWidget>
#include <QFileDialog>

#define AlphaBlend(alpha, start, end) ( ((1-alpha) * start) + (alpha * end) )
QTimer * timer = NULL;

#include "ParticleMesh.h"
#include "myglobals.h"

#include "Raytracing.h"

#include "kmeans.h"

void particles::create()
{
	if( widget ) return;

	//drawArea()->setAxisIsDrawn(true);

    ModePluginDockWidget * dockwidget = new ModePluginDockWidget("Particles", mainWindow());

	ParticlesWidget * pw = new ParticlesWidget();
    widget = pw;

    dockwidget->setWidget( widget );
    mainWindow()->addDockWidget(Qt::RightDockWidgetArea, dockwidget);

	// General Tests
	connect(pw->ui->testButton, &QPushButton::released, [=]{
		//for(auto p : sphere_fibonacci_points( 100 ))drawArea()->drawPoint(p, 5);
		//drawArea()->update();

	});

	// Load and process shapes:
	connect(pw->ui->loadShapes, &QPushButton::released, [=]{
		QStringList files = QFileDialog::getOpenFileNames(nullptr, "Open Shapes", "", "All Supported (*.obj *.off)");
		for(auto filename : files){
			SurfaceMeshModel fromMesh;
			fromMesh.read( filename.toStdString() );
			pw->pmeshes.push_back( new ParticleMesh( &fromMesh, 256 ) );
			//document()->addModel(pw->pmeshes.back()->surface_mesh);
		}
		emit( pw->shapesLoaded() );
	});

	// Post-processing
	connect(pw, &ParticlesWidget::shapesLoaded, [=]{
		//for(auto s : pw->pmeshes) document()->addModel( s->surface_mesh );

		mainWindow()->setStatusBarMessage("Shapes loaded, now processing..");
		qApp->processEvents();

        int perSampleRaysCount = 360;
        std::vector< Eigen::Vector3d > sampledRayDirections = sphere_fibonacci_points( perSampleRaysCount );

		typedef clustering::l2norm_squared< std::vector<double> > dist_fn;

		if( true )
		{
			for(auto & s : pw->pmeshes)
			{
				std::vector< Eigen::Vector3f > rayOrigins;
				std::vector< Eigen::Vector3f > rayDirections;

				for(auto & p : s->particles)
				{
					for(auto d : sampledRayDirections)
					{
						rayOrigins.push_back( p.pos.cast<float>() );
						rayDirections.push_back( d.cast<float>() );
					}
				}

				// Smooth mesh
				//SurfaceMeshHelper h(s->surface_mesh);
				//h.smoothVertexProperty<Vector3>(VPOINT, 3, Vector3(0,0,0));

				raytracing::Raytracing<Eigen::Vector3f> rt(s->surface_mesh, rayOrigins, rayDirections);

				mainWindow()->setStatusBarMessage( QString("Ray tracing: rays (%1) / time (%2 ms)").arg( rayOrigins.size() ).arg( rt.time ) );

				double maxDistSurface = -DBL_MAX;

				std::vector< std::vector<double> > & desc = (s->desc = std::vector< std::vector<double> >( s->particles.size() ));

				for(auto & p : s->particles)
				{
					std::vector<double> descriptor( perSampleRaysCount );
					for(int i = 0; i < perSampleRaysCount; i++){
						descriptor[i] = rt.hits[p.id * perSampleRaysCount + i].distance;
					}

					desc[p.id] = descriptor;
				}

				clustering::kmeans< std::vector< std::vector<double> >, dist_fn > km(desc, 6);
				km.run(100, 0.01);

				for(auto & p : s->particles)
				{
					// Test clustering
					p.flag = (int) km.clusters()[p.id];

					std::vector<double> descriptor = desc[p.id];

					p.alpha = *std::min_element( descriptor.begin(), descriptor.end() );
					maxDistSurface = std::max( maxDistSurface, p.alpha );

					//p.measure = ;
					int idx = std::max_element( descriptor.begin(), descriptor.end() ) - descriptor.begin();
					p.direction = sampledRayDirections[idx].normalized();
				}

				for(auto & p : s->particles)
				{
					p.alpha = pow(p.alpha / maxDistSurface, 2);
				}

				// Debug
				{
					starlab::LineSegments * vs = new starlab::LineSegments;
					for(auto particle : s->particles) 
					{		
						Eigen::Vector4d color(abs(particle.direction[0]), abs(particle.direction[1]), abs(particle.direction[2]), particle.alpha);
						color[0] *= color[0];color[1] *= color[1];color[2] *= color[2];

						vs->addLine(particle.pos,  Vector3(particle.pos + particle.direction * 0.005), QColor::fromRgbF(color[0],color[1],color[2],1));
					}
					//s->debug.push_back(vs);
				}
			}

		}

		/// [ Correspondence ] match particles
		{
			for(size_t i = 0; i < pw->pmeshes.size(); i++){
				for(size_t j = i+1; j < pw->pmeshes.size(); j++){
					NanoKdTree * itree = pw->pmeshes[i]->kdtree;
					NanoKdTree * jtree = pw->pmeshes[j]->kdtree;

					for(auto & iparticle : pw->pmeshes[i]->particles)
					{
						if( true )
						{
							iparticle.correspondence = jtree->closest( iparticle.relativePos );
						}
						else
						{
							// Experiment
							KDResults matches;
							jtree->ball_search( iparticle.relativePos, 0.2, matches );

							QMap<double, int> measures;
							for(auto p : matches) 
							{
								double weight = p.second;
								//double dist = dist_fn()( pw->pmeshes[i]->desc[iparticle.id], pw->pmeshes[j]->desc[p.first] );
								double dist = 1;
								measures[ weight * dist ] = p.first;
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
							itree->ball_search( jparticle.relativePos, 0.2, matches );

							QMap<double, int> measures;
							for(auto p : matches) 
							{
								double weight = p.second;
								//double dist = dist_fn()( pw->pmeshes[j]->desc[jparticle.id], pw->pmeshes[i]->desc[p.first] );
								double dist = 1;
								measures[ weight * dist ] = p.first;
							}
							jparticle.correspondence = measures[ measures.keys().front() ];
						}
					}
				}
			}
			pw->isReady = true;
		}

		drawArea()->update();
	});
}

void particles::decorate()
{
	ParticlesWidget * pwidget = (ParticlesWidget*) widget;
	if(!pwidget || !pwidget->isReady || pwidget->pmeshes.size() < 1) return;

	if(pwidget->pmeshes.size() < 2) 
	{
		// Evaluation
		for( auto s : pwidget->pmeshes )
		{
			s->drawParticles();
			s->drawDebug( *drawArea() );
		}

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
	for(auto pmesh : pwidget->pmeshes) largeBox.extend(pmesh->bbox);
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
