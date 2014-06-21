#include "ParticleMesh.h"

#include "voxelization.h"

#include "RenderObjectExt.h"
#include <QGLWidget>

inline QVector<QColor> rndColors(int count){
	QVector<QColor> c;
	for(int i = 0; i < count; i++) c << starlab::qRandomColor3();
	return c;
}

QVector<QColor> ParticleMesh::rndcolors = rndColors(512);

ParticleMesh::ParticleMesh(SurfaceMeshModel * mesh, int gridsize, double particle_raidus) : surface_mesh(NULL),
	raidus(particle_raidus), tranlsation(Eigen::Vector3d(0,0,0))
{
	// Voxelization
	grid = ComputeVoxelization<VoxelVector>(mesh, gridsize, true, true);

	// Build voxelization mesh
	QString objectName = mesh->name;
	surface_mesh = new SurfaceMeshModel(objectName + ".obj", objectName);
	{
		int voffset = 0;
		for(auto q : grid.quads){
			std::vector<Vertex> quad_verts;
			for(int i = 0; i < 4; i++){
				surface_mesh->add_vertex( q[i].cast<double>() );
				quad_verts.push_back( Vertex(voffset + i) );
			}
			surface_mesh->add_face( quad_verts );
			voffset += 4;
		}

		surface_mesh->triangulate();
		meregeVertices( surface_mesh );

		for(auto v : surface_mesh->vertices())
			if(surface_mesh->is_isolated(v))
				surface_mesh->remove_vertex(v);

		surface_mesh->garbage_collection();
	}

	// Remove outer most voxel
	//container.data.erase(std::remove_if(container.data.begin(), container.data.end(), 
	//	[](const VoxelData<Eigen::Vector3f> & vd) { return !vd.isOuter; }), container.data.end());

	// Insert particles
    for( auto voxel : grid.data )
    {
		Eigen::Vector3f point = grid.voxelPos(voxel.morton);

        Particle particle( point.cast<double>() );

        particle.id = particles.size();
		particle.morton = voxel.morton;
		mortonToParticleID[voxel.morton] = particle.id;

        particles.push_back( particle );
		bbox.extend( point.cast<double>() );
    }

	grid.findOccupied();

	// KD-tree
	{
		kdtree = new NanoKdTree;
		Vector3 sizes = bbox.sizes();
		for( auto & particle : particles )
		{
			Vector3 mapped = (particle.pos - bbox.min());
			for(int i = 0; i < 3; i++) mapped[i] /= sizes[i];

			particle.relativePos = mapped;
			kdtree->addPoint( particle.relativePos );
		}
		kdtree->build();
	}

	process();
}

void ParticleMesh::process()
{
	double bmin = bbox.min().z(), bmax = bbox.max().z();

	for(auto & particle : particles)
	{
		particle.measure = (particle.pos.z() - bmin) / (bmax - bmin);
	}
}

void ParticleMesh::drawParticles( qglviewer::Camera * camera )
{
	// Light setup
	{
		GLfloat lightColor[] = {0.9f, 0.9f, 0.9f, 1.0f};
		glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor);

		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHTING);

		// Specular
		glEnable(GL_COLOR_MATERIAL);
		glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
		float specReflection[] = { 0.8f, 0.8f, 0.8f, 1.0f };
		GLfloat high_shininess [] = { 100.0 }; 
		glMaterialfv(GL_FRONT, GL_SPECULAR, specReflection);
		glMaterialfv (GL_FRONT, GL_SHININESS, high_shininess); 
		glMateriali(GL_FRONT, GL_SHININESS, 56);
	}

	// Setup
	qglviewer::Vec pos = camera->position();
	qglviewer::Vec dir = camera->viewDirection();
	qglviewer::Vec revolve = camera->revolveAroundPoint();
	Vector3 eye(pos[0],pos[1],pos[2]);
	Vector3 direction(dir[0],dir[1],dir[2]);
	Vector3 center(revolve[0],revolve[1],revolve[2]);

	// Sort particles
	std::map<size_t,double> distances;
	double minDist = DBL_MAX, maxDist = -DBL_MAX;

	for(size_t i = 0; i < particles.size(); i++)
	{
		distances[i] = abs((particles[i].pos - eye).dot(direction));
		minDist = std::min(minDist, distances[i]);
		maxDist = std::max(maxDist, distances[i]);
	}

	// Sort
	std::vector<std::pair<size_t,double> > myVec(distances.begin(), distances.end());
	std::sort( myVec.begin(), myVec.end(), [](std::pair<size_t,double> a, std::pair<size_t,double> b){ return a.second < b.second; } );

	glDisable( GL_LIGHTING );
	//glEnable(GL_LIGHTING);
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glDepthMask(GL_FALSE);

	double camDist = (eye - center).norm();
	double ratio = 1.0 / camDist;

	glPointSize( std::max(2, std::min(10, int(8.0 * ratio))) );

	glBegin(GL_POINTS);

	for(auto pi : myVec)
	{
		auto & particle = particles[pi.first];
		//QColor c = starlab::qtJetColor(particle.measure);
		//Eigen::Vector4d color( c.redF(), c.greenF(), c.blueF(), particle.alpha );
		//Eigen::Vector4d color(abs(particle.direction[0]), abs(particle.direction[1]), abs(particle.direction[2]), particle.alpha);
		//color[0] *= color[0];color[1] *= color[1];color[2] *= color[2];
		//glColor4dv( color.data() );

		if(camDist > 0.5)
			particle.alpha = 1.0;
		else
			particle.alpha = std::max(0.3, 1.0 - ((pi.second - minDist) / (maxDist-minDist)));
		
		// Fake normals
		//Vector3 d = (particle.pos - eye).normalized();
		//Vector3 n = -direction;
		//Vector3 pn = (d - 2*(d.dot(n)) * n).normalized();
		//glNormal3dv( pn.data() );

		QColor c = rndcolors[ particle.flag ];
		Eigen::Vector4d color(c.redF(),c.greenF(),c.blueF(), particle.alpha);
		glColor4dv( color.data() );
		glVertex3dv( particle.pos.data() );
	}

	glEnd();

	//glDepthMask(GL_TRUE);
	glEnable(GL_LIGHTING);
}

void ParticleMesh::drawDebug(QGLWidget & widget)
{
	for(auto d : debug) d->draw( widget );
}

ParticleMesh::~ParticleMesh()
{
	if(surface_mesh) delete surface_mesh;
	if(kdtree) delete kdtree;
}
