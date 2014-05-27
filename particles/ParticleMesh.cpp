#include "RenderObjectExt.h"
#include "ParticleMesh.h"
#include <QGLWidget>

ParticleMesh::ParticleMesh(const std::vector<Eigen::Vector3d> &fromPoints, double raidus) : 
	raidus(raidus), tranlsation(Eigen::Vector3d(0,0,0))
{
	kdtree = new NanoKdTree;

    for(auto point : fromPoints)
    {
        particles.push_back( Particle(point) );
		bbox.extend( point );

		kdtree->addPoint( point );
    }

	kdtree->build();

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

void ParticleMesh::drawParticles()
{
	glPointSize(10);
	glBegin(GL_POINTS);

	for(auto particle : particles)
	{
		glColorQt(starlab::qtJetColor(particle.measure));
		glVertex3dv( particle.pos.data() );
	}

	glEnd();
}
