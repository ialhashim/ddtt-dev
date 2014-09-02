#include "ParticleMesh.h"
#include "myglobals.h"

#include "RenderObjectExt.h"
#include <QGLWidget>

#include "bluenoise.h"

#include "kmeans.h"

Q_DECLARE_METATYPE(Eigen::Vector3d)

QVector<QColor> ParticleMesh::rndcolors = rndColors2(10000);

ParticleMesh::ParticleMesh(SurfaceMeshModel * mesh, int gridsize) : surface_mesh(NULL)
{
	if(!mesh) return;

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

		surface_mesh->garbage_collection();
		surface_mesh->triangulate();
		meregeVertices( surface_mesh );
	}

	// Remove outer most voxel
	//container.data.erase(std::remove_if(container.data.begin(), container.data.end(), 
	//	[](const VoxelData<Eigen::Vector3f> & vd) { return !vd.isOuter; }), container.data.end());

	// Insert particles
    for( auto voxel : grid.data )
    {
		Eigen::Vector3f point = grid.voxelPos(voxel.morton);

        Particle<Vector3> particle( point.cast<double>() );

        particle.id = particles.size();
		particle.morton = voxel.morton;
		mortonToParticleID[voxel.morton] = particle.id;

        particles.push_back( particle );
    }

	grid.findOccupied();

	// Cache adjacency
	cachedAdj.resize(particles.size());

	process();
}

Eigen::AlignedBox3d ParticleMesh::bbox()
{
	Eigen::AlignedBox3d box;
	for(auto & p : particles) box.extend(p.pos);
	return box;
}

void ParticleMesh::process()
{
	/// Relative z-value:
	//Eigen::AlignedBox3d box = bbox();
	//double bmin = box.min().z(), bmax = box.max().z();
	//for(auto & particle : particles) particle.measure = (particle.pos.z() - bmin) / (bmax - bmin);
	
	/// Normalized distance to floor:
	computeDistanceToFloor();
}

void ParticleMesh::drawParticles( qglviewer::Camera * camera )
{
	// Light setup
	if(false)
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

		QColor c = rndcolors[ particle.segment ];
		if(particle.flag == VIZ_WEIGHT) c = starlab::qtJetColor(particle.weight); 
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

std::vector< std::vector< std::vector<float> > > ParticleMesh::toGrid()
{
	size_t gridsize = grid.gridsize;

	std::vector< std::vector< std::vector<float> > > g;
	g.resize(gridsize, std::vector<std::vector<float> >(gridsize, std::vector<float>(gridsize, 1)));

	double gridlength = gridsize * grid.unitlength;

	for(auto & particle : particles)
	{
		Vector3 p = particle.pos - grid.translation.cast<double>();
		Vector3 delta = p / gridlength;
		Eigen::Vector3i gridpnt( delta.x() * gridsize, delta.y() * gridsize, delta.z() * gridsize );

		// skip outside of grid..
		if(gridpnt.x() < 0 || gridpnt.y() < 0 || gridpnt.z() < 0) continue; 
		if(gridpnt.x() > gridsize-1 || gridpnt.y() > gridsize-1 || gridpnt.z() > gridsize-1) continue;

		g[gridpnt[2]][gridpnt[1]][gridpnt[0]] = -1;
	}

	return g;
}

void ParticleMesh::distort()
{
	Eigen::AlignedBox3d box = bbox();

	for(auto & particle : particles)
	{
		double t = (particle.pos.x() - bbox().min().x()) / box.sizes().x();
		double d = sin( t * 10 );
		Vector3 delta( 0, 0, d );
		particle.pos += delta;
	}
}

ParticleMesh::~ParticleMesh()
{
	if(surface_mesh) delete surface_mesh;
}

SegmentGraph ParticleMesh::toGraph( SegmentGraph::vertices_set selected )
{
	SegmentGraph graph = SegmentGraph();

	int eidx = 0;
	std::vector<size_t> pindices;

	NanoKdTree tree;
	for(auto & p : particles) 
	{
		bool isInclude = true;

		if(!selected.empty() && selected.find(p.id) == selected.end()) 
			isInclude = false;

		if( isInclude ){
			tree.addPoint(p.pos);
			pindices.push_back(p.id);
		}
	}
	tree.build();

	for (auto & pid : pindices)
	{
		auto & p = particles[pid];

		KDResults matches;
		tree.ball_search(p.pos, grid.unitlength*1.01, matches);
		matches.erase(matches.begin()); // remove self

		for(auto match : matches)
		{
			double edge_weight = 1.0;
			double d2 = match.second;
			edge_weight = d2 /*std::sqrt(d2)*/;
			graph.AddEdge( uint(p.id), uint(match.first), edge_weight, eidx++ );
		}
	}

	return graph;
}

std::vector< double > ParticleMesh::agd( int numStartPoints )
{
	auto graph = toGraph();

	std::vector<double> sum_distances( particles.size(), 0.0 );

	// Random staring points
	std::vector<int> v;

	if( numStartPoints > 0 )
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dis(0, int(particles.size()-1));
		for(int i = 0; i < numStartPoints; i++)
			v.push_back( dis(gen) );

		// Avoid duplicates
		std::sort(v.begin(), v.end());
		auto last = std::unique(v.begin(), v.end());
		v.erase(last, v.end());
	}
	else
		for(size_t i = 0; i < particles.size(); i++) v.push_back(int(i));

	// Sum distances to other
	#pragma omp parallel for
	for(int pi = 0; pi < (int)v.size(); pi++){
		auto curGraph = graph;
		curGraph.DijkstraComputePaths(v[pi]);
		for(size_t pj = 0; pj < particles.size(); pj++)
			sum_distances[pj] += curGraph.min_distance[pj];
	}

	// Average
	auto avg_distances = sum_distances;
	for(auto & p : particles) avg_distances[p.id] = sum_distances[p.id] / particles.size();

	// Bounds
	double minDist = *std::min_element(avg_distances.begin(),avg_distances.end());
	double maxDist = *std::max_element(avg_distances.begin(),avg_distances.end());

	// Normalize
	for(auto & v : avg_distances) v = (v-minDist) / (maxDist-minDist);

	return avg_distances;
}

SpatialHash< Vector3, Vector3::Scalar > ParticleMesh::spatialHash()
{
	std::vector<Vector3> positions;
	for(auto & p : particles) positions.push_back(p.pos);

	return SpatialHash< Vector3, Vector3::Scalar >( positions, grid.unitlength );
}

std::vector<size_t> ParticleMesh::randomSamples( int numSamples, bool isSpread )
{
	std::vector<size_t> samples;
	std::set<size_t> set;

	if( numSamples >= particles.size() )
	{
		for(size_t i = 0; i < particles.size(); i++) set.insert(i);
	}
	else if( isSpread )
	{
		Eigen::AlignedBox3d box = bbox();

		double spreadFactor = pow(std::max(1.0, double(particles.size()) / numSamples), 1.0 / 3.0);

		std::vector<Vector3> samples;
		bluenoise_sample<3,double,Vector3>( grid.unitlength * spreadFactor, box.min(), box.max(), samples );

		size_t gridsize = grid.gridsize;
		double gridlength = gridsize * grid.unitlength;

		for(auto s : samples)
		{
			Vector3 p = s - grid.translation.cast<double>();
			Vector3 delta = p / gridlength;
			Eigen::Vector3i gridpnt( delta.x() * gridsize, delta.y() * gridsize, delta.z() * gridsize );

			// skip outside of grid..
			if(gridpnt.x() < 0 || gridpnt.y() < 0 || gridpnt.z() < 0) continue; 
			if(gridpnt.x() > gridsize-1 || gridpnt.y() > gridsize-1 || gridpnt.z() > gridsize-1) continue;

			uint64_t m = mortonEncode_LUT(gridpnt.z(),gridpnt.y(),gridpnt.x());
			if( grid.occupied[m] ) set.insert(mortonToParticleID[m]);
		}
	}
	else
	{
		std::vector<int> idx;
		for(size_t i = 0; i < particles.size(); i++) idx.push_back(i);
		for(int idx : random_sampling<int>(idx, numSamples))
			set.insert(idx);
	}

	for(auto i : set) samples.push_back(i);

	return samples;
}

void ParticleMesh::computeDistanceToFloor()
{
	std::set<uint> sources;
	for(auto & p : particles){
		unsigned int x,y,z;
		mortonDecode(p.morton,z,y,x);
		if( z > 1 ) continue;
		sources.insert( uint(p.id) );

		p.flag = ParticleFlags::FLOOR;
	}

	auto g = toGraph();
	g.DijkstraComputePathsMany( sources );

	double minVal = *std::min_element(g.min_distance.begin(),g.min_distance.end());
	auto tipPoint = std::max_element(g.min_distance.begin(),g.min_distance.end());
	double maxVal = *tipPoint;
	double range = maxVal - minVal;

	// Normalize
	for(auto & p : particles)
		p.measure = (g.min_distance[p.id] - minVal) / range;

	// Keep a path from ground to furthest tip point
	auto path = g.DijkstraGetShortestPathsTo( tipPoint - g.min_distance.begin() );

	for(auto p : path) if(p < particles.size()) pathFromFloor.push_back(p);
}

QMap< unsigned int, SegmentGraph > ParticleMesh::segmentToComponents( SegmentGraph fromGraph, SegmentGraph & neiGraph )
{
	QMap< unsigned int, SegmentGraph > result;
	
	// Make a copy of the graph we will modify
	if(fromGraph.IsEmpty()) return result;

	std::vector<SegmentGraph::Edge> cutEdges;

	// Remove edges between two nodes having different segments
	for(auto e : fromGraph.GetEdgesSet())
	{
		int s1 = particles[e.index].segment;
		int s2 = particles[e.target].segment;

		if(s1 != s2) 
		{
			fromGraph.removeEdge( e.index, e.target );
			cutEdges.push_back( e );
		}
	}

	auto all_parts = fromGraph.toConnectedParts();
	auto edgeToParts = [&]( SegmentGraph::Edge & e ){
		std::pair<uint,uint> result;
		for(auto & part : all_parts) if(part.IsHasVertex(e.index)) {result.first = part.uid; part.sid = particles[e.index].segment;}
		for(auto & part : all_parts) if(part.IsHasVertex(e.target)) {result.second = part.uid; part.sid = particles[e.target].segment;}
		return result;
	};

	auto edgeCenter = [&]( SegmentGraph::Edge & e ){
		return Vector3(0.5 * ( particles[e.index].pos + particles[e.target].pos ));
	};

	auto edgeDirection = [&]( SegmentGraph::Edge & e ){
		return Vector3(( particles[e.index].pos - particles[e.target].pos ).normalized());
	};

	// Boundary edges
	typedef std::pair<uint,uint> EdgeID;
	QMap< EdgeID, std::vector<SegmentGraph::Edge> > boundaryEdges;

	for(auto & e : cutEdges)
	{
		// First find parts containing an original edge
		auto parts = edgeToParts( e );

		auto key = EdgeID(std::min(parts.first, parts.second), std::max(parts.first, parts.second));
		boundaryEdges[ key ].push_back( e );
	}

	for(auto key : boundaryEdges.keys())
	{
		auto count = boundaryEdges[key].size();

		neiGraph.AddEdge( key.first, key.second, count);

		/// Assign properties:

		// Boundary plane and center:
		auto e = boundaryEdges[key].front();
		Vector3 b_normal = edgeDirection(e);
		Vector3 b_center = edgeCenter(e);

		// Better estimation
		if(count > 3)
		{
			Eigen::MatrixXd points( count, 3 ); int r = 0;
			for(auto & e : boundaryEdges[key]) points.row(r++) = edgeCenter(e);
			
			// Fit a plane to a set of points with help of SVD
			b_center = Vector3(points.colwise().mean());
			points = points.rowwise() - b_center.transpose();
			Eigen::JacobiSVD<Eigen::MatrixXd> svd(points, Eigen::ComputeThinU | Eigen::ComputeThinV);
			b_normal = svd.matrixV().col(2).normalized();
		}

		QVariant normalVar; normalVar.setValue( b_normal );
		QVariant centerVar; centerVar.setValue( b_center );

		neiGraph.SetEdgeProperty( key.first, key.second, "normal", normalVar );
		neiGraph.SetEdgeProperty( key.first, key.second, "center", centerVar );
	}

	// Disconnected graphs
	if( boundaryEdges.empty() ){
		for(auto & part : all_parts) 
			neiGraph.AddVertex( part.uid );
	}

	// Collect resulting parts indexed by their unique identifier
	for(auto & part : all_parts)
		result[part.uid] = part;

	return result;
}

std::vector<size_t> ParticleMesh::neighbourhood( Particle<Vector3> & p, int step )
{
	if(!cachedAdj[p.id][step].empty())
		return cachedAdj[p.id][step];

	std::vector<size_t> result;
	std::queue<size_t> toSee;
	std::set<size_t> visisted;

	unsigned int x,y,z;
	mortonDecode(p.morton, x, y, z);
	Eigen::Vector3i c0(x, y, z);

	toSee.push( p.id );

	while( !toSee.empty() )
	{
		size_t current = toSee.front();
		toSee.pop();
		visisted.insert(current);

		mortonDecode(particles[current].morton, x, y, z);

		for(int u = -1; u <= 1; u++){
			for(int v = -1; v <= 1; v++){
				for(int w = -1; w <= 1; w++){
					Eigen::Vector3i c(x + u, y + v, z + w);
					if(c.x() < 0 || c.y() < 0 || c.z() < 0) continue;
					if(c.x() > grid.gridsize-1 || c.y() > grid.gridsize-1 || c.z() > grid.gridsize-1) continue;;

					uint64_t m = mortonEncode_LUT( c.x(), c.y(), c.z() );
					if( m == p.morton) continue;
					if( !grid.occupied[m] ) continue;

					auto pid = mortonToParticleID[m];

					auto distance = (c-c0).norm();

					if( visisted.find(pid) == visisted.end() && distance <= step )
					{
						result.push_back( pid );
						toSee.push( pid );
						visisted.insert( pid );
					}
				}
			}
		}
	}

	cachedAdj[p.id][step] = result;

	return result;
}

std::vector< std::pair< double, size_t > > ParticleMesh::closestParticles(const Vector3 & point, double threshold)
{
	typedef std::pair< double, size_t > DistParticleID;

	std::vector< DistParticleID > result;

	for(auto & p : particles){
		double dist = (p.pos-point).norm();
		if(dist <= threshold) result.push_back( std::make_pair(dist, p.id) );
	}

	std::sort( result.begin(), result.end(), [](const DistParticleID& a, const DistParticleID & b){ return a.first < b.first; } );

	return result;
}

void ParticleMesh::cluster( int K, const std::vector<size_t> & seeds, bool use_l1_norm, bool showSeeds )
{
	if(!particles.size()) return;

	// Clustering engine:
	clustering::kmeans< std::vector< VectorFloat >, clustering::lpnorm< VectorFloat > > km(desc, K, clustering::KmeansInitPlusPlus);

	// Distance measure:
	if(use_l1_norm) clustering::lpnorm_p = 1;
	else clustering::lpnorm_p = 2;

	// Custom seeding:
	if( seeds.size() )
	{
		km._centers.clear();
		for(auto pid : seeds) km._centers.push_back( desc[pid] );
	}

	// [DEBUG] seed points
	if( showSeeds )
	{
		starlab::PointSoup * ps = new starlab::PointSoup(20);
		if(!seeds.empty()) for(auto pid : seeds) ps->addPoint(particles[pid].pos, Qt::black);
		else for(auto pid : km.initindices) ps->addPoint(particles[pid].pos, Qt::black);
		debug.push_back(ps);
	}

	// Parameters:
	int numIterations = 1000;
	double minchangesfraction = 0.005;

	km.run(numIterations, minchangesfraction);

	// Record centers
	this->cluster_centers = km.centers();

	// Assign particles to classes:
	#pragma omp parallel for
	for(int pi = 0; pi < (int)particles.size(); pi++){
		auto & p = particles[pi];
		p.segment = (int) km.clusters()[p.id];
	}
	
	/*bool isClusterShapesTogether = false;
	if( isClusterShapesTogether )
	{
		// Collect descriptors
		std::vector< std::vector<float> > allDesc;
		for(auto & s : pw->pmeshes){
			if( pw->ui->bandvalue() > 0 )
				allDesc.insert(allDesc.end(), sig.begin(), sig.end());
			else
				allDesc.insert(allDesc.end(), desc.begin(), desc.end());
		}

		// Cluster
		clustering::kmeans< std::vector< std::vector<float> >, dist_fn > km(allDesc, K);
		km.run(numIterations, minchangesfraction);

		// Assign classes
		int pi = 0;
		for(auto & s : pw->pmeshes)
			for(auto & p : particles)
				p.segment = (int) km.clusters()[pi++];
	}*/
}

void ParticleMesh::shrinkSmallerClusters()
{
	std::vector<int> newSegmentAssignment(particles.size());

	#pragma omp parallel for
	for(int pi = 0; pi < (int)particles.size(); pi++){
		auto & p = particles[pi];

		std::map<size_t,size_t> clusterHistogram;
		for( auto pj : neighbourhood(p) ) 
			clusterHistogram[ particles[pj].segment ]++;

		std::vector< std::pair<size_t,size_t> > ch( clusterHistogram.begin(), clusterHistogram.end() );
		std::sort(ch.begin(),ch.end(),[](std::pair<size_t,size_t> a, std::pair<size_t,size_t> b){ return a.second < b.second; });

		// majority rule
		newSegmentAssignment[p.id] = (int)ch.back().first;
		//p.segment = ch.back().first; // iterative arbitrary growing
	}

	for(auto & p : particles) p.segment = newSegmentAssignment[p.id];
}

Particle<Vector3> ParticleMesh::pointToParticle( const Vector3 & point )
{
	size_t gridsize = grid.gridsize;
	double gridlength = gridsize * grid.unitlength;

	Vector3 p = point - grid.translation.cast<double>();
	Vector3 delta = p / gridlength;
	Eigen::Vector3i gridpnt( delta.x() * gridsize, delta.y() * gridsize, delta.z() * gridsize );
	uint64_t m = mortonEncode_LUT(gridpnt.z(),gridpnt.y(),gridpnt.x());

	// Check: not occupied or outside of grid
	bool isOccupied = grid.occupied[m];
	bool isOutsideGrid = (gridpnt.x() < 0 || gridpnt.y() < 0 || gridpnt.z() < 0) 
			|| (gridpnt.x() > gridsize-1 || gridpnt.y() > gridsize-1 || gridpnt.z() > gridsize-1);
	if( !isOccupied || isOutsideGrid ) return Particle<Vector3>(Vector3(DBL_MAX,DBL_MAX,DBL_MAX));

	return particles[ mortonToParticleID[m] ];
}

std::vector< Vector3 > ParticleMesh::particlesCorners( SegmentGraph::vertices_set selected )
{
	double gridunit = grid.unitlength;
	double halfgridunit = 0.5 * gridunit;

	std::vector< Vector3 > points;
	for(auto & particle : particles) 
	{
		bool isInclude = true;
		if(!selected.empty() && selected.find(particle.id) == selected.end()) isInclude = false;
		if( isInclude ) 
		{
			Vector3 corner = particle.pos - Vector3(halfgridunit,halfgridunit,halfgridunit);
			QVector<Vector3> cornersBottom, cornersTop;

			// Bottom row
			cornersBottom.push_back(corner);
			cornersBottom.push_back(corner + Vector3(gridunit,0,0));
			cornersBottom.push_back(corner + Vector3(gridunit,gridunit,0));
			cornersBottom.push_back(corner + Vector3(0,gridunit,0));

			// Top row
			for(auto & p : cornersBottom) cornersTop << (p + Vector3(0,0,gridunit));

			// Combined
			for(auto & p : (cornersBottom + cornersTop))
				points.push_back( p );
		}
	}
	return points;
}

std::vector<size_t> ParticleMesh::specialSeeding( SeedType seedType, int K, SegmentGraph::vertices_set selected )
{
	std::set<size_t> seeds;

	if( seedType == GROUND )
	{
		// Seed based on ground distance
		for(int i = 0; i < K; i++)
		{
			double t = double(i) / (K-1);
			int idx = t * (pathFromFloor.size()-1);
			seeds.insert( pathFromFloor[idx] );
		}
	}

	if( seedType == DESCRIPTOR )
	{
		if(selected.empty())
			for(size_t i = 0; i < particles.size();i++)
				selected.insert(i);

		// Collect sorted list of magnitudes
		QMap<double,size_t> descParticle;
		for(auto particleID : selected)
		{
			auto d = Eigen::Map<Eigen::VectorXf>(&desc[particleID][0], desc[particleID].size());
			descParticle[ d.norm() * this->particles[particleID].flat ] = particleID;
		}
		QVector<size_t> sorted = descParticle.values().toVector();

		// Seed based on magnitude of descriptor
		for(int i = 0; i < K; i++)
		{
			double t = double(i) / (K-1);
			int idx = t * (sorted.size()-1);
			seeds.insert( sorted[idx] );
		}
	}

	std::vector<size_t> seedsVector;
	for(auto seed : seeds) seedsVector.push_back(seed);

	return seedsVector;
}

SurfaceMeshModel * ParticleMesh::meshPoints( const std::vector<Eigen::Vector3f> & points ) const
{
	SurfaceMeshModel * m = new SurfaceMeshModel("meshed.obj","meshed");

	size_t gridsize = grid.gridsize;
	double gridlength = gridsize * grid.unitlength;

	std::vector<char> occupied(grid.occupied.size(), EMPTY_VOXEL);
	std::set<uint64_t> voxels;

	for(const auto & point : points)
	{
		Vector3 p = point.cast<double>() - grid.translation.cast<double>();
		Vector3 delta = p / gridlength;
		Eigen::Vector3i gridpnt( delta.x() * gridsize, delta.y() * gridsize, delta.z() * gridsize );
		uint64_t m = mortonEncode_LUT(gridpnt.z(),gridpnt.y(),gridpnt.x());

		occupied[m] = FULL_VOXEL;
		voxels.insert(m);
	}

	// Prepare set of quads
	std::vector< std::pair<uint64_t, Eigen::Vector3i> > allQuads;

	for(auto voxel_morton : voxels)
	{					
		unsigned int x,y,z;
		mortonDecode(voxel_morton, x, y, z);

		std::vector<uint64_t> neigh;
		if(x < gridsize - 1) neigh.push_back( mortonEncode_LUT(x+1,y,z) );
		if(y < gridsize - 1) neigh.push_back( mortonEncode_LUT(x,y+1,z) );
		if(z < gridsize - 1) neigh.push_back( mortonEncode_LUT(x,y,z+1) );

		if(x > 0) neigh.push_back( mortonEncode_LUT(x-1,y,z) );
		if(y > 0) neigh.push_back( mortonEncode_LUT(x,y-1,z) );
		if(z > 0) neigh.push_back( mortonEncode_LUT(x,y,z-1) );

		// Inside / outside
		for(auto n : neigh)
		{
			if(occupied[n] != occupied[voxel_morton])
			{
				unsigned int v[3], w[3];
				mortonDecode(voxel_morton, v[0], v[1], v[2]);
				mortonDecode(n, w[0], w[1], w[2]);
				Eigen::Vector3i direction (int(w[2])-int(v[2]), int(w[1])-int(v[1]), int(w[0])-int(v[0]));
				allQuads.push_back( std::make_pair(voxel_morton, direction) );
			}
		}

		// Edge of grid
		unsigned int v[3];
		mortonDecode(voxel_morton, v[2], v[1], v[0]);
		bool isBoundary = (v[0] == 0 || v[1] == 0 || v[2] == 0) || 
			(v[0] == gridsize-1 || v[1] == gridsize-1 || v[2] == gridsize-1);
		if( !isBoundary ) continue;

		for(int i = 0; i < 3; i++)
		{
			Eigen::Vector3i d(0,0,0);
			if(v[i] == 0) d[i] = -1;
			else if(v[i] == gridsize-1) d[i] = 1;
			if(d[0]!=0||d[1]!=0||d[2]!=0) 
				allQuads.push_back( std::make_pair(voxel_morton, d) );
		}
	}

	// Generate surface quads in world coordinates
	double unitlength = grid.unitlength;	
	Vector3 delta = grid.translation.cast<double>() + ( 0.5 * Vector3(unitlength,unitlength,unitlength) );

	for(auto p : allQuads)
	{			
		unsigned int v[3];
		mortonDecode(p.first, v[0], v[1], v[2]);
		std::vector<Vector3> quad = voxelQuad<Vector3>( p.second, unitlength );
		std::vector<Vertex> quad_verts;
		for(auto p : quad){
			p += Vector3(v[2] * unitlength, v[1] * unitlength, v[0] * unitlength) + delta;
			quad_verts.push_back( Vertex(m->n_vertices()) );
			m->add_vertex( p );
		}

		m->add_face(quad_verts);
	}

	return m;
}

Vector3 ParticleMesh::mainDirection( size_t particleID )
{
	std::vector< bool > isVisisted( usedDirections.size(), false );

	float maxDiam = -FLT_MAX;
	size_t maxIdx = 0;

	for(size_t idx = 0; idx < usedDirections.size(); idx++)
	{
		if(isVisisted[antiRays[idx]]) continue;
		isVisisted[antiRays[idx]] = true;
		auto & pdesc = desc[particleID];

		auto diam = pdesc[idx] + pdesc[antiRays[idx]];
		if(diam > maxDiam)
		{
			maxDiam = diam;
			maxIdx = idx;
		}
	}

	return usedDirections[maxIdx];
}

std::vector< Vector3 > ParticleMesh::particlesPositions(const std::set<unsigned int> & P)
{
	std::vector<Vector3> points;
	for(auto p : P) points.push_back(particles[p].pos);
	return points;
}

void ParticleMesh::serialize(std::ostream& os) const
{ 
	// Particles
	os << particles.size() << std::endl;
	for(auto & p : particles) os << p;

	// Descriptors
	os << std::endl;
	for(auto & p : particles) os << desc[p.id]; 

	// Grid
	os << std::endl;
	os << grid.gridsize << " " << grid.unitlength << " ";
}

void ParticleMesh::deserialize(std::istream& is)
{
	size_t particleCount;
	is >> particleCount;
	particles.resize(particleCount);
	desc.resize(particleCount);
	for(size_t i = 0; i < particleCount; i++) is >> particles[i];
	for(size_t i = 0; i < particleCount; i++) is >> desc[i];
	is >> grid.gridsize >> grid.unitlength;
}
