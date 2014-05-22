// Adapted from https://github.com/Forceflow/ooc_svo_builder
#pragma once
#include <stack>

#include <Eigen/Core>
#include "SurfaceMeshModel.h"

#include "morton.h"

template <typename T>
struct AABox {
	T min;
	T max;
	AABox(): min(T()), max(T()){}
	AABox(T min, T max): min(min), max(max){}

	// Assumes 'T' is floating point
	AABox(AABox a, AABox b){
		min = T(DBL_MAX, DBL_MAX, DBL_MAX);
		max = -min;
		for(int i = 0; i < 3; i++){
			min[i] = std::min( min[i], std::min( a.min[i], b.min[i] ) );
			max[i] = std::max( max[i], std::max( a.max[i], b.max[i] ) );
		}
	}
};

// Intersection methods
inline AABox<Vector3> computeBoundingBox(const Vector3 &v0, const Vector3 &v1, const Vector3 &v2){
	AABox<Vector3> answer; 
	answer.min[0] = std::min(v0[0],std::min(v1[0],v2[0]));
	answer.min[1] = std::min(v0[1],std::min(v1[1],v2[1]));
	answer.min[2] = std::min(v0[2],std::min(v1[2],v2[2]));
	answer.max[0] = std::max(v0[0],std::max(v1[0],v2[0]));
	answer.max[1] = std::max(v0[1],std::max(v1[1],v2[1]));
	answer.max[2] = std::max(v0[2],std::max(v1[2],v2[2]));
	return answer;
}

inline Vector3 average3Vec(const Vector3 v0, const Vector3 v1, const Vector3 v2){
	Vector3 answer;
	for (size_t i = 0; i < 3; i++){
		answer[i] = (v0[i] + v1[i] + v2[i]) / 3.0;
	}
	return answer;
}

template <typename T> T clampval(const T& value, const T& low, const T& high) {
  return value < low ? low : (value > high ? high : value); 
}

// sizeof(VoxelData would work too)
const size_t VOXELDATA_SIZE = sizeof(uint64_t)+2 * (3 * sizeof(double));

#define EMPTY_VOXEL 0
#define FULL_VOXEL 1
#define OUTER_VOXEL 2

#define X 0
#define Y 1
#define Z 2

// This struct defines VoxelData for our voxelizer.
// This is the main memory hogger: the less data you store here, the better.
struct VoxelData{
	uint64_t morton;
	Vector3 color;
	Vector3 normal;
	
	VoxelData() : morton(0), normal(Vector3()), color(Vector3()){}
	VoxelData(uint64_t morton, Vector3 normal = Vector3(0,0,0), Vector3 color = Vector3(0,0,0)) : morton(morton), normal(normal), color(color){}

	bool operator > (const VoxelData &a) const{
		return morton > a.morton;
	}

	bool operator < (const VoxelData &a) const{
		return morton < a.morton;
	}
};

struct VoxelContainer{
	std::vector<VoxelData> data;
	Vector3 translation;
};

struct BasicTriangle{
	BasicTriangle() { counter = 0; v0_color = v1_color = v2_color = Vector3(0,0,0); }
	void setPoint(Vector3 p){ if(counter == 0) v0 = p; if(counter == 1) v1 = p; if(counter == 2) v2 = p; counter++; }
	Vector3 v0, v1, v2;
	Vector3 v0_color, v1_color, v2_color;
	int counter;
};

// Implementation of algorithm from http://research.michael-schwarz.com/publ/2010/vox/ (Schwarz & Seidel)
// Adapted for mortoncode -based subgrids

inline void voxelize_schwarz_method(SurfaceMeshModel * mesh, const uint64_t morton_start, const uint64_t morton_end, 
		const double unitlength, std::vector<char> & voxels, vector<VoxelData> &data, size_t &nfilled) 
{
	voxels.clear();
	voxels.resize(morton_end - morton_start, EMPTY_VOXEL);

	data.clear();

	// compute partition min and max in grid coords
	AABox< Eigen::Matrix<unsigned int, 3, 1> > p_bbox_grid;
	mortonDecode(morton_start, p_bbox_grid.min[2], p_bbox_grid.min[1], p_bbox_grid.min[0]);
	mortonDecode(morton_end - 1, p_bbox_grid.max[2], p_bbox_grid.max[1], p_bbox_grid.max[0]);

	// COMMON PROPERTIES FOR ALL TRIANGLES
	double unit_div = 1.0 / unitlength;
	Vector3 delta_p(unitlength, unitlength, unitlength);

	Vector3VertexProperty points = mesh->vertex_coordinates();
	Vector3FaceProperty fnormals = mesh->face_normals();

	// voxelize every triangle
	for(auto f : mesh->faces())
	{
		BasicTriangle t;
		for(auto vi : mesh->vertices(f)) t.setPoint(points[vi]);

		// compute triangle bbox in world and grid
		AABox<Vector3> t_bbox_world = computeBoundingBox(t.v0, t.v1, t.v2);
		AABox<Eigen::Vector3i> t_bbox_grid;
		t_bbox_grid.min[0] = (int)(t_bbox_world.min[0] * unit_div);
		t_bbox_grid.min[1] = (int)(t_bbox_world.min[1] * unit_div);
		t_bbox_grid.min[2] = (int)(t_bbox_world.min[2] * unit_div);
		t_bbox_grid.max[0] = (int)(t_bbox_world.max[0] * unit_div);
		t_bbox_grid.max[1] = (int)(t_bbox_world.max[1] * unit_div);
		t_bbox_grid.max[2] = (int)(t_bbox_world.max[2] * unit_div);

		// clamp
		t_bbox_grid.min[0] = clampval<int>(t_bbox_grid.min[0], p_bbox_grid.min[0], p_bbox_grid.max[0]);
		t_bbox_grid.min[1] = clampval<int>(t_bbox_grid.min[1], p_bbox_grid.min[1], p_bbox_grid.max[1]);
		t_bbox_grid.min[2] = clampval<int>(t_bbox_grid.min[2], p_bbox_grid.min[2], p_bbox_grid.max[2]);
		t_bbox_grid.max[0] = clampval<int>(t_bbox_grid.max[0], p_bbox_grid.min[0], p_bbox_grid.max[0]);
		t_bbox_grid.max[1] = clampval<int>(t_bbox_grid.max[1], p_bbox_grid.min[1], p_bbox_grid.max[1]);
		t_bbox_grid.max[2] = clampval<int>(t_bbox_grid.max[2], p_bbox_grid.min[2], p_bbox_grid.max[2]);

		// COMMON PROPERTIES FOR THE TRIANGLE
		Vector3 e0 = t.v1 - t.v0;
		Vector3 e1 = t.v2 - t.v1;
		Vector3 e2 = t.v0 - t.v2;
		Vector3 to_normalize = (e0).cross(e1);
		Vector3 n = (to_normalize).normalized(); // triangle normal
		// PLANE TEST PROPERTIES
		Vector3 c = Vector3(0.0, 0.0, 0.0); // critical point
		if (n[X] > 0) { c[X] = unitlength; }
		if (n[Y] > 0) { c[Y] = unitlength; }
		if (n[Z] > 0) { c[Z] = unitlength; }
		double d1 = n.dot(c - t.v0);
		double d2 = n.dot((delta_p - c) - t.v0);
		// PROJECTION TEST PROPERTIES
		// XY plane
		Eigen::Vector2d n_xy_e0 = Eigen::Vector2d(-1.0*e0[Y], e0[X]);
		Eigen::Vector2d n_xy_e1 = Eigen::Vector2d(-1.0*e1[Y], e1[X]);
		Eigen::Vector2d n_xy_e2 = Eigen::Vector2d(-1.0*e2[Y], e2[X]);
		if (n[Z] < 0.0) {
			n_xy_e0 = -1.0 * n_xy_e0;
			n_xy_e1 = -1.0 * n_xy_e1;
			n_xy_e2 = -1.0 * n_xy_e2;
		}
		double d_xy_e0 = (-1.0 * (n_xy_e0.dot(Eigen::Vector2d(t.v0[X], t.v0[Y])))) + std::max(0.0, unitlength*n_xy_e0[0]) + std::max(0.0, unitlength*n_xy_e0[1]);
		double d_xy_e1 = (-1.0 * (n_xy_e1.dot(Eigen::Vector2d(t.v1[X], t.v1[Y])))) + std::max(0.0, unitlength*n_xy_e1[0]) + std::max(0.0, unitlength*n_xy_e1[1]);
		double d_xy_e2 = (-1.0 * (n_xy_e2.dot(Eigen::Vector2d(t.v2[X], t.v2[Y])))) + std::max(0.0, unitlength*n_xy_e2[0]) + std::max(0.0, unitlength*n_xy_e2[1]);
		// YZ plane
		Eigen::Vector2d n_yz_e0 = Eigen::Vector2d(-1.0*e0[Z], e0[Y]);
		Eigen::Vector2d n_yz_e1 = Eigen::Vector2d(-1.0*e1[Z], e1[Y]);
		Eigen::Vector2d n_yz_e2 = Eigen::Vector2d(-1.0*e2[Z], e2[Y]);
		if (n[X] < 0.0) {
			n_yz_e0 = -1.0 * n_yz_e0;
			n_yz_e1 = -1.0 * n_yz_e1;
			n_yz_e2 = -1.0 * n_yz_e2;
		}
		double d_yz_e0 = (-1.0 * (n_yz_e0.dot(Eigen::Vector2d(t.v0[Y], t.v0[Z])))) + std::max(0.0, unitlength*n_yz_e0[0]) + std::max(0.0, unitlength*n_yz_e0[1]);
		double d_yz_e1 = (-1.0 * (n_yz_e1.dot(Eigen::Vector2d(t.v1[Y], t.v1[Z])))) + std::max(0.0, unitlength*n_yz_e1[0]) + std::max(0.0, unitlength*n_yz_e1[1]);
		double d_yz_e2 = (-1.0 * (n_yz_e2.dot(Eigen::Vector2d(t.v2[Y], t.v2[Z])))) + std::max(0.0, unitlength*n_yz_e2[0]) + std::max(0.0, unitlength*n_yz_e2[1]);
		// ZX plane
		Eigen::Vector2d n_zx_e0 = Eigen::Vector2d(-1.0*e0[X], e0[Z]);
		Eigen::Vector2d n_zx_e1 = Eigen::Vector2d(-1.0*e1[X], e1[Z]);
		Eigen::Vector2d n_zx_e2 = Eigen::Vector2d(-1.0*e2[X], e2[Z]);
		if (n[Y] < 0.0) {
			n_zx_e0 = -1.0 * n_zx_e0;
			n_zx_e1 = -1.0 * n_zx_e1;
			n_zx_e2 = -1.0 * n_zx_e2;
		}
		double d_xz_e0 = (-1.0 * (n_zx_e0.dot(Eigen::Vector2d(t.v0[Z], t.v0[X])))) + std::max(0.0, unitlength*n_zx_e0[0]) + std::max(0.0, unitlength*n_zx_e0[1]);
		double d_xz_e1 = (-1.0 * (n_zx_e1.dot(Eigen::Vector2d(t.v1[Z], t.v1[X])))) + std::max(0.0, unitlength*n_zx_e1[0]) + std::max(0.0, unitlength*n_zx_e1[1]);
		double d_xz_e2 = (-1.0 * (n_zx_e2.dot(Eigen::Vector2d(t.v2[Z], t.v2[X])))) + std::max(0.0, unitlength*n_zx_e2[0]) + std::max(0.0, unitlength*n_zx_e2[1]);

		// test possible grid boxes for overlap
		for (int x = t_bbox_grid.min[0]; x <= t_bbox_grid.max[0]; x++){
			for (int y = t_bbox_grid.min[1]; y <= t_bbox_grid.max[1]; y++){
				for (int z = t_bbox_grid.min[2]; z <= t_bbox_grid.max[2]; z++){

					uint64_t index = mortonEncode_LUT(z, y, x);

					if (voxels[index - morton_start] == FULL_VOXEL){ continue; } // already marked, continue

					// TRIANGLE PLANE THROUGH BOX TEST
					Vector3 p = Vector3(x*unitlength, y*unitlength, z*unitlength);
					double nDOTp = n.dot(p);
					if ((nDOTp + d1) * (nDOTp + d2) > 0.0){ continue; }

					// PROJECTION TESTS
					// XY
					Eigen::Vector2d p_xy = Eigen::Vector2d(p[X], p[Y]);
					if ((n_xy_e0.dot(p_xy) + d_xy_e0) < 0.0){ continue; }
					if ((n_xy_e1.dot(p_xy) + d_xy_e1) < 0.0){ continue; }
					if ((n_xy_e2.dot(p_xy) + d_xy_e2) < 0.0){ continue; }

					// YZ
					Eigen::Vector2d p_yz = Eigen::Vector2d(p[Y], p[Z]);
					if ((n_yz_e0.dot(p_yz) + d_yz_e0) < 0.0){ continue; }
					if ((n_yz_e1.dot(p_yz) + d_yz_e1) < 0.0){ continue; }
					if ((n_yz_e2.dot(p_yz) + d_yz_e2) < 0.0){ continue; }

					// XZ	
					Eigen::Vector2d p_zx = Eigen::Vector2d(p[Z], p[X]);
					if ((n_zx_e0.dot(p_zx) + d_xz_e0) < 0.0){ continue; }
					if ((n_zx_e1.dot(p_zx) + d_xz_e1) < 0.0){ continue; }
					if ((n_zx_e2.dot(p_zx) + d_xz_e2) < 0.0){ continue; }

					voxels[index - morton_start] = FULL_VOXEL;
					data.push_back(VoxelData(index, fnormals[f], average3Vec(t.v0_color, t.v1_color, t.v2_color))); // we ignore data limits for colored voxelization

					nfilled++;
					continue;
				}
			}
		}
	}
}

inline AABox<Vector3> createMeshBBCube( SurfaceMeshModel * mesh )
{
	Eigen::AlignedBox3d mesh_bbox = mesh->bbox();

	Vector3 mesh_min = mesh_bbox.min();
	Vector3 mesh_max = mesh_bbox.max();

	Vector3 lengths = mesh_max - mesh_min;

	for(int i=0; i<3;i++){
		double delta = lengths.maxCoeff() - lengths[i];
		if(delta != 0){
			mesh_min[i] = mesh_min[i] - (delta / 2.0);
			mesh_max[i] = mesh_max[i] + (delta / 2.0);
		}
	}
	return AABox<Vector3>(mesh_min, mesh_max);
}

inline VoxelContainer ComputeVoxelization( SurfaceMeshModel * mesh, double & unitlength, size_t gridsize = 512, bool isMakeSolid = false )
{
	VoxelContainer container;

	// Move mesh to positive world
	Vector3VertexProperty points = mesh->vertex_coordinates();
	Vector3 corner = mesh->bbox().min();
	Vector3 delta = mesh->bbox().center() - corner;
	for(auto v : mesh->vertices()) points[v] -= corner;
	
	AABox<Vector3> mesh_bbox = createMeshBBCube( mesh );

	unitlength = (mesh_bbox.max[0] - mesh_bbox.min[0]) / (float)gridsize;
	uint64_t morton_part = (gridsize * gridsize * gridsize);

	// Storage for voxel on/off
	std::vector<char> voxels; 

	// morton codes for this partition
	uint64_t start = 0;
	uint64_t end = morton_part;
	size_t nfilled = 0;

	// VOXELIZATION
	voxelize_schwarz_method(mesh, start, end, unitlength, voxels, container.data, nfilled);

	container.translation = corner;

	// Move mesh back to original position
	for(auto v : mesh->vertices()) points[v] += corner;

	if( isMakeSolid )
	{	
		// Original set of surface voxels
		std::vector<char> surface_voxels = voxels; 

		// Flood fill from the outside walls
		std::stack<uint64_t> stack;
		for(int u = 0; u < gridsize; u++){
			for(int v = 0; v < gridsize; v++){
				stack.push( mortonEncode_LUT(u, v, 0) );
				stack.push( mortonEncode_LUT(0, u, v) );
				stack.push( mortonEncode_LUT(v, 0, u) );

				stack.push( mortonEncode_LUT(u, v, (unsigned int)gridsize-1) );
				stack.push( mortonEncode_LUT((unsigned int)gridsize-1, u, v) );
				stack.push( mortonEncode_LUT(v, (unsigned int)gridsize-1, u) );
			}
		}

		// Flood fill (do we need to visit corners as well??)
		while( !stack.empty() )
		{
			uint64_t curVox = stack.top();
			stack.pop();

			if( voxels.at(curVox) == EMPTY_VOXEL )
			{
				voxels.at(curVox) = FULL_VOXEL;

				unsigned int x,y,z;
				mortonDecode(curVox, x, y, z);

				if(x < gridsize - 1) stack.push( mortonEncode_LUT(x+1,y,z) );
				if(y < gridsize - 1) stack.push( mortonEncode_LUT(x,y+1,z) );
				if(z < gridsize - 1) stack.push( mortonEncode_LUT(x,y,z+1) );

				if(x > 0) stack.push( mortonEncode_LUT(x-1,y,z) );
				if(y > 0) stack.push( mortonEncode_LUT(x,y-1,z) );
				if(z > 0) stack.push( mortonEncode_LUT(x,y,z-1) );
			}
		}

		std::vector<VoxelData> inside;

		for(uint64_t m = 0; m < morton_part; m++){
			if(voxels[m] == EMPTY_VOXEL || surface_voxels[m] == FULL_VOXEL)
				inside.push_back( VoxelData(m) );
		}

		container.data = inside;
	}

	return container;
}

enum BooleanOperation{ BOOL_UNION, BOOL_DIFFERENCE, BOOL_INTERSECTION, BOOL_XOR };

inline vector<VoxelData> ComputeVoxelizationCSG( SurfaceMeshModel * meshA, SurfaceMeshModel * meshB, 
										 double & unitlength, size_t gridsize = 1024, BooleanOperation operation = BOOL_UNION )
{	
	// Find bounding box of the two
	Eigen::AlignedBox3d all_bbox = meshA->bbox().merged(meshB->bbox());
	Vector3 corner = all_bbox.min();

	Vector3VertexProperty pointsA = meshA->vertex_coordinates();
	Vector3VertexProperty pointsB = meshB->vertex_coordinates();

	// Move meshes to positive world
	{
		for(auto v : meshA->vertices()) pointsA[v] -= corner;
		for(auto v : meshB->vertices()) pointsB[v] -= corner;
		meshA->updateBoundingBox();
		meshB->updateBoundingBox();
	}

	// Padding - modifies geometry both scaling and translating
	double s = 1.0;
	Vector3 single_voxel;
	{
		all_bbox = meshA->bbox().merged(meshB->bbox());

		unitlength = all_bbox.max().maxCoeff() / gridsize;

		size_t padding = 1;
		single_voxel = Vector3(1,1,1) * (unitlength * padding);

		s = (all_bbox.max().maxCoeff() - (2 * padding * unitlength)) / all_bbox.max().maxCoeff();

		for(auto v : meshA->vertices()) pointsA[v] = (pointsA[v] * s) + single_voxel;
		for(auto v : meshB->vertices()) pointsB[v] = (pointsB[v] * s) + single_voxel;

		meshA->updateBoundingBox();
		meshB->updateBoundingBox();
	}

	// Get size of voxel
	std::vector<char> voxelsA, voxelsB; 
	vector<VoxelData> dataA, dataB; 

	vector<VoxelData> intersection;

	// morton codes
	uint64_t morton_part = (gridsize * gridsize * gridsize);
	uint64_t start = 0;
	uint64_t end = morton_part;

	size_t nfilledA = 0, nfilledB = 0;

	// VOXELIZATION
	voxelize_schwarz_method(meshA, start, end, unitlength, voxelsA, dataA, nfilledA);
	voxelize_schwarz_method(meshB, start, end, unitlength, voxelsB, dataB, nfilledB);

	std::vector< std::vector<char>* > shells;

	shells.push_back( &voxelsA );
	shells.push_back( &voxelsB );

	for(auto shell : shells)
	{
		// Flood fill from outside
		std::stack<VoxelData> stack;
		stack.push( VoxelData( mortonEncode_LUT(0,0,0) ) );

		while( !stack.empty() )
		{
			VoxelData curVox = stack.top();
			stack.pop();

			if( shell->at(curVox.morton) == EMPTY_VOXEL )
			{
				shell->at(curVox.morton) = FULL_VOXEL;

				unsigned int x,y,z;
				mortonDecode(curVox.morton, x, y, z);

				if(x < gridsize - 1) stack.push( VoxelData(mortonEncode_LUT(x+1,y,z)) );
				if(y < gridsize - 1) stack.push( VoxelData(mortonEncode_LUT(x,y+1,z)) );
				if(z < gridsize - 1) stack.push( VoxelData(mortonEncode_LUT(x,y,z+1)) );

				if(x > 0) stack.push( VoxelData(mortonEncode_LUT(x-1,y,z)) );
				if(y > 0) stack.push( VoxelData(mortonEncode_LUT(x,y-1,z)) );
				if(z > 0) stack.push( VoxelData(mortonEncode_LUT(x,y,z-1)) );
			}
		}
	}
	
	for(uint64_t m = 0; m < morton_part; m++)
	{
		switch (operation)
		{
		case BOOL_UNION:
			if(voxelsA[m] == EMPTY_VOXEL || voxelsB[m] == EMPTY_VOXEL)
				intersection.push_back( VoxelData(m) );
			break;
		case BOOL_DIFFERENCE:
			if(voxelsA[m] == EMPTY_VOXEL && voxelsB[m] != EMPTY_VOXEL)
				intersection.push_back( VoxelData(m) );
			break;
		case BOOL_INTERSECTION:
			if(voxelsA[m] == EMPTY_VOXEL && voxelsB[m] == EMPTY_VOXEL)
				intersection.push_back( VoxelData(m) );
			break;
		case BOOL_XOR:
			if((voxelsA[m] == EMPTY_VOXEL || voxelsB[m] == EMPTY_VOXEL) && 
					!(voxelsA[m] == EMPTY_VOXEL && voxelsB[m] == EMPTY_VOXEL))
				intersection.push_back( VoxelData(m) );
			break;
		default:
			break;
		}
	}

	// Restore geometry
	{
		for(auto v : meshA->vertices()) pointsA[v] = ((pointsA[v] * (1.0 / s)) - single_voxel) + corner;
		for(auto v : meshB->vertices()) pointsB[v] = ((pointsB[v] * (1.0 / s)) - single_voxel) + corner;

		meshA->updateBoundingBox();
		meshB->updateBoundingBox();
	}

	return intersection;
}
