#pragma once
#include <unordered_map>
#include <vector>

template<typename Vector3, typename Scalar>
struct SpatialHash
{
	SpatialHash( const std::vector<Vector3> & particles, Scalar cellsize = Scalar(1) / 128, size_t tablesize = 1e6 ) 
		: cellsize(cellsize), tablesize(tablesize)
	{
		Scalar maxScalarVal = std::numeric_limits<Scalar>::max();
		minbound = Vector3( maxScalarVal, maxScalarVal, maxScalarVal );
		maxbound = -minbound;

		// Find particles bounds
		for(size_t i = 0; i < particles.size(); i++){
			for(int j = 0; j < 3; j++){
				minbound[j] = std::min(minbound[j], particles[i][j]);
				maxbound[j] = std::max(maxbound[j], particles[i][j]);
			}
		}

		// Add particles to grid
		for(size_t i = 0; i < particles.size(); i++){
			int c[3];
			for(int j = 0; j < 3; j++) c[j] = (particles[i][j] - minbound[j]) / cellsize;

			auto key = hash(c[0],c[1],c[2]);
			auto got = grid.find( key );

			if(got == grid.end()){
				std::vector<size_t> list;
				list.push_back( i );
				grid.insert( std::make_pair(key, list) );
			}
			else
				grid[ key ].push_back( i );
		}
	}

	inline std::vector<size_t> intersection( const Vector3 & p ) const{
		int c[3];
		for(int j = 0; j < 3; j++){
			if(p[j] < minbound[j] || p[j] > maxbound[j]) return std::vector<size_t>();
			c[j] = (p[j] - minbound[j]) / cellsize;
		}
		auto key = hash(c[0],c[1],c[2]);
		if(grid.find( key ) == grid.end()) return std::vector<size_t>();
		else return grid.at( key );
	}

	std::vector<size_t> around( const Vector3 & p, Scalar radius ) const
	{
		std::vector<size_t> isect;

		int c[3];
		for(int j = 0; j < 3; j++) c[j] = (p[j] - minbound[j]) / cellsize;
		
		std::vector<size_t> sphere;
		int r = radius / cellsize;
		if(r < 1) return isect;

		radius = pow(radius, 2);

		for(int z = -r; z <= r; z++){
			for(int y = -r; y <= r; y++){
				for(int x = -r; x <= r; x++){
					int cj[3] = { c[0] + x, c[1] + y, c[2] + z };
					if(cj[0] < 0 || cj[1] < 0 || cj[2] < 0) continue;

					Vector3 pj( (cj[0] * cellsize) + (minbound[0]), 
								(cj[1] * cellsize) + (minbound[1]), 
								(cj[2] * cellsize) + (minbound[2]) );

					if( (p - pj).squaredNorm() > radius ) continue;

					auto it = intersection( pj );
					isect.insert( isect.end(), it.begin(), it.end() );
				}
			}
		}

		return isect;
	}

	inline size_t hash(int ix, int iy, int iz) const
	{ return size_t ((ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791)) % tablesize; }

	size_t tablesize;
	Scalar cellsize;
	Vector3 minbound, maxbound;
	std::unordered_map< size_t, std::vector<size_t> > grid;
};
