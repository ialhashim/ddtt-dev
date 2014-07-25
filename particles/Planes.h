#pragma once

#include <vector>
#include <QMap>
#include <QSet>

#include <Eigen/Core>

struct Plane{ 
	typedef Eigen::Vector3d Vector3;

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
