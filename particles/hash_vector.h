#pragma once

#include <vector>
#include <numeric>
#include <functional>

namespace std{

	//
	// hash_range
	//
	template<typename InputIterator>
	inline std::size_t hash_range(std::size_t seed, InputIterator first, InputIterator last) 
	{
		return std::accumulate(first, last, seed, 
			[&](size_t sum, int v){
				return sum + (seed ^ (v + 0x9e3779b9 + (seed << 6) + (seed >> 2)));
			}
		);
	}

	template<typename InputIterator>
	inline std::size_t hash_range(InputIterator first, InputIterator last) {
		return std::hash_range(0, first, last);
	}

	template<>
	struct hash< std::vector<int> >{
		std::size_t operator()(std::vector<int> s) const{
			return hash_range(s.begin(), s.end());
		}
	};
}
