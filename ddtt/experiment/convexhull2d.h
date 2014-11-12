#pragma once
#include <deque>

template<class Point>
inline bool equal_points(const Point& a, const Point& b, double eps = 1e-16){
	return abs(a[0] - b[0]) < eps && abs(a[1] - b[1]) < eps;
};

// Melkman's algorithm for finding convex hull of 2D points
// From CGAL
template <class Point, class InputIterator, class OutputIterator>
static inline OutputIterator convexhull2d(InputIterator first, InputIterator last, OutputIterator result)
{
	// Subroutine
	auto left_turn = [&](const Point& a, const Point& b, const Point& c){
		return ((b[0] - a[0])*(c[1] - a[1]) - (c[0] - a[0])*(b[1] - a[1]) > 0);
	};

	std::deque<Point> Q;
	if (first == last) return result;								// 0 elements
	Point p = *first;
	if (++first == last) { *result = p; ++result; return result; }	// 1 element
	Point q = *first;
	if (++first == last) {											// 2 elements
		*result = p; ++result;
		if (!equal_points(p, q)){ *result = q; ++result; }
		return result;
	}
	Q.push_back(p);

	Point r;
	while (first != last){
		r = *first;
		// visited input sequence =  p,..., q, r
		if (left_turn(p, q, r)) { Q.push_back(q);  break; }
		if (left_turn(q, p, r)) { Q.push_front(q); break; }
		q = r;
		++first;
	}

	Point current = q;
	if (first != last)           // non-collinear point r exists
	{
		current = r;
		// current, Q.front(), ..., Q.back()
		// ccw convex hull of visited points
		Point s;
		while (++first != last)
		{
			r = *first;
			if (left_turn(current, r, Q.front()) || left_turn(Q.back(), r, current))
			{
				s = current;
				while (!Q.empty() && !left_turn(r, s, Q.front())) { s = Q.front(); Q.pop_front(); }
				Q.push_front(s);
				s = current;
				while (!Q.empty() && !left_turn(s, r, Q.back())) { s = Q.back(); Q.pop_back(); }
				Q.push_back(s);
				current = r;
			}
		}
	}

	Q.push_back(current);
	std::copy(Q.begin(), Q.end(), result);
	return result;
}

template < class Container >
static inline Container convexhull2d(Container & c){
	Container result;
	convexhull2d<Container::value_type>(c.begin(), c.end(), std::back_inserter(result));
	return result;
}

// Warning: only use for small number of points
template < class Container >
static inline std::vector<size_t> convexhull2d_indices(Container & c){
	std::vector<size_t> indices;
	Container result;
	convexhull2d<Container::value_type>(c.begin(), c.end(), std::back_inserter(result));
	auto it = c.begin();
	while (it != c.end()){
		auto & cur = *it;
		for (size_t i = 0; i < result.size(); i++){
			if (equal_points(result[i], cur)){
				indices.push_back(it - c.begin());
				break;
			}
		}
		it++;
	}
	return indices;
}
