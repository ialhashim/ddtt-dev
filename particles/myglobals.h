#pragma once

#include <vector>

#include <Eigen/Core>
typedef Eigen::Vector3d Vector3Type;

//    SPHERE_FIBONACCI_POINTS computes sphere points on a Fibonacci spiral.
//
//    John Burkardt - 21 October 2013 / This code is distributed under the GNU LGPL license.
//
//  Reference:
//
//    Richard Swinbank, James Purser,
//    Fibonacci grids: A novel approach to global modelling. July 2006
//
inline std::vector< Vector3Type > sphere_fibonacci_points ( int n = 100 )
{
	double cphi;
	double i_r8,n_r8;
	double phi,sphi,theta;
	const double pi = 3.141592653589793;

	phi = ( 1.0 + std::sqrt ( 5.0 ) ) / 2.0;
	n_r8 = ( double ) ( n );

	std::vector< Vector3Type > points;

	for (int j = 0; j < n; j++ )
	{
		i_r8 = ( double ) ( - n + 1 + 2 * j );
		theta = 2.0 * pi * i_r8 / phi;
		sphi = i_r8 / n_r8;
		cphi = std::sqrt ( ( n_r8 + i_r8 ) * ( n_r8 - i_r8 ) ) / n_r8;

		points.push_back( Vector3Type(cphi * std::sin ( theta ), cphi * std::cos ( theta ), sphi) );
	}

	return points;
}

template<typename Vector3d>
static inline Vector3d orthogonalVector(const Vector3d& n) {
	if ((abs(n.y()) >= 0.9 * abs(n.x())) &&
		abs(n.z()) >= 0.9 * abs(n.x())) return Vector3d(0.0, -n.z(), n.y());
	else if ( abs(n.x()) >= 0.9 * abs(n.y()) &&
		abs(n.z()) >= 0.9 * abs(n.y()) ) return Vector3d(-n.z(), 0.0, n.x());
	else return Vector3d(-n.y(), n.x(), 0.0);
}

template<typename container>
static inline double median(const container & vec){
	typedef container::size_type vec_sz;
	vec_sz size = vec.size();
	sort(vec.begin(), vec.end());
	vec_sz mid = size/2;
	return size % 2 == 0 ? (vec[mid] + vec[mid-1]) / 2 : vec[mid];
}

static inline double deg_to_rad(double deg) {
	double rad = (M_PI * deg) / 180.0;
	return rad;
}

static inline QVector<QColor> rndColors(int count){
	QVector<QColor> c;
	for(int i = 0; i < count; i++) c << starlab::qRandomColor3();
	return c;
}

inline QVector<QColor> rndColors2(int count){
	QVector<QColor> colors;
	float currentHue = 0.0;
	for (int i = 0; i < count; i++){
		colors.push_back( QColor::fromHslF(currentHue, 1.0, 0.5) );
		currentHue += 0.618033988749895f;
		currentHue = std::fmod(currentHue, 1.0f);
	}
	return colors;
}
