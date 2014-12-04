#pragma once
#include <limits>

template<typename Scalar>
struct Bounds
{
	Bounds( const Scalar & minVal = std::numeric_limits<Scalar>::max(), 
		const Scalar & maxVal = std::numeric_limits<Scalar>::lowest() ) 
		: minimum(minVal), maximum(maxVal), count(0), sum(Scalar(0)){}

	Bounds( const Bounds & other ){
		minimum = other.minimum;
		maximum = other.maximum;
		count = other.count;
		sum = other.sum;
	}

	template<typename container>
	static Bounds from( const container & values ){
		Bounds b;
		for(auto & v : values) b.extend(v);
		return b;
	}

	Scalar minimum, maximum;
	Scalar sum;
	size_t count;

	inline void extend(const Scalar & val){
		minimum = std::min(minimum, val);
		maximum = std::max(maximum, val);
		sum += val;
		count++;
	}

	inline Scalar average() const { return sum / count; }
	inline Scalar range() const { return maximum - minimum; }
	inline Scalar intersection( const Bounds & other ) const {
		Scalar overlap = std::max(Scalar(0), 
			std::min(maximum, other.maximum, std::less<Scalar>()) - 
			std::max(minimum, other.minimum, std::greater<Scalar>()));
		return overlap;
	}	
	inline Scalar normalized(const Scalar & val) const {
		return (val - minimum) / (maximum - minimum);
	}
	template<typename container>
	inline container normalize( const container & values ) const {
		container c;
		for(auto & v : values) c.push_back( normalized(v) );
		return c;
	}
};

typedef Bounds<double> Boundsd;
typedef Bounds<float> Boundsf;
