#pragma once
#include <fstream>

class Serializable{
public:
	virtual void serialize(std::ostream& os) const = 0;
	virtual void deserialize(std::istream& is) = 0;
};

inline std::ostream& operator<< (std::ostream& os, const Serializable& s) {
	s.serialize(os);
	return os;
}
inline std::istream& operator>> (std::istream& is, Serializable& s) {
	s.deserialize(is);
	return is;
}

// My special types
inline std::ostream& operator<< (std::ostream& os, const Eigen::Vector3d& v) {
	os << v[0] <<" "<< v[1] <<" "<< v[2];
	return os;
}
inline std::istream& operator>> (std::istream& is, Eigen::Vector3d& v) {
	is >> v[0] >> v[1] >> v[2];
	return is;
}

inline std::ostream& operator<< (std::ostream& os, const std::vector<float>& V) {
	os << V.size() << " ";
	for(auto v : V) os << v << " ";
	os << std::endl;
	return os;
}
inline std::istream& operator>> (std::istream& is, std::vector<float>& v) {
	v.clear();
	size_t count; is >> count;
	for(size_t i = 0; i < count; i++){
		float temp;
		is >> temp;
		v.push_back(temp);
	}
	return is;
}
