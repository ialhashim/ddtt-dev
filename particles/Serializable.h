#pragma once
#include <QDataStream>

class Serializable{
public:
	virtual void serialize(QDataStream& os) const = 0;
	virtual void deserialize(QDataStream& is) = 0;
};

inline QDataStream& operator<< (QDataStream& os, const Serializable& s) {
	s.serialize(os);
	return os;
}
inline QDataStream& operator>> (QDataStream& is, Serializable& s) {
	s.deserialize(is);
	return is;
}

// My special types
inline QDataStream& operator<< (QDataStream& os, const Eigen::Vector3d& v) {
	os << v[0] << v[1] << v[2];
	return os;
}
inline QDataStream& operator>> (QDataStream& is, Eigen::Vector3d& v) {
	is >> v[0] >> v[1] >> v[2];
	return is;
}

inline QDataStream& operator<< (QDataStream& os, const std::vector<float>& V) {
	os << V.size();
	for(auto v : V) os << v;
	return os;
}
inline QDataStream& operator>> (QDataStream& is, std::vector<float>& v) {
	v.clear();
	size_t count; is >> count;
	for(size_t i = 0; i < count; i++){
		float temp;	is >> temp;
		v.push_back(temp);
	}
	return is;
}
