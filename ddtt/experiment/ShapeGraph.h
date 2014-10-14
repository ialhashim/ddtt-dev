#pragma once
#include "StructureGraph.h"

namespace Structure{
	struct Landmark : public Eigen::Vector3d{
		Landmark(size_t id = -1, const Eigen::Vector3d vec = Eigen::Vector3d(0, 0, 0)) :
			id(id), Eigen::Vector3d(vec){ u = v = -1; partid = "none"; }
		double u, v;
		QString partid;
		size_t id;
		void serialize(QDataStream& os) const {
			os << id;
			os << partid;
			os << u << v;
			os << this->x() << this->y() << this->z();
		}
		void deserialize(QDataStream& is){
			is >> id;
			is >> partid;
			is >> u >> v;
			is >> this->x() >> this->y() >> this->z();
		}
	};
	struct ShapeGraph : public Graph{
		ShapeGraph(QString path) :Graph(path){}
		QVector<Landmark> landmarks;

		void saveLandmarks(QString filename){
			QFile file(filename); file.open(QIODevice::WriteOnly);
			QDataStream os(&file);
			os << (int)landmarks.size();
			for (auto & l : landmarks) l.serialize(os);
		}

		void loadLandmarks(QString filename){
			QFile file(filename); file.open(QIODevice::ReadOnly);
			QDataStream is(&file);
			int count;
			is >> count;
			landmarks.resize(count);
			for (int i = 0; i < count; i++) landmarks[i].deserialize(is);
		}
	};
}
