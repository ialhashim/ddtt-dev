#pragma once
#include "StructureGraph.h"

namespace Structure{
	struct Landmark : public Eigen::Vector3d{
		Landmark(size_t id = -1, const Eigen::Vector3d vec = Eigen::Vector3d(0, 0, 0)) :
			id(id), Eigen::Vector3d(vec), constraint_id(-1){ u = v = -1; partid = "none"; }
		double u, v;
		QString partid;
		size_t id;
		int constraint_id;
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

	typedef QVector<Landmark> Landmarks;

	struct ShapeGraph : public Graph{
		ShapeGraph(QString path) : Graph(path){}
		ShapeGraph(const ShapeGraph& other) : Graph(other){
			this->landmarks = other.landmarks;
		}
		QVector<Landmarks> landmarks;

		void saveLandmarks(QString filename){
			QFile file(filename); file.open(QIODevice::WriteOnly);
			QDataStream os(&file);
			os << (int)landmarks.size();
			for (auto & landmark : landmarks) 
			{
				os << (int)landmark.size();
				for (auto & l : landmark) l.serialize(os);
			}
		}

		void loadLandmarks(QString filename){
			QFile file(filename); file.open(QIODevice::ReadOnly);
			QDataStream is(&file);
			int count;
			is >> count;
			landmarks.resize(count);
			for (int i = 0; i < count; i++) 
			{
				int landmarkCount;
				is >> landmarkCount;
				landmarks[i].resize(landmarkCount);
				for (int j = 0; j < landmarkCount; j++)
					landmarks[i][j].deserialize(is);
			}
		}

		static void correspondTwoCurves(Structure::Curve *sCurve, Structure::Curve *tCurve)
		{
			std::vector<Vector3> sCtrlPoint = sCurve->controlPoints();
			std::vector<Vector3> tCtrlPoint = tCurve->controlPoints();

			// Euclidean for now, could use Geodesic distance instead if need
			Vector3 scenter = sCurve->center();
			Vector3 sfront = sCtrlPoint.front() - scenter;
			Vector3 tcenter = tCurve->center();
			Vector3 tfront = tCtrlPoint.front() - tcenter;
			Vector3 tback = tCtrlPoint.back() - tcenter;

			float f2f = (sfront - tfront).norm();
			float f2b = (sfront - tback).norm();

			float diff = std::abs(f2f - f2b);

			float threshold = 0.1f;

			if (f2f > f2b && diff > threshold)
			{
				// Flip the target
				std::vector<Scalar> tCtrlWeight = tCurve->controlWeights();
				std::reverse(tCtrlPoint.begin(), tCtrlPoint.end());
				std::reverse(tCtrlWeight.begin(), tCtrlWeight.end());

				NURBS::NURBSCurved newCurve(tCtrlPoint, tCtrlWeight);
				tCurve->curve = newCurve;
			}
		}

		static void correspondTwoSheets(Structure::Sheet *sSheet, Structure::Sheet *tSheet)
		{
			// Old properties
			NURBS::NURBSRectangled &oldRect = tSheet->surface;
			int uDegree = oldRect.GetDegree(0);
			int vDegree = oldRect.GetDegree(1);
			bool uLoop = oldRect.IsLoop(0);
			bool vLoop = oldRect.IsLoop(1);
			bool uOpen = oldRect.IsOpen(0);
			bool vOpen = oldRect.IsOpen(1);
			bool isModified = false;
			bool isUVFlipped = false;

			// Control points and weights
			Array2D_Vector3 sCtrlPoint = sSheet->surface.mCtrlPoint;
			Array2D_Real sCtrlWeight = sSheet->surface.mCtrlWeight;

			Array2D_Vector3 tCtrlPoint = tSheet->surface.mCtrlPoint;
			Array2D_Real tCtrlWeight = tSheet->surface.mCtrlWeight;

			Array2D_Vector3 tCtrlPointNew;
			Array2D_Real tCtrlWeightNew;

			Vector3 scenter = sSheet->center();
			Vector3 tcenter = tSheet->center();

			// Get the extreme points.
			Vector3 s00 = sCtrlPoint.front().front();
			Vector3 s01 = sCtrlPoint.front().back();
			Vector3 s10 = sCtrlPoint.back().front();
			Vector3 sU = s10 - s00;
			Vector3 sV = s01 - s00;

			Vector3 t00 = tCtrlPoint.front().front();
			Vector3 t01 = tCtrlPoint.front().back();
			Vector3 t10 = tCtrlPoint.back().front();
			Vector3 tU = t10 - t00;
			Vector3 tV = t01 - t00;

			// Flip if need
			Vector3 sUV = cross(sU, sV);
			Vector3 tUV = cross(tU, tV);
			if (dot(sUV, tUV) < 0)
			{
				// Reverse the target along u direction
				std::reverse(tCtrlPoint.begin(), tCtrlPoint.end());
				std::reverse(tCtrlWeight.begin(), tCtrlWeight.end());

				// Update tU
				tU = -tU;
				tUV = -tUV;
				isModified = true;
			}

			// Rotate if need
			Scalar cosAngle = dot(sU.normalized(), tU.normalized());
			Scalar cos45 = sqrtf(2.0) / 2;

			// Do Nothing
			if (cosAngle > cos45)
			{
				tCtrlPointNew = tCtrlPoint;
				tCtrlWeightNew = tCtrlWeight;
			}
			// Rotate 180 degrees
			else if (cosAngle < -cos45)
			{
				//  --> sV				tU
				// |					|
				// sU             tV <--
				// By flipping along both directions
				std::reverse(tCtrlPoint.begin(), tCtrlPoint.end());
				std::reverse(tCtrlWeight.begin(), tCtrlWeight.end());

				for (int i = 0; i < (int)tCtrlPoint.size(); i++)
				{
					std::reverse(tCtrlPoint[i].begin(), tCtrlPoint[i].end());
					std::reverse(tCtrlWeight[i].begin(), tCtrlWeight[i].end());
				}

				// The new control points and weights
				tCtrlPointNew = tCtrlPoint;
				tCtrlWeightNew = tCtrlWeight;
				isModified = true;
			}
			// Rotate 90 degrees 
			else
			{
				Vector3 stU = cross(sU, tU);
				if (dot(stU, sUV) >= 0)
				{
					//  --> sV		tV
					// |			|
					// sU           --> tU
					// Transpose and reverse along U
					tCtrlPointNew = transpose<Vector3>(tCtrlPoint);
					tCtrlWeightNew = transpose<Scalar>(tCtrlWeight);

					std::reverse(tCtrlPointNew.begin(), tCtrlPointNew.end());
					std::reverse(tCtrlWeightNew.begin(), tCtrlWeightNew.end());
				}
				else
				{
					//  --> sV	tU<--
					// |			 |
					// sU			tV
					// Reverse along U and Transpose
					std::reverse(tCtrlPoint.begin(), tCtrlPoint.end());
					std::reverse(tCtrlWeight.begin(), tCtrlWeight.end());

					tCtrlPointNew = transpose<Vector3>(tCtrlPoint);
					tCtrlWeightNew = transpose<Scalar>(tCtrlWeight);
				}

				isModified = true;
				isUVFlipped = true;
			}

			// Create a new sheet if need
			if (isModified)
			{
				NURBS::NURBSRectangled newRect;
				if (isUVFlipped)
					newRect = NURBS::NURBSRectangled(tCtrlPointNew, tCtrlWeightNew, vDegree, uDegree, vLoop, uLoop, vOpen, uOpen);
				else
					newRect = NURBS::NURBSRectangled(tCtrlPointNew, tCtrlWeightNew, uDegree, vDegree, uLoop, vLoop, uOpen, vOpen);

				tSheet->surface = newRect;
			}
		}
	};
}
