// ========================================================================================
//  ApproxMVBB
//  Copyright (C) 2014 by Gabriel Nützi <nuetzig (at) imes (d0t) mavt (d0t) ethz (døt) ch>
//
//  This Source Code Form is subject to the terms of the Mozilla Public
//  License, v. 2.0. If a copy of the MPL was not distributed with this
//  file, You can obtain one at http://mozilla.org/MPL/2.0/.
// ========================================================================================

#pragma once

#include <Eigen/Geometry>
#include <vector>

namespace ApproxMVBB{

	typedef Eigen::Vector3d Vector3;
	typedef Eigen::Vector2d Vector2;
	using  Vector3List = std::vector<Vector3, Eigen::aligned_allocator<Vector3> >;
	using  Vector2List = std::vector<Vector2, Eigen::aligned_allocator<Vector2> >;
	using  Matrix3Dyn = Eigen::Matrix<double, 3, Eigen::Dynamic >;
	using  Matrix2Dyn = Eigen::Matrix<double, 2, Eigen::Dynamic >;
	using  Matrix22 = Eigen::Matrix2d;
	using  Matrix33 = Eigen::Matrix3d;

	using Array3 = Eigen::Array<double, 3, 1>;
	using Array2 = Eigen::Array<double, 2, 1>;

	class OOBB{
	public:

		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

			OOBB() {
			this->reset();
		}

		OOBB(const Vector3 & l,
			const Vector3 & u,
			const Matrix33 & A_IK)
			:m_q_KI(A_IK) {
			m_q_KI.normalize();
			m_minPoint = Vector3(std::min(l(0), u(0)), std::min(l(1), u(1)), std::min(l(2), u(2)));
			m_maxPoint = Vector3(std::max(l(0), u(0)), std::max(l(1), u(1)), std::max(l(2), u(2)));
		}

		/** Switch the z-Axis to the axis with index i (default x-Axis)*/
		void switchZAxis(unsigned int i) {
			if (i > 1) {
				return;
			}
			if (i == 0) {
				// Make new x-Axis the z-Axis
				// R_NK = Rotate around 90∞ around Y, and 90∞ around X (always in K frame) )
				// A_IN = A_IK * A_KN = R_KI * R_NK
				m_q_KI = m_q_KI * Eigen::Quaterniond(0.5, 0.5, 0.5, 0.5);
				// Change points  Coordinates I_[x,y,z] -> K_[y,z,x]
				std::swap(m_minPoint(0), m_minPoint(1));
				std::swap(m_minPoint(1), m_minPoint(2));

				std::swap(m_maxPoint(0), m_maxPoint(1));
				std::swap(m_maxPoint(1), m_maxPoint(2));
			}
			else {
				// Make new y-Axis the z-Axis
				// R_NK = Rotate around 90∞ around -X, and 90∞ around -Y (always in K frame) )
				// A_IN = A_IK * A_KN = R_KI * R_NK
				m_q_KI = m_q_KI * Eigen::Quaterniond(0.5, -0.5, -0.5, -0.5);
				// Change points  Coordinates I_[x,y,z] -> K_[z,x,y]
				std::swap(m_minPoint(0), m_minPoint(2));
				std::swap(m_minPoint(1), m_minPoint(2));

				std::swap(m_maxPoint(0), m_maxPoint(2));
				std::swap(m_maxPoint(1), m_maxPoint(2));
			}

		}


		void expandZeroExtent(double percentageOfLongestAxis, double eps, double defaultExtent){
			//Expand all axes with almost zero extend
			Array3 e = extent();
			Array3::Index idx;
			double max = e.maxCoeff(&idx);

			// if longest axis is also smaller then eps -> make default extent!
			if (max < eps){
				expand(defaultExtent);
				return;
			}
			// otherwise
			double l = 0.5*max*percentageOfLongestAxis;
			for (int i = 0; i < 2; ++i){
				if (i != idx && e(i) < eps){
					m_minPoint(i) -= l;
					m_maxPoint(i) += l;
				}
			}
		}

		void reset() {
			// Violating the constraint min<max for making a completey empty box!
			m_minPoint(0) = std::numeric_limits<double>::max();
			m_maxPoint(0) = std::numeric_limits<double>::min();
			m_minPoint(1) = std::numeric_limits<double>::max();
			m_maxPoint(1) = std::numeric_limits<double>::min();
			m_minPoint(2) = std::numeric_limits<double>::max();
			m_maxPoint(2) = std::numeric_limits<double>::min();
		}

		inline void setZAxisLongest(){
			Vector3::Index i;
			maxExtent(i);
			if (i < 2){
				switchZAxis(i);
			}
		}

		/** Add point expressed in OOBB's K frame to the OOBB */
		template<typename Derived>
		OOBB & unite(const Eigen::MatrixBase<Derived> & p){
			m_maxPoint(0) = std::max(m_maxPoint(0), p(0));
			m_maxPoint(1) = std::max(m_maxPoint(1), p(1));
			m_maxPoint(2) = std::max(m_maxPoint(2), p(2));
			m_minPoint(0) = std::min(m_minPoint(0), p(0));
			m_minPoint(1) = std::min(m_minPoint(1), p(1));
			m_minPoint(2) = std::min(m_minPoint(2), p(2));
			return *this;
		}

		inline Eigen::Array3d extent() const{
			return (m_maxPoint - m_minPoint).array();
		}

		inline double maxExtent() const{
			return (m_maxPoint - m_minPoint).maxCoeff();
		}

		inline double maxExtent(Vector3::Index & i) const{
			return (m_maxPoint - m_minPoint).maxCoeff(&i);
		}

		inline bool isEmpty() const {
			return m_maxPoint(0) <= m_minPoint(0) || m_maxPoint(1) <= m_minPoint(1) || m_maxPoint(2) <= m_minPoint(2);
		}

		inline void expand(double d) {
			//ApproxMVBB_ASSERTMSG(d>=0,"d>=0")
			m_minPoint -= Vector3(d, d, d);
			m_maxPoint += Vector3(d, d, d);
		}

		inline void expand(Vector3 d) {
			//ApproxMVBB_ASSERTMSG(d(0)>=0 && d(1)>=0 && d(2)>=0,"d>=0")
			m_minPoint -= d;
			m_maxPoint += d;
		}

		inline double volume() const {
			Vector3 d = m_maxPoint - m_minPoint;
			return d(0) * d(1) * d(2);
		}

		/** Get direction vectors in I Frame */
		inline Vector3 getDirection(unsigned int i){
			//ApproxMVBB_ASSERTMSG(i<3,"Index wrong: " << i)
			Vector3 d; d.setZero(); d(i) = 1.0;
			return m_q_KI * d; // A_IK* d;
		}

		Eigen::Quaterniond m_q_KI;  ///< Rotation of frame I to frame K, corresponds to a transformation A_IK;
		Vector3 m_minPoint; ///< in K Frame
		Vector3 m_maxPoint; ///< in K Frame
	};

	template<typename Derived>
	OOBB approximateMVBB(const Eigen::MatrixBase<Derived> & points,
						 const double epsilon,
						 const unsigned int pointSamples = 400,
						 const unsigned int gridSize = 5,
						 const unsigned int mvbbDiamOptLoops = 0,
						 const unsigned int gridSearchOptLoops = 6
	);
}
