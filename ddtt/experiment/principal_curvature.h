// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PRINCIPAL_CURVATURE_H
#define IGL_PRINCIPAL_CURVATURE_H

#include <Eigen/SparseCholesky>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <vector>
#include <stdio.h>
#include <map>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <list>
#include <cmath>
#include <limits>

// Lib IGL includes
//#include <igl/igl_inline.h>
//#include <igl/cotmatrix.h>
//#include <igl/writeOFF.h>

//#include <igl/adjacency_list.h>
//#include <igl/per_face_normals.h>
//#include <igl/per_vertex_normals.h>
//#include <igl/avg_edge_length.h>
//#include <igl/vertex_triangle_adjacency.h>

#ifndef IGL_STATIC_LIBRARY
#  define IGL_INLINE inline
#else
#  define IGL_INLINE
#endif

typedef enum
{
	SPHERE_SEARCH,
	K_RING_SEARCH
} searchType;

typedef enum
{
	AVERAGE,
	PROJ_PLANE
} normalType;


namespace igl
{
	// Utility:
	template <typename FaceArray, typename Index>
	inline void adjacency_list(const FaceArray & F, std::vector<std::vector<Index> >& A)
	{
		A.clear();
		A.resize(F.maxCoeff() + 1);

		// Loop over faces
		for (int i = 0; i < F.rows(); i++)
		{
			// Loop over this face (triangle)
			for (int j = 0; j < 3; j++)
			{
				// Get indices of edge: s --> d
				int s = F(i, j);
				int d = F(i, (j + 1) % 3);
				A.at(s).push_back(d);
				A.at(d).push_back(s);
			}
		}

		// Remove duplicates
		for (int i = 0; i < (int)A.size(); ++i)
		{
			std::sort(A[i].begin(), A[i].end());
			A[i].erase(std::unique(A[i].begin(), A[i].end()), A[i].end());
		}
	}

	template <typename DerivedV, typename DerivedF, typename IndexType>
	inline void vf(const Eigen::PlainObjectBase<DerivedV>& V, const Eigen::PlainObjectBase<DerivedF>& F,
		std::vector<std::vector<IndexType> >& VF, std::vector<std::vector<IndexType> >& VFi)
	{
		VF.clear(); VF.resize(V.rows());
		VFi.clear(); VFi.resize(V.rows());

		for (int fi = 0; fi < F.rows(); ++fi){
			for (int i = 0; i < F.cols(); ++i){
				VF[F(fi, i)].push_back(fi);
				VFi[F(fi, i)].push_back(i);
			}
		}
	}

	template <typename DerivedV, typename DerivedF, typename IndexType>
	inline void vertex_triangle_adjacency(
		const Eigen::PlainObjectBase<DerivedV>& V,
		const Eigen::PlainObjectBase<DerivedF>& F,
		std::vector<std::vector<IndexType> >& VF,
		std::vector<std::vector<IndexType> >& VFi)
	{
		VF.clear();
		VFi.clear();

		VF.resize(V.rows());
		VFi.resize(V.rows());

		for (int fi = 0; fi < F.rows(); ++fi)
		{
			for (int i = 0; i < F.cols(); ++i)
			{
				VF[F(fi, i)].push_back(fi);
				VFi[F(fi, i)].push_back(i);
			}
		}
	}

#define SQRT_ONE_OVER_THREE 0.57735026918962573
	template <typename DerivedV, typename DerivedF, typename DerivedZ, typename DerivedN>
	inline void per_face_normals(
		const Eigen::PlainObjectBase<DerivedV>& V,
		const Eigen::PlainObjectBase<DerivedF>& F,
		const Eigen::PlainObjectBase<DerivedZ> & Z,
		Eigen::PlainObjectBase<DerivedN> & N)
	{
		N.resize(F.rows(), 3);
		// loop over faces
		int Frows = F.rows();
		//#pragma omp parallel for if (Frows>10000)
		for (int i = 0; i < Frows; i++)
		{
			const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> v1 = V.row(F(i, 1)) - V.row(F(i, 0));
			const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> v2 = V.row(F(i, 2)) - V.row(F(i, 0));
			N.row(i) = v1.cross(v2);//.normalized();
			typename DerivedV::Scalar r = N.row(i).norm();
			if (r == 0)
			{
				N.row(i) = Z;
			}
			else
			{
				N.row(i) /= r;
			}
		}
	}

	template <typename DerivedV, typename DerivedF, typename DerivedN>
	inline void per_face_normals(
		const Eigen::PlainObjectBase<DerivedV>& V,
		const Eigen::PlainObjectBase<DerivedF>& F,
		Eigen::PlainObjectBase<DerivedN> & N)
	{
		using namespace Eigen;
		Matrix<typename DerivedN::Scalar, 3, 1> Z(0, 0, 0);
		return per_face_normals(V, F, Z, N);
	}

	template <typename DerivedV, typename DerivedF, typename DeriveddblA>
	inline void doublearea(
		const Eigen::PlainObjectBase<DerivedV> & V,
		const Eigen::PlainObjectBase<DerivedF> & F,
		Eigen::PlainObjectBase<DeriveddblA> & dblA)
	{
		const int dim = V.cols();
		// Only support triangles
		assert(F.cols() == 3);
		const int m = F.rows();
		// Compute edge lengths
		Eigen::PlainObjectBase<DerivedV> l;
		// "Lecture Notes on Geometric Robustness" Shewchuck 09, Section 3.1
		// http://www.cs.berkeley.edu/~jrs/meshpapers/robnotes.pdf

		// Projected area helper
		const auto & proj_doublearea =
			[&V, &F](const int x, const int y, const int f)->double
		{
			auto rx = V(F(f, 0), x) - V(F(f, 2), x);
			auto sx = V(F(f, 1), x) - V(F(f, 2), x);
			auto ry = V(F(f, 0), y) - V(F(f, 2), y);
			auto sy = V(F(f, 1), y) - V(F(f, 2), y);
			return rx*sy - ry*sx;
		};

		switch (dim)
		{
		case 3:
		{
			dblA = Eigen::PlainObjectBase<DeriveddblA>::Zero(m, 1);
			for (int f = 0; f<F.rows(); f++)
			{
				for (int d = 0; d<3; d++)
				{
					double dblAd = proj_doublearea(d, (d + 1) % 3, f);
					dblA(f) += dblAd*dblAd;
				}
			}
			dblA = dblA.array().sqrt().eval();
			break;
		}
		case 2:
		{
			dblA.resize(m, 1);
			for (int f = 0; f<m; f++)
			{
				dblA(f) = proj_doublearea(0, 1, f);
			}
			break;
		}
		/*default:
		{
			edge_lengths(V, F, l);
			return doublearea(l, dblA);
		}*/
		}
	}

	template <
		typename DerivedA,
		typename DerivedB,
		typename DerivedC,
		typename DerivedD>
		inline void doublearea(
		const Eigen::PlainObjectBase<DerivedA> & A,
		const Eigen::PlainObjectBase<DerivedB> & B,
		const Eigen::PlainObjectBase<DerivedC> & C,
		Eigen::PlainObjectBase<DerivedD> & D)
	{
		assert(A.cols() == 2 && "corners should be 2d");
		assert(B.cols() == 2 && "corners should be 2d");
		assert(C.cols() == 2 && "corners should be 2d");
		assert(A.rows() == B.rows() && "corners should have same length");
		assert(A.rows() == C.rows() && "corners should have same length");
		const auto & R = A - C;
		const auto & S = B - C;
		D = R.col(0).array()*S.col(1).array() - R.col(1).array()*S.col(0).array();
	}

	template <
		typename DerivedA,
		typename DerivedB,
		typename DerivedC>
		inline typename DerivedA::Scalar doublearea_single(
		const Eigen::PlainObjectBase<DerivedA> & A,
		const Eigen::PlainObjectBase<DerivedB> & B,
		const Eigen::PlainObjectBase<DerivedC> & C)
	{
		auto r = A - C;
		auto s = B - C;
		return r(0)*s(1) - r(1)*s(0);
	}

	template <typename Derivedl, typename DeriveddblA>
	inline void doublearea(
		const Eigen::PlainObjectBase<Derivedl> & ul,
		Eigen::PlainObjectBase<DeriveddblA> & dblA)
	{
		using namespace Eigen;
		using namespace std;
		// Only support triangles
		assert(ul.cols() == 3);
		// Number of triangles
		const int m = ul.rows();
		Eigen::PlainObjectBase<Derivedl> l;
		MatrixXi _;
		sort(ul, 2, false, l, _);
		// semiperimeters
		Matrix<typename Derivedl::Scalar, Dynamic, 1> s = l.rowwise().sum()*0.5;
		assert(s.rows() == m);
		// resize output
		dblA.resize(l.rows(), 1);
		for (int i = 0; i<m; i++)
		{
			//// Heron's formula for area
			//const typename Derivedl::Scalar arg =
			//  s(i)*(s(i)-l(i,0))*(s(i)-l(i,1))*(s(i)-l(i,2));
			//assert(arg>=0);
			//dblA(i) = 2.0*sqrt(arg);
			// Kahan's Heron's formula
			const typename Derivedl::Scalar arg =
				(l(i, 0) + (l(i, 1) + l(i, 2)))*
				(l(i, 2) - (l(i, 0) - l(i, 1)))*
				(l(i, 2) + (l(i, 0) - l(i, 1)))*
				(l(i, 0) + (l(i, 1) - l(i, 2)));
			dblA(i) = 2.0*0.25*sqrt(arg);
			assert(l(i, 2) - (l(i, 0) - l(i, 1)) && "FAILED KAHAN'S ASSERTION");
			assert(dblA(i) == dblA(i) && "DOUBLEAREA() PRODUCED NaN");
		}
	}

	enum PerVertexNormalsWeightingType
	{
		// Incident face normals have uniform influence on vertex normal
		PER_VERTEX_NORMALS_WEIGHTING_TYPE_UNIFORM = 0,
		// Incident face normals are averaged weighted by area
		PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA = 1,
		// Incident face normals are averaged weighted by incident angle of vertex
		PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE = 2,
		// Area weights
		PER_VERTEX_NORMALS_WEIGHTING_TYPE_DEFAULT = 3,
		NUM_PER_VERTEX_NORMALS_WEIGHTING_TYPE = 4
	};

	template <typename DerivedV, typename DerivedF>
	inline void per_vertex_normals(
		const Eigen::PlainObjectBase<DerivedV>& V,
		const Eigen::PlainObjectBase<DerivedF>& F,
		const igl::PerVertexNormalsWeightingType weighting,
		const Eigen::PlainObjectBase<DerivedV>& FN,
		Eigen::PlainObjectBase<DerivedV> & N)
	{
		// Resize for output
		N.setZero(V.rows(), 3);

		Eigen::MatrixXd W(F.rows(), 3);
		switch (weighting)
		{
		case PER_VERTEX_NORMALS_WEIGHTING_TYPE_UNIFORM:
			W.setConstant(1.);
			break;
		default:
			assert(false && "Unknown weighting type");
		case PER_VERTEX_NORMALS_WEIGHTING_TYPE_DEFAULT:
		case PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA:
		{
			Eigen::VectorXd A;
			doublearea(V, F, A);
			W = A.replicate(1, 3);
			break;
		}
		case PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE:
			//internal_angles(V, F, W);
			break;
		}

		// loop over faces
		const int Frows = F.rows();
		//// Minimum number of iterms per openmp thread
		//#ifndef IGL_OMP_MIN_VALUE
		//#  define IGL_OMP_MIN_VALUE 1000
		//#endif
		//#pragma omp parallel for if (Frows>IGL_OMP_MIN_VALUE)
		for (int i = 0; i < Frows; i++)
		{
			// throw normal at each corner
			for (int j = 0; j < 3; j++)
			{
				// Does this need to be critical?
				//#pragma omp critical
				N.row(F(i, j)) += W(i, j)*FN.row(i);
			}
		}
		// take average via normalization
		N.rowwise().normalize();
	}

	template <typename DerivedV, typename DerivedF>
	inline void per_vertex_normals(
		const Eigen::PlainObjectBase<DerivedV>& V,
		const Eigen::PlainObjectBase<DerivedF>& F,
		const Eigen::PlainObjectBase<DerivedV>& FN,
		Eigen::PlainObjectBase<DerivedV> & N)
	{
		return per_vertex_normals(V, F, PER_VERTEX_NORMALS_WEIGHTING_TYPE_DEFAULT, FN, N);
	}


	// Compute the principal curvature directions and magnitude of the given triangle mesh
	//   DerivedV derived from vertex positions matrix type: i.e. MatrixXd
	//   DerivedF derived from face indices matrix type: i.e. MatrixXi
	// Inputs:
	//   V       eigen matrix #V by 3
	//   F       #F by 3 list of mesh faces (must be triangles)
	//   radius  controls the size of the neighbourhood used, 1 = average edge lenght
	//
	// Outputs:
	//   PD1 #V by 3 maximal curvature direction for each vertex.
	//   PD2 #V by 3 minimal curvature direction for each vertex.
	//   PV1 #V by 1 maximal curvature value for each vertex.
	//   PV2 #V by 1 minimal curvature value for each vertex.
	//
	// See also: average_onto_faces, average_onto_vertices
	//
	// This function has been developed by: Nikolas De Giorgis, Luigi Rocca and Enrico Puppo.
	// The algorithm is based on:
	// Efficient Multi-scale Curvature and Crease Estimation
	// Daniele Panozzo, Enrico Puppo, Luigi Rocca
	// GraVisMa, 2010

	template <
		typename DerivedV,
		typename DerivedF,
		typename DerivedPD1,
		typename DerivedPD2,
		typename DerivedPV1,
		typename DerivedPV2>
		IGL_INLINE void principal_curvature(
		const Eigen::PlainObjectBase<DerivedV>& V,
		const Eigen::PlainObjectBase<DerivedF>& F,
		Eigen::PlainObjectBase<DerivedPD1>& PD1,
		Eigen::PlainObjectBase<DerivedPD2>& PD2,
		Eigen::PlainObjectBase<DerivedPV1>& PV1,
		Eigen::PlainObjectBase<DerivedPV2>& PV2,
		unsigned radius = 5,
		bool useKring = true);
}

class CurvatureCalculator
{
public:
	/* Row number i represents the i-th vertex, whose columns are:
	 curv[i][0] : K2
	 curv[i][1] : K1
	 curvDir[i][0] : PD1
	 curvDir[i][1] : PD2
	 */
	std::vector< std::vector<double> > curv;
	std::vector< std::vector<Eigen::Vector3d> > curvDir;
	bool curvatureComputed;
	class Quadric
	{
	public:

		IGL_INLINE Quadric()
		{
			a() = b() = c() = d() = e() = 1.0;
		}

		IGL_INLINE Quadric(double av, double bv, double cv, double dv, double ev)
		{
			a() = av;
			b() = bv;
			c() = cv;
			d() = dv;
			e() = ev;
		}

		IGL_INLINE double& a() { return data[0]; }
		IGL_INLINE double& b() { return data[1]; }
		IGL_INLINE double& c() { return data[2]; }
		IGL_INLINE double& d() { return data[3]; }
		IGL_INLINE double& e() { return data[4]; }

		double data[5];

		IGL_INLINE double evaluate(double u, double v)
		{
			return a()*u*u + b()*u*v + c()*v*v + d()*u + e()*v;
		}

		IGL_INLINE double du(double u, double v)
		{
			return 2.0*a()*u + b()*v + d();
		}

		IGL_INLINE double dv(double u, double v)
		{
			return 2.0*c()*v + b()*u + e();
		}

		IGL_INLINE double duv(double u, double v)
		{
			return b();
		}

		IGL_INLINE double duu(double u, double v)
		{
			return 2.0*a();
		}

		IGL_INLINE double dvv(double u, double v)
		{
			return 2.0*c();
		}


		IGL_INLINE static Quadric fit(std::vector<Eigen::Vector3d> &VV, bool zeroDetCheck, bool svd)
		{
			using namespace std;
			assert(VV.size() >= 5);
			if (VV.size() < 5)
			{
				cerr << "ASSERT FAILED!" << endl;
				exit(0);
			}

			Eigen::MatrixXd A(VV.size(), 5);
			Eigen::MatrixXd b(VV.size(), 1);
			Eigen::MatrixXd sol(5, 1);

			for (unsigned int c = 0; c < VV.size(); ++c)
			{
				double u = VV[c][0];
				double v = VV[c][1];
				double n = VV[c][2];

				A(c, 0) = u*u;
				A(c, 1) = u*v;
				A(c, 2) = v*v;
				A(c, 3) = u;
				A(c, 4) = v;

				b(c) = n;
			}

			sol = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

			return Quadric(sol(0), sol(1), sol(2), sol(3), sol(4));
		}
	};

public:

	Eigen::MatrixXd vertices;
	// Face list of current mesh    (#F x 3) or (#F x 4)
	// The i-th row contains the indices of the vertices that forms the i-th face in ccw order
	Eigen::MatrixXi faces;

	std::vector<std::vector<int> > vertex_to_vertices;
	std::vector<std::vector<int> > vertex_to_faces;
	std::vector<std::vector<int> > vertex_to_faces_index;
	Eigen::MatrixXd face_normals;
	Eigen::MatrixXd vertex_normals;

	/* Size of the neighborhood */
	double sphereRadius;
	int kRing;

	bool localMode; /* Use local mode */
	bool projectionPlaneCheck; /* Check collected vertices on tangent plane */
	bool montecarlo;
	bool svd; /* Use svd calculation instead of pseudoinverse */
	bool zeroDetCheck; /* Check if the determinant is close to zero */
	unsigned int montecarloN;

	searchType st; /* Use either a sphere search or a k-ring search */
	normalType nt;

	double lastRadius;
	double scaledRadius;
	std::string lastMeshName;

	/* Benchmark related variables */
	bool expStep; /* True if we want the radius to increase exponentially */
	int step;  /* If expStep==false, by how much rhe radius increases on every step */
	int maxSize; /* The maximum limit of the radius in the benchmark */

	IGL_INLINE CurvatureCalculator();
	IGL_INLINE void init(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

	IGL_INLINE void finalEigenStuff(int, std::vector<Eigen::Vector3d>, Quadric);
	IGL_INLINE void fitQuadric(Eigen::Vector3d, std::vector<Eigen::Vector3d> ref, const  std::vector<int>&, Quadric *);
	IGL_INLINE void applyProjOnPlane(Eigen::Vector3d, std::vector<int>, std::vector<int>&);
	IGL_INLINE void getSphere(const int, const double, std::vector<int>&, int min);
	IGL_INLINE void getKRing(const int, const double, std::vector<int>&);
	IGL_INLINE Eigen::Vector3d project(Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d);
	IGL_INLINE void computeReferenceFrame(int, Eigen::Vector3d, std::vector<Eigen::Vector3d>&);
	IGL_INLINE void getAverageNormal(int, std::vector<int>, Eigen::Vector3d&);
	IGL_INLINE void getProjPlane(int, std::vector<int>, Eigen::Vector3d&);
	IGL_INLINE void applyMontecarlo(std::vector<int>&, std::vector<int>*);
	IGL_INLINE void computeCurvature();
	IGL_INLINE void printCurvature(std::string outpath);
	IGL_INLINE double getAverageEdge();

	IGL_INLINE static int rotateForward(float *v0, float *v1, float *v2)
	{
		float t;

		if (abs(*v2) >= abs(*v1) && abs(*v2) >= abs(*v0))
			return 0;

		t = *v0;
		*v0 = *v2;
		*v2 = *v1;
		*v1 = t;

		return 1 + rotateForward(v0, v1, v2);
	}

	IGL_INLINE static void rotateBackward(int nr, float *v0, float *v1, float *v2)
	{
		float t;

		if (nr == 0)
			return;

		t = *v2;
		*v2 = *v0;
		*v0 = *v1;
		*v1 = t;

		rotateBackward(nr - 1, v0, v1, v2);
	}

	IGL_INLINE static Eigen::Vector3d chooseMax(Eigen::Vector3d n, Eigen::Vector3d abc, float ab)
	{
		int i, max_i;
		float max_sp;
		Eigen::Vector3d nt[8];

		n.normalize();
		abc.normalize();

		max_sp = -std::numeric_limits<float>::max();

		for (i = 0; i < 4; i++)
		{
			nt[i] = n;
			if (ab > 0)
			{
				switch (i)
				{
				case 0:
					break;

				case 1:
					nt[i][2] = -n[2];
					break;

				case 2:
					nt[i][0] = -n[0];
					nt[i][1] = -n[1];
					break;

				case 3:
					nt[i][0] = -n[0];
					nt[i][1] = -n[1];
					nt[i][2] = -n[2];
					break;
				}
			}
			else
			{
				switch (i)
				{
				case 0:
					nt[i][0] = -n[0];
					break;

				case 1:
					nt[i][1] = -n[1];
					break;

				case 2:
					nt[i][0] = -n[0];
					nt[i][2] = -n[2];
					break;

				case 3:
					nt[i][1] = -n[1];
					nt[i][2] = -n[2];
					break;
				}
			}

			if (nt[i].dot(abc) > max_sp)
			{
				max_sp = nt[i].dot(abc);
				max_i = i;
			}
		}

		return nt[max_i];
	}

};

class comparer
{
public:
	IGL_INLINE bool operator() (const std::pair<int, double>& lhs, const std::pair<int, double>&rhs) const
	{
		return lhs.second > rhs.second;
	}
};

IGL_INLINE CurvatureCalculator::CurvatureCalculator()
{
	this->localMode = true;
	this->projectionPlaneCheck = true;
	this->sphereRadius = 5;
	this->st = SPHERE_SEARCH;
	this->nt = AVERAGE;
	this->montecarlo = false;
	this->montecarloN = 0;
	this->kRing = 3;
	this->svd = true;
	this->zeroDetCheck = true;
	this->curvatureComputed = false;
	this->expStep = true;
}

IGL_INLINE void CurvatureCalculator::init(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
	// Normalize vertices
	vertices = V;

	//  vertices = vertices.array() - vertices.minCoeff();
	//  vertices = vertices.array() / vertices.maxCoeff();
	//  vertices = vertices.array() * (1.0/igl::avg_edge_length(V,F));

	faces = F;
	igl::adjacency_list(F, vertex_to_vertices);
	igl::vertex_triangle_adjacency(V, F, vertex_to_faces, vertex_to_faces_index);
	igl::per_face_normals(V, F, face_normals);
	igl::per_vertex_normals(V, F, face_normals, vertex_normals);
}

IGL_INLINE void CurvatureCalculator::fitQuadric(Eigen::Vector3d v, std::vector<Eigen::Vector3d> ref, const std::vector<int>& vv, Quadric *q)
{
	std::vector<Eigen::Vector3d> points;
	points.reserve(vv.size());

	for (unsigned int i = 0; i < vv.size(); ++i) {

		Eigen::Vector3d  cp = vertices.row(vv[i]);

		// vtang non e` il v tangente!!!
		Eigen::Vector3d  vTang = cp - v;

		double x = vTang.dot(ref[0]);
		double y = vTang.dot(ref[1]);
		double z = vTang.dot(ref[2]);
		points.push_back(Eigen::Vector3d(x, y, z));
	}
	*q = Quadric::fit(points, zeroDetCheck, svd);
}

IGL_INLINE void CurvatureCalculator::finalEigenStuff(int i, std::vector<Eigen::Vector3d> ref, Quadric q)
{

	double a = q.a();
	double b = q.b();
	double c = q.c();
	double d = q.d();
	double e = q.e();

	//  if (fabs(a) < 10e-8 || fabs(b) < 10e-8)
	//  {
	//    std::cout << "Degenerate quadric: " << i << std::endl;
	//  }

	double E = 1.0 + d*d;
	double F = d*e;
	double G = 1.0 + e*e;

	Eigen::Vector3d n = Eigen::Vector3d(-d, -e, 1.0).normalized();

	double L = 2.0 * a * n[2];
	double M = b * n[2];
	double N = 2 * c * n[2];


	// ----------------- Eigen stuff
	Eigen::Matrix2d m;
	m << L*G - M*F, M*E - L*F, M*E - L*F, N*E - M*F;
	m = m / (E*G - F*F);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(m);

	Eigen::Vector2d c_val = eig.eigenvalues();
	Eigen::Matrix2d c_vec = eig.eigenvectors();

	// std::cerr << "c_val:" << c_val << std::endl;
	// std::cerr << "c_vec:" << c_vec << std::endl;

	// std::cerr << "c_vec:" << c_vec(0) << " "  << c_vec(1) << std::endl;

	c_val = -c_val;

	Eigen::Vector3d v1, v2;
	v1[0] = c_vec(0);
	v1[1] = c_vec(1);
	v1[2] = 0; //d * v1[0] + e * v1[1];

	v2[0] = c_vec(2);
	v2[1] = c_vec(3);
	v2[2] = 0; //d * v2[0] + e * v2[1];


	// v1 = v1.normalized();
	// v2 = v2.normalized();

	Eigen::Vector3d v1global = ref[0] * v1[0] + ref[1] * v1[1] + ref[2] * v1[2];
	Eigen::Vector3d v2global = ref[0] * v2[0] + ref[1] * v2[1] + ref[2] * v2[2];

	v1global.normalize();
	v2global.normalize();

	v1global *= c_val(0);
	v2global *= c_val(1);

	if (c_val[0] > c_val[1])
	{
		curv[i] = std::vector<double>(2);
		curv[i][0] = c_val(1);
		curv[i][1] = c_val(0);
		curvDir[i] = std::vector<Eigen::Vector3d>(2);
		curvDir[i][0] = v2global;
		curvDir[i][1] = v1global;
	}
	else
	{
		curv[i] = std::vector<double>(2);
		curv[i][0] = c_val(0);
		curv[i][1] = c_val(1);
		curvDir[i] = std::vector<Eigen::Vector3d>(2);
		curvDir[i][0] = v1global;
		curvDir[i][1] = v2global;
	}
	// ---- end Eigen stuff
}

IGL_INLINE void CurvatureCalculator::getKRing(const int start, const double r, std::vector<int>&vv)
{
	int bufsize = vertices.rows();
	vv.reserve(bufsize);
	std::list<std::pair<int, int> > queue;
	bool* visited = (bool*)calloc(bufsize, sizeof(bool));
	queue.push_back(std::pair<int, int>(start, 0));
	visited[start] = true;
	while (!queue.empty())
	{
		int toVisit = queue.front().first;
		int distance = queue.front().second;
		queue.pop_front();
		vv.push_back(toVisit);
		if (distance < (int)r)
		{
			for (unsigned int i = 0; i < vertex_to_vertices[toVisit].size(); i++)
			{
				int neighbor = vertex_to_vertices[toVisit][i];
				if (!visited[neighbor])
				{
					queue.push_back(std::pair<int, int>(neighbor, distance + 1));
					visited[neighbor] = true;
				}
			}
		}
	}
	free(visited);
	return;
}

IGL_INLINE Eigen::Vector3d CurvatureCalculator::project(Eigen::Vector3d v, Eigen::Vector3d  vp, Eigen::Vector3d ppn)
{
	return (vp - (ppn * ((vp - v).dot(ppn))));
}

IGL_INLINE void CurvatureCalculator::computeReferenceFrame(int i, Eigen::Vector3d normal, std::vector<Eigen::Vector3d>& ref)
{

	Eigen::Vector3d longest_v = Eigen::Vector3d::Zero();
	longest_v = Eigen::Vector3d(vertices.row(vertex_to_vertices[i][0]));

	longest_v = (project(vertices.row(i), longest_v, normal) - Eigen::Vector3d(vertices.row(i))).normalized();

	/* L'ultimo asse si ottiene come prodotto vettoriale tra i due
	 * calcolati */
	Eigen::Vector3d y_axis = (normal.cross(longest_v)).normalized();
	ref[0] = longest_v;
	ref[1] = y_axis;
	ref[2] = normal;
}

IGL_INLINE void CurvatureCalculator::getAverageNormal(int j, std::vector<int> vv, Eigen::Vector3d& normal)
{
	normal = (vertex_normals.row(j)).normalized();
	if (localMode)
		return;

	for (unsigned int i = 0; i < vv.size(); i++)
	{
		normal += vertex_normals.row(vv[i]).normalized();
	}
	normal.normalize();
}

IGL_INLINE void CurvatureCalculator::getProjPlane(int j, std::vector<int> vv, Eigen::Vector3d& ppn)
{
	int nr;
	float a, b, c;
	float nx, ny, nz;
	float abcq;

	a = b = c = 0;

	if (localMode)
	{
		for (unsigned int i = 0; i < vertex_to_faces.at(j).size(); ++i)
		{
			Eigen::Vector3d faceNormal = face_normals.row(vertex_to_faces.at(j).at(i));
			a += faceNormal[0];
			b += faceNormal[1];
			c += faceNormal[2];
		}
	}
	else
	{
		for (unsigned int i = 0; i < vv.size(); ++i)
		{
			a += vertex_normals.row(vv[i])[0];
			b += vertex_normals.row(vv[i])[1];
			c += vertex_normals.row(vv[i])[2];
		}
	}
	nr = rotateForward(&a, &b, &c);
	abcq = a*a + b*b + c*c;
	nx = sqrt(a*a / abcq);
	ny = sqrt(b*b / abcq);
	nz = sqrt(1 - nx*nx - ny*ny);
	rotateBackward(nr, &a, &b, &c);
	rotateBackward(nr, &nx, &ny, &nz);

	ppn = chooseMax(Eigen::Vector3d(nx, ny, nz), Eigen::Vector3d(a, b, c), a * b);
	ppn.normalize();
}


IGL_INLINE double CurvatureCalculator::getAverageEdge()
{
	double sum = 0;
	int count = 0;

	for (int i = 0; i < faces.rows(); i++)
	{
		for (short unsigned j = 0; j < 3; j++)
		{
			Eigen::Vector3d p1 = vertices.row(faces.row(i)[j]);
			Eigen::Vector3d p2 = vertices.row(faces.row(i)[(j + 1) % 3]);

			double l = (p1 - p2).norm();

			sum += l;
			++count;
		}
	}

	return (sum / (double)count);
}


IGL_INLINE void CurvatureCalculator::applyProjOnPlane(Eigen::Vector3d ppn, std::vector<int> vin, std::vector<int> &vout)
{
	for (std::vector<int>::iterator vpi = vin.begin(); vpi != vin.end(); ++vpi)
		if (vertex_normals.row(*vpi) * ppn > 0.0f)
			vout.push_back(*vpi);
}

IGL_INLINE void CurvatureCalculator::applyMontecarlo(std::vector<int>& vin, std::vector<int> *vout)
{
	if (montecarloN >= vin.size())
	{
		*vout = vin;
		return;
	}

	float p = ((float)montecarloN) / (float)vin.size();
	for (std::vector<int>::iterator vpi = vin.begin(); vpi != vin.end(); ++vpi)
	{
		float r;
		if ((r = ((float)rand() / RAND_MAX)) < p)
		{
			vout->push_back(*vpi);
		}
	}
}

IGL_INLINE void CurvatureCalculator::computeCurvature()
{
	using namespace std;

	//CHECK che esista la mesh
	size_t vertices_count = vertices.rows();

	if (vertices_count <= 0)
		return;

	curvDir = std::vector< std::vector<Eigen::Vector3d> >(vertices_count);
	curv = std::vector<std::vector<double> >(vertices_count);



	scaledRadius = getAverageEdge()*sphereRadius;

	std::vector<int> vv;
	std::vector<int> vvtmp;
	Eigen::Vector3d normal;

	//double time_spent;
	//double searchtime=0, ref_time=0, fit_time=0, final_time=0;

	for (size_t i = 0; i < vertices_count; ++i)
	{
		vv.clear();
		vvtmp.clear();
		Eigen::Vector3d me = vertices.row(i);
		switch (st)
		{
		case SPHERE_SEARCH:
			//getSphere(i,scaledRadius,vv,6);
			break;
		case K_RING_SEARCH:
			getKRing(i, kRing, vv);
			break;
		default:
			fprintf(stderr, "Error: search type not recognized");
			return;
		}

		std::vector<Eigen::Vector3d> ref(3);
		if (vv.size() < 6)
		{
			std::cerr << "Could not compute curvature of radius " << scaledRadius << endl;
			return;
		}


		if (projectionPlaneCheck)
		{
			vvtmp.reserve(vv.size());
			applyProjOnPlane(vertex_normals.row(i), vv, vvtmp);
			if (vvtmp.size() >= 6 && vvtmp.size() < vv.size())
				vv = vvtmp;

		}


		switch (nt)
		{
		case AVERAGE:
			getAverageNormal(i, vv, normal);
			break;
		case PROJ_PLANE:
			getProjPlane(i, vv, normal);
			break;
		default:
			fprintf(stderr, "Error: normal type not recognized");
			return;
		}
		if (vv.size() < 6)
		{
			std::cerr << "Could not compute curvature of radius " << scaledRadius << endl;
			return;
		}
		if (montecarlo)
		{
			if (montecarloN < 6)
				break;
			vvtmp.reserve(vv.size());
			applyMontecarlo(vv, &vvtmp);
			vv = vvtmp;
		}

		if (vv.size() < 6)
			return;
		computeReferenceFrame(i, normal, ref);

		Quadric q;
		fitQuadric(me, ref, vv, &q);
		finalEigenStuff(i, ref, q);
	}

	lastRadius = sphereRadius;
	curvatureComputed = true;
}

IGL_INLINE void CurvatureCalculator::printCurvature(std::string outpath)
{
	using namespace std;
	if (!curvatureComputed)
		return;

	std::ofstream of;
	of.open(outpath.c_str());

	if (!of)
	{
		fprintf(stderr, "Error: could not open output file %s\n", outpath.c_str());
		return;
	}

	int vertices_count = vertices.rows();
	of << vertices_count << endl;
	for (int i = 0; i < vertices_count; i++)
	{
		of << curv[i][0] << " " << curv[i][1] << " " << curvDir[i][0][0] << " " << curvDir[i][0][1] << " " << curvDir[i][0][2] << " " <<
			curvDir[i][1][0] << " " << curvDir[i][1][1] << " " << curvDir[i][1][2] << endl;
	}

	of.close();

}

template <
	typename DerivedV,
	typename DerivedF,
	typename DerivedPD1,
	typename DerivedPD2,
	typename DerivedPV1,
	typename DerivedPV2>
	IGL_INLINE void igl::principal_curvature(
	const Eigen::PlainObjectBase<DerivedV>& V,
	const Eigen::PlainObjectBase<DerivedF>& F,
	Eigen::PlainObjectBase<DerivedPD1>& PD1,
	Eigen::PlainObjectBase<DerivedPD2>& PD2,
	Eigen::PlainObjectBase<DerivedPV1>& PV1,
	Eigen::PlainObjectBase<DerivedPV2>& PV2,
	unsigned radius,
	bool useKring)
{
	using namespace std;

	// Preallocate memory
	PD1.resize(V.rows(), 3);
	PD2.resize(V.rows(), 3);

	// Preallocate memory
	PV1.resize(V.rows(), 1);
	PV2.resize(V.rows(), 1);

	// Precomputation
	CurvatureCalculator cc;
	cc.init(V.template cast<double>(), F.template cast<int>());
	cc.sphereRadius = radius;

	if (useKring)
	{
		cc.kRing = radius;
		cc.st = K_RING_SEARCH;
	}

	// Compute
	cc.computeCurvature();

	// Copy it back
	for (unsigned i = 0; i < V.rows(); i++)
	{
		Eigen::Vector3d d1;
		Eigen::Vector3d d2;
		PD1.row(i) << cc.curvDir[i][0][0], cc.curvDir[i][0][1], cc.curvDir[i][0][2];
		PD2.row(i) << cc.curvDir[i][1][0], cc.curvDir[i][1][1], cc.curvDir[i][1][2];
		PD1.row(i).normalize();
		PD2.row(i).normalize();

		if (std::isnan(PD1(i, 0)) || std::isnan(PD1(i, 1)) || std::isnan(PD1(i, 2)) || std::isnan(PD2(i, 0)) || std::isnan(PD2(i, 1)) || std::isnan(PD2(i, 2)))
		{
			PD1.row(i) << 0, 0, 0;
			PD2.row(i) << 0, 0, 0;
		}

		PV1(i) = cc.curv[i][0];
		PV2(i) = cc.curv[i][1];

		if (PD1.row(i) * PD2.row(i).transpose() > 10e-6)
		{
			cerr << "PRINCIPAL_CURVATURE: Something is wrong with vertex: i" << endl;
			PD1.row(i) *= 0;
			PD2.row(i) *= 0;
		}
	}

}


/*
#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
// generated by autoexplicit.sh
template void igl::principal_curvature<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, unsigned int, bool);
#endif
*/

#endif
