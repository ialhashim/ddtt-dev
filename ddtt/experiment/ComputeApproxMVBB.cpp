// ========================================================================================
//  ApproxMVBB
//  Copyright (C) 2014 by Gabriel Nützi <nuetzig (at) imes (d0t) mavt (d0t) ethz (døt) ch>
//
//  This Source Code Form is subject to the terms of the Mozilla Public
//  License, v. 2.0. If a copy of the MPL was not distributed with this
//  file, You can obtain one at http://mozilla.org/MPL/2.0/.
// ========================================================================================

#include "ComputeApproxMVBB.h"
#include <random>

#include <limits>

namespace Diameter{

	typedef struct {
		double *extremity1;
		double *extremity2;
		double squareDiameter;
		int reduction_mode;
	} typeSegment;

	/* util */


	typedef struct {
		int n;
		int nalloc;
		typeSegment *seg;
	} typeListOfSegments;


	void *_AllocateListOfPoints(const int n, const int dim)
	{
		void *b = NULL;
		double **list, **dd;
		double *d;
		int i;

		if (n <= 0 || dim <= 0) return(NULL);

		b = (void*)malloc(n         * sizeof(double*)
			+ n * dim * sizeof(double));
		if (b == (void *)NULL)  return((void *)NULL);

		dd = list = (double **)b;
		dd += n;

		d = (double*)dd;
		for (i = 0; i<n; i++, d += dim) list[i] = d;

		return(b);
	}








	void *_AllocateListOfSegments(const int n)
	{
		void *b = NULL;
		typeSegment *d;
		int i;


		if (n <= 0) return(NULL);

		b = (void*)malloc(n * sizeof(typeSegment));

		if (b == (void *)NULL)  return((void *)NULL);

		d = (typeSegment *)b;
		for (i = 0; i<n; i++) {
			d[i].extremity1 = (double*)NULL;
			d[i].extremity2 = (double*)NULL;
			d[i].squareDiameter = 0.0;
			d[i].reduction_mode = 0;
		}

		return(b);
	}







#define _NALLOC_ 20


	int _AddSegmentToList(typeSegment *s, typeListOfSegments *list)
	{
		typeSegment *d;

		if (list->nalloc <= 0)
			list->n = list->nalloc = 0;

		if (list->n >= list->nalloc) {

			d = (typeSegment *)_AllocateListOfSegments(list->nalloc + _NALLOC_);
			if (d == NULL) {
				return(0);
			}

			if (list->nalloc > 0) {
				memcpy(d, list->seg, list->nalloc*sizeof(typeSegment));
				free(list->seg);
			}

			list->nalloc += _NALLOC_;
			list->seg = d;

		}

		list->seg[list->n] = *s;
		list->n++;

		return(1);
	}


	extern void _VerboseWhenReducing();
	extern void _NoVerboseWhenReducing();
	extern  int _GetVerboseWhenReducing();

	extern void _SetTryToReducePInIterative(int t);
	extern void _DoTryToReducePInIterative();
	extern void _DoNotTryToReducePInIterative();

	extern void _SetReductionModeInIterative(int m);
	extern  int _GetReductionModeInIterative();
	extern void _SetReductionModeOfDiameter(int m);
	extern  int _GetReductionModeOfDiameter();
	extern void _SetReductionModeOfDbleNorm(int m);
	extern  int _GetReductionModeOfDbleNorm();


	extern void _DoTryToReduceQ();
	extern void _DoNotTryToReduceQ();
	extern  int _GetTryToReduceQ();


	extern void _SetQscanToForward();
	extern void _SetQscanToBackward();
	extern  int _GetQscan();

	extern void _DoTryToGetTightBounds();
	extern void _DoNotTryToGetTightBounds();
	extern  int _GetTightBounds();

	extern int _LastPointOutsideSphereWithDiameter(typeSegment *theSeg,
		const double squareDiameter,
		double **theList,
		const int first,
		int *last,
		const int dim,
		const int _reduction_mode_);

	extern int _LastPointOutsideSphereAndBoundWithDiameter(typeSegment *theSeg,
		const double squareDiameter,
		double **theList,
		const int first,
		int *last,
		const int dim,
		const int _reduction_mode_,
		double *bound);

	extern int _FarthestPointFromSphere(typeSegment *theSeg,
		double **theList,
		const int first,
		int *last,
		const int dim,
		const int _reduction_mode_);




	extern void _CountPointsInSpheres(typeSegment *theSeg,
		const double squareDiameter,
		double **theList,
		const int first,
		const int last,
		const int dim);




	extern double _MaximalSegmentInTwoLists(typeSegment *theSeg,
		const int index1,
		double **theList1,
		int *first1,
		int *last1,
		double **theList2,
		int *first2,
		int *last2,
		int dim);

	extern double _MaximalSegmentInOneList(typeSegment *theSeg,
		const int index,
		double **theList,
		int *first,
		int *last,
		const int dim);

	extern double _MaximalDistanceFromPoint(int *index,
		const double *ref,
		double **theList,
		const int first,
		const int last,
		const int dim);







	extern double _QuadraticDiameterInOneList(typeSegment *theDiam,
		double **theList,
		const int first,
		const int last,
		const int dim);

	extern double _QuadraticDiameterInTwoLists(typeSegment *theDiam,
		int   *index1,
		int   *index2,
		double **theList1,
		const int first1,
		const int last1,
		double **theList2,
		const int first2,
		const int last2,
		const int dim);








	extern void _SwapPoints(double **theList, const int i, const int j);





	typedef struct {
		int c1;
		int c2;
	} typeCounter;

	extern void _InitCounter(typeCounter *c);
	extern void _AddToCounter(typeCounter *c, const int i);
	extern double _GetCounterAverage(typeCounter *c, const int i);

#ifdef _STATS_
	extern void _InitScalarProductCounter();
	extern double _GetScalarProductAverage(int n);
#endif




	/* square distances
	*/
	extern double _SquareDistance(const double *a, const double *b, const int dim);
	extern double _SquareDistance3D(const double *a, const double *b);
	extern double _SquareDistance2D(const double *a, const double *b);

	/* dot products
	ab.cd = ( (b[i]-a[i]) . (d[i]-c[i]) )
	*/
	extern double _ScalarProduct(const double *a, const double *b,
		const double *c, const double *d, const int dim);
	extern double _ScalarProduct3D(const double *a, const double *b,
		const double *c, const double *d);
	extern double _ScalarProduct2D(const double *a, const double *b,
		const double *c, const double *d);


	extern int _FindPointInList(double **theList,
		const int first,
		const int last,
		double x0,
		double x1);


	/* verbose when reducing
	*/
	static int _verbose_when_reducing_ = 0;

	void _VerboseWhenReducing()
	{
		_verbose_when_reducing_ = 1;
	}
	void _NoVerboseWhenReducing()
	{
		_verbose_when_reducing_ = 0;
	}
	int _GetVerboseWhenReducing()
	{
		return(_verbose_when_reducing_);
	}



	/* reduction in the iterative search of the double normal
	*/
	static int _reduction_mode_in_iterative_ = 1;
	void _SetReductionModeInIterative(int m)
	{
		switch (m) {
		case 0:
		case 1:
			_reduction_mode_in_iterative_ = m;
			break;
		case 2:
		default:
			break;
		}
	}
	int _GetReductionModeInIterative()
	{
		return(_reduction_mode_in_iterative_);
	}




	/* 'reduction' of diameter
	*/
	static int _reduction_mode_of_diameter_ = 1;
	void _SetReductionModeOfDiameter(int m)
	{
		switch (m) {
		case 0:
		case 1:
		case 2:
			_reduction_mode_of_diameter_ = m;
			break;
		default:
			break;
		}
	}
	int _GetReductionModeOfDiameter()
	{
		return(_reduction_mode_of_diameter_);
	}



	/* 'reduction' of double normals
	*/
	static int _reduction_mode_of_dbleNorm_ = 1;
	void _SetReductionModeOfDbleNorm(int m)
	{
		switch (m) {
		case 0:
		case 1:
		case 2:
			_reduction_mode_of_dbleNorm_ = m;
			break;
		default:
			break;
		}
	}
	int _GetReductionModeOfDbleNorm()
	{
		return(_reduction_mode_of_dbleNorm_);
	}




	/* reduction by processing couple of double normals
	*/
	static int _try_to_reduce_Q_ = 1;

	void _DoTryToReduceQ()
	{
		_try_to_reduce_Q_ = 1;
	}
	void _DoNotTryToReduceQ()
	{
		_try_to_reduce_Q_ = 0;
	}
	int _GetTryToReduceQ()
	{
		return(_try_to_reduce_Q_);
	}





	static int _Q_scan_ = 0;
	void _SetQscanToForward()
	{
		_Q_scan_ = 1;
	}
	void _SetQscanToBackward()
	{
		_Q_scan_ = 0;
	}
	int _GetQscan()
	{
		return(_Q_scan_);
	}




	static int _tight_bounds_ = 0;
	void _DoTryToGetTightBounds()
	{
		_tight_bounds_ = 1;
	}
	void _DoNotTryToGetTightBounds()
	{
		_tight_bounds_ = 0;
	}
	int _GetTightBounds()
	{
		return(_tight_bounds_);
	}






	/* partially sort a list of points

	given a "diameter", points which are 'outside'
	the sphere with a given threshold are put at the beginning
	of the list,
	the returned value is the index of the last point 'outside'
	ie
	points from #first    to #index are outside
	"     "    #index+1 to #last  are inside

	If there are no points outside, #first-1 is returned.

	Given a segment [AB], its centre C is (A+B)/2.
	The dot product MA.MB is equal to MC^2 - |AB|^2 / 4

	1. Being naive, ie to test if points are outside
	the ball of diameter [AB], yield to condition
	MA.MB > 0

	2. Being a little smarter, ie to  test if points are outside
	the ball of center (A+B)/2 and diameter D (assume D>|AB|)
	yield to condition
	MA.MB > ( D^2 - |AB|^2 ) / 4

	*/

	int _LastPointOutsideSphereWithDiameter(typeSegment *theSeg,
		const double squareDiameter,
		double **theList,
		const int first,
		int *last,
		const int dim,
		const int _reduction_mode_)
	{
		int i;
		int index = first - 1;
		int l = *last;

		double mamb, am2;

		double minThreshold;
		double medThreshold;
		double maxThreshold;

		double ab = sqrt(theSeg->squareDiameter);
		double ab2 = theSeg->squareDiameter;

		double R = sqrt(squareDiameter);
		double R2 = squareDiameter;

		if (first > *last) return(first - 1);

		if (squareDiameter <= theSeg->squareDiameter) {

			maxThreshold = medThreshold = 0.0;
			minThreshold = -0.23205080756887729352 * theSeg->squareDiameter;

		}
		else {

			minThreshold = (R - .86602540378443864676*ab)
				* (R - .86602540378443864676*ab)
				- 0.25 * ab2;

			medThreshold = 0.5 * ab2 + 0.25 * R2
				- .43301270189221932338 * R * sqrt(4 * ab2 - R2);

			maxThreshold = 0.25 * (R2 - ab2);
		}



		/* 0 : no reduction
		1 : just inside the smallest sphere
		2 : complete reduction
		*/
		switch (_reduction_mode_) {

		case 0:

			/* NO REDUCTION CASE
			*/
			if (dim == 2) {
				for (i = first; i <= l; i++) {
					mamb = _ScalarProduct2D(theList[i], theSeg->extremity1,
						theList[i], theSeg->extremity2);
					if (mamb > maxThreshold) {
						index++;
						_SwapPoints(theList, index, i);
					}
				}
				return(index);
			}

			if (dim == 3) {
				for (i = first; i <= l; i++) {
					mamb = _ScalarProduct3D(theList[i], theSeg->extremity1,
						theList[i], theSeg->extremity2);
					if (mamb > maxThreshold) {
						index++;
						_SwapPoints(theList, index, i);
					}
				}
				return(index);
			}

			for (i = first; i <= l; i++) {
				mamb = _ScalarProduct(theList[i], theSeg->extremity1,
					theList[i], theSeg->extremity2, dim);
				if (mamb > maxThreshold) {
					index++;
					_SwapPoints(theList, index, i);
				}
			}
			return(index);

			/* END
			NO REDUCTION CASE
			*/
			break;

		default:
		case 1:

			/* REDUCTION INSIDE THE SMALLEST SPHERE
			*/

			if (dim == 2) {
				for (i = first; i <= l; i++) {
					mamb = _ScalarProduct2D(theList[i], theSeg->extremity1,
						theList[i], theSeg->extremity2);
					if (mamb > maxThreshold) {
						index++;
						_SwapPoints(theList, index, i);
					}
					else if (mamb <= minThreshold) {
						_SwapPoints(theList, i, l);
						i--;   l--;
					}
				}
				*last = l;
				return(index);
			}

			if (dim == 3) {
				for (i = first; i <= l; i++) {
					mamb = _ScalarProduct3D(theList[i], theSeg->extremity1,
						theList[i], theSeg->extremity2);
					if (mamb > maxThreshold) {
						index++;
						_SwapPoints(theList, index, i);
					}
					else if (mamb <= minThreshold) {
						_SwapPoints(theList, i, l);
						i--;   l--;
					}
				}
				*last = l;
				return(index);
			}

			for (i = first; i <= l; i++) {
				mamb = _ScalarProduct(theList[i], theSeg->extremity1,
					theList[i], theSeg->extremity2, dim);
				if (mamb > maxThreshold) {
					index++;
					_SwapPoints(theList, index, i);
				}
				else if (mamb <= minThreshold) {
					_SwapPoints(theList, i, l);
					i--;   l--;
				}
			}
			*last = l;
			return(index);

			/* END
			REDUCTION INSIDE THE SMALLEST SPHERE
			*/
			break;

		case 2:

			/* COMPLETE REDUCTION

			On suppose implicitement
			que l'ensemble des points est
			compris dans l'intersection des spheres centrees
			sur A et B et de rayon |AB|.
			En effet, les conditions (pour un rayon dans
			l'intervalle [min,med]) pour l'elimination ne le
			verifient pas.
			*/
			if (dim == 2) {

				if (squareDiameter <= theSeg->squareDiameter) {

					for (i = first; i <= l; i++) {
						mamb = _ScalarProduct2D(theList[i], theSeg->extremity1,
							theList[i], theSeg->extremity2);
						if (mamb > maxThreshold) {
							index++;
							_SwapPoints(theList, index, i);
						}
						else if (mamb > minThreshold) {
							am2 = _SquareDistance2D(theList[i], theSeg->extremity1);
							if (3.0 * (am2 * ab2 - (am2 - mamb)*(am2 - mamb)) - mamb*mamb < 0) {
								_SwapPoints(theList, i, l);
								i--;   l--;
							}
						}
						else {
							/* if ( mamb <= minThreshold ) */
							_SwapPoints(theList, i, l);
							i--;   l--;
						}
					}

				}
				else {

					for (i = first; i <= l; i++) {
						mamb = _ScalarProduct2D(theList[i], theSeg->extremity1,
							theList[i], theSeg->extremity2);
						if (mamb > maxThreshold) {
							index++;
							_SwapPoints(theList, index, i);
						}
						else if (mamb > medThreshold) {
							continue;
						}
						else if (mamb > minThreshold) {
							am2 = _SquareDistance2D(theList[i], theSeg->extremity1);
							if (3.0 * (am2 * ab2 - (am2 - mamb)*(am2 - mamb))
								- (R2 - ab2 - mamb)*(R2 - ab2 - mamb) < 0) {
								_SwapPoints(theList, i, l);
								i--;   l--;
							}
						}
						else {
							/* if ( mamb <= minThreshold ) */
							_SwapPoints(theList, i, l);
							i--;   l--;
						}
					}

				}
				*last = l;
				return(index);
			}

			if (dim == 3) {

				if (squareDiameter <= theSeg->squareDiameter) {

					for (i = first; i <= l; i++) {
						mamb = _ScalarProduct3D(theList[i], theSeg->extremity1,
							theList[i], theSeg->extremity2);
						if (mamb > maxThreshold) {
							index++;
							_SwapPoints(theList, index, i);
						}
						else if (mamb > minThreshold) {
							am2 = _SquareDistance3D(theList[i], theSeg->extremity1);
							if (3.0 * (am2 * ab2 - (am2 - mamb)*(am2 - mamb)) - mamb*mamb < 0) {
								_SwapPoints(theList, i, l);
								i--;   l--;
							}
						}
						else {
							/* if ( mamb <= minThreshold ) */
							_SwapPoints(theList, i, l);
							i--;   l--;
						}
					}

				}
				else {

					for (i = first; i <= l; i++) {
						mamb = _ScalarProduct3D(theList[i], theSeg->extremity1,
							theList[i], theSeg->extremity2);
						if (mamb > maxThreshold) {
							index++;
							_SwapPoints(theList, index, i);
						}
						else if (mamb > medThreshold) {
							continue;
						}
						else if (mamb > minThreshold) {
							am2 = _SquareDistance3D(theList[i], theSeg->extremity1);
							if (3.0 * (am2 * ab2 - (am2 - mamb)*(am2 - mamb))
								- (R2 - ab2 - mamb)*(R2 - ab2 - mamb) < 0) {
								_SwapPoints(theList, i, l);
								i--;   l--;
							}
						}
						else {
							/* if ( mamb <= minThreshold ) */
							_SwapPoints(theList, i, l);
							i--;   l--;
						}
					}

				}
				*last = l;
				return(index);
			}

			if (squareDiameter <= theSeg->squareDiameter) {

				for (i = first; i <= l; i++) {
					mamb = _ScalarProduct(theList[i], theSeg->extremity1,
						theList[i], theSeg->extremity2, dim);
					if (mamb > maxThreshold) {
						index++;
						_SwapPoints(theList, index, i);
					}
					else if (mamb > minThreshold) {
						am2 = _SquareDistance(theList[i], theSeg->extremity1, dim);
						if (3.0 * (am2 * ab2 - (am2 - mamb)*(am2 - mamb)) - mamb*mamb < 0) {
							_SwapPoints(theList, i, l);
							i--;   l--;
						}
					}
					else {
						/* if ( mamb <= minThreshold ) */
						_SwapPoints(theList, i, l);
						i--;   l--;
					}
				}

			}
			else {

				for (i = first; i <= l; i++) {
					mamb = _ScalarProduct(theList[i], theSeg->extremity1,
						theList[i], theSeg->extremity2, dim);
					if (mamb > maxThreshold) {
						index++;
						_SwapPoints(theList, index, i);
					}
					else if (mamb > medThreshold) {
						continue;
					}
					else if (mamb > minThreshold) {
						am2 = _SquareDistance(theList[i], theSeg->extremity1, dim);
						if (3.0 * (am2 * ab2 - (am2 - mamb)*(am2 - mamb))
							- (R2 - ab2 - mamb)*(R2 - ab2 - mamb) < 0) {
							_SwapPoints(theList, i, l);
							i--;   l--;
						}
					}
					else {
						/* if ( mamb <= minThreshold ) */
						_SwapPoints(theList, i, l);
						i--;   l--;
					}
				}

			}
			*last = l;
			return(index);

			/* END
			COMPLETE REDUCTION
			*/
		}
		return(-1);
	}


















	int _LastPointOutsideSphereAndBoundWithDiameter(typeSegment *theSeg,
		const double squareDiameter,
		double **theList,
		const int first,
		int *last,
		const int dim,
		const int _reduction_mode_,
		double *bound)
	{
		int i;
		int index = first - 1;
		int l = *last;

		double mamb, am2;

		double minThreshold;
		double medThreshold;
		double maxThreshold;

		double ab = sqrt(theSeg->squareDiameter);
		double ab2 = theSeg->squareDiameter;

		double R = sqrt(squareDiameter);
		double R2 = squareDiameter;

		/* bound
		c'est la plus grande valeur pour les points
		l'interieur de la sphere
		si squareDiameter > theSeg->squareDiameter
		cette valeur peut etre positive
		*/
		double b = *bound = (-theSeg->squareDiameter * 0.25);


		if (squareDiameter <= theSeg->squareDiameter) {

			maxThreshold = medThreshold = 0.0;
			minThreshold = -0.23205080756887729352 * theSeg->squareDiameter;

			return(_LastPointOutsideSphereWithDiameter(theSeg, squareDiameter,
				theList, first, last, dim,
				_reduction_mode_));

		}
		else {

			minThreshold = (R - .86602540378443864676*ab)
				* (R - .86602540378443864676*ab)
				- 0.25 * ab2;

			medThreshold = 0.5 * ab2 + 0.25 * R2
				- .43301270189221932338 * R * sqrt(4 * ab2 - R2);

			maxThreshold = 0.25 * (R2 - ab2);
		}



		/* 0 : no reduction
		1 : just inside the smallest sphere
		2 : complete reduction
		*/
		switch (_reduction_mode_) {

		case 0:

			/* NO REDUCTION CASE
			*/
			if (dim == 2) {
				for (i = first; i <= l; i++) {
					mamb = _ScalarProduct2D(theList[i], theSeg->extremity1,
						theList[i], theSeg->extremity2);
					if (mamb > maxThreshold) {
						index++;
						_SwapPoints(theList, index, i);
						continue;
					}
					if (b < mamb) b = mamb;
				}
				*bound = b;
				return(index);
			}

			if (dim == 3) {
				for (i = first; i <= l; i++) {
					mamb = _ScalarProduct3D(theList[i], theSeg->extremity1,
						theList[i], theSeg->extremity2);
					if (mamb > maxThreshold) {
						index++;
						_SwapPoints(theList, index, i);
						continue;
					}
					if (b < mamb) b = mamb;
				}
				*bound = b;
				return(index);
			}

			for (i = first; i <= l; i++) {
				mamb = _ScalarProduct(theList[i], theSeg->extremity1,
					theList[i], theSeg->extremity2, dim);
				if (mamb > maxThreshold) {
					index++;
					_SwapPoints(theList, index, i);
					continue;
				}
				if (b < mamb) b = mamb;
			}
			*bound = b;
			return(index);

			/* END
			NO REDUCTION CASE
			*/
			break;

		default:
		case 1:

			/* REDUCTION INSIDE THE SMALLEST SPHERE
			*/

			if (dim == 2) {
				for (i = first; i <= l; i++) {
					mamb = _ScalarProduct2D(theList[i], theSeg->extremity1,
						theList[i], theSeg->extremity2);
					if (mamb > maxThreshold) {
						index++;
						_SwapPoints(theList, index, i);
						continue;
					}
					if (mamb <= minThreshold) {
						_SwapPoints(theList, i, l);
						i--;   l--;
						continue;
					}
					if (b < mamb) b = mamb;
				}
				*last = l;
				*bound = b;
				return(index);
			}

			if (dim == 3) {
				for (i = first; i <= l; i++) {
					mamb = _ScalarProduct3D(theList[i], theSeg->extremity1,
						theList[i], theSeg->extremity2);
					if (mamb > maxThreshold) {
						index++;
						_SwapPoints(theList, index, i);
						continue;
					}
					if (mamb <= minThreshold) {
						_SwapPoints(theList, i, l);
						i--;   l--;
						continue;
					}
					if (b < mamb) b = mamb;
				}
				*last = l;
				*bound = b;
				return(index);
			}

			for (i = first; i <= l; i++) {
				mamb = _ScalarProduct(theList[i], theSeg->extremity1,
					theList[i], theSeg->extremity2, dim);
				if (mamb > maxThreshold) {
					index++;
					_SwapPoints(theList, index, i);
					continue;
				}
				if (mamb <= minThreshold) {
					_SwapPoints(theList, i, l);
					i--;   l--;
					continue;
				}
				if (b < mamb) b = mamb;
			}
			*last = l;
			*bound = b;
			return(index);

			/* END
			REDUCTION INSIDE THE SMALLEST SPHERE
			*/
			break;

		case 2:

			/* COMPLETE REDUCTION

			On suppose implicitement
			que l'ensemble des points est
			compris dans l'intersection des spheres centrees
			sur A et B et de rayon |AB|.
			En effet, les conditions (pour un rayon dans
			l'intervalle [min,med]) pour l'elimination ne le
			verifient pas.
			*/
			if (dim == 2) {

				if (squareDiameter <= theSeg->squareDiameter) {

					for (i = first; i <= l; i++) {
						mamb = _ScalarProduct2D(theList[i], theSeg->extremity1,
							theList[i], theSeg->extremity2);
						if (mamb > maxThreshold) {
							index++;
							_SwapPoints(theList, index, i);
							continue;
						}
						if (mamb > minThreshold) {
							am2 = _SquareDistance2D(theList[i], theSeg->extremity1);
							if (3.0 * (am2 * ab2 - (am2 - mamb)*(am2 - mamb)) - mamb*mamb < 0) {
								_SwapPoints(theList, i, l);
								i--;   l--;
								continue;
							}
							if (b < mamb) b = mamb;
							continue;
						}
						/* if ( mamb <= minThreshold ) */
						_SwapPoints(theList, i, l);
						i--;   l--;
					}

				}
				else {

					for (i = first; i <= l; i++) {
						mamb = _ScalarProduct2D(theList[i], theSeg->extremity1,
							theList[i], theSeg->extremity2);
						if (mamb > maxThreshold) {
							index++;
							_SwapPoints(theList, index, i);
							continue;
						}
						if (mamb > medThreshold) {
							if (b < mamb) b = mamb;
							continue;
						}
						if (mamb > minThreshold) {
							am2 = _SquareDistance2D(theList[i], theSeg->extremity1);
							if (3.0 * (am2 * ab2 - (am2 - mamb)*(am2 - mamb))
								- (R2 - ab2 - mamb)*(R2 - ab2 - mamb) < 0) {
								_SwapPoints(theList, i, l);
								i--;   l--;
								continue;
							}
							if (b < mamb) b = mamb;
							continue;
						}
						/* if ( mamb <= minThreshold ) */
						_SwapPoints(theList, i, l);
						i--;   l--;
					}

				}
				*last = l;
				*bound = b;
				return(index);
			}

			if (dim == 3) {

				if (squareDiameter <= theSeg->squareDiameter) {

					for (i = first; i <= l; i++) {
						mamb = _ScalarProduct3D(theList[i], theSeg->extremity1,
							theList[i], theSeg->extremity2);
						if (mamb > maxThreshold) {
							index++;
							_SwapPoints(theList, index, i);
							continue;
						}
						if (mamb > minThreshold) {
							am2 = _SquareDistance3D(theList[i], theSeg->extremity1);
							if (3.0 * (am2 * ab2 - (am2 - mamb)*(am2 - mamb)) - mamb*mamb < 0) {
								_SwapPoints(theList, i, l);
								i--;   l--;
								continue;
							}
							if (b < mamb) b = mamb;
							continue;
						}
						/* if ( mamb <= minThreshold ) */
						_SwapPoints(theList, i, l);
						i--;   l--;
					}

				}
				else {

					for (i = first; i <= l; i++) {
						mamb = _ScalarProduct3D(theList[i], theSeg->extremity1,
							theList[i], theSeg->extremity2);
						if (mamb > maxThreshold) {
							index++;
							_SwapPoints(theList, index, i);
							continue;
						}
						if (mamb > medThreshold) {
							if (b < mamb) b = mamb;
							continue;
						}
						if (mamb > minThreshold) {
							am2 = _SquareDistance3D(theList[i], theSeg->extremity1);
							if (3.0 * (am2 * ab2 - (am2 - mamb)*(am2 - mamb))
								- (R2 - ab2 - mamb)*(R2 - ab2 - mamb) < 0) {
								_SwapPoints(theList, i, l);
								i--;   l--;
								continue;
							}
							if (b < mamb) b = mamb;
							continue;
						}
						/* if ( mamb <= minThreshold ) */
						_SwapPoints(theList, i, l);
						i--;   l--;
					}

				}
				*last = l;
				*bound = b;
				return(index);
			}

			if (squareDiameter <= theSeg->squareDiameter) {

				for (i = first; i <= l; i++) {
					mamb = _ScalarProduct(theList[i], theSeg->extremity1,
						theList[i], theSeg->extremity2, dim);
					if (mamb > maxThreshold) {
						index++;
						_SwapPoints(theList, index, i);
						continue;
					}
					if (mamb > minThreshold) {
						am2 = _SquareDistance(theList[i], theSeg->extremity1, dim);
						if (3.0 * (am2 * ab2 - (am2 - mamb)*(am2 - mamb)) - mamb*mamb < 0) {
							_SwapPoints(theList, i, l);
							i--;   l--;
							continue;
						}
						if (b < mamb) b = mamb;
						continue;
					}
					/* if ( mamb <= minThreshold ) */
					_SwapPoints(theList, i, l);
					i--;   l--;
				}

			}
			else {

				for (i = first; i <= l; i++) {
					mamb = _ScalarProduct(theList[i], theSeg->extremity1,
						theList[i], theSeg->extremity2, dim);
					if (mamb > maxThreshold) {
						index++;
						_SwapPoints(theList, index, i);
						continue;
					}
					if (mamb > medThreshold) {
						if (b < mamb) b = mamb;
						continue;
					}
					if (mamb > minThreshold) {
						am2 = _SquareDistance(theList[i], theSeg->extremity1, dim);
						if (3.0 * (am2 * ab2 - (am2 - mamb)*(am2 - mamb))
							- (R2 - ab2 - mamb)*(R2 - ab2 - mamb) < 0) {
							_SwapPoints(theList, i, l);
							i--;   l--;
							continue;
						}
						if (b < mamb) b = mamb;
						continue;
					}
					/* if ( mamb <= minThreshold ) */
					_SwapPoints(theList, i, l);
					i--;   l--;
				}

			}
			*last = l;
			*bound = b;
			return(index);

			/* END
			COMPLETE REDUCTION
			*/
		}
		return(-1);
	}
















	void _CountPointsInSpheres(typeSegment *theSeg,
		const double squareDiameter,
		double **theList,
		const int first,
		const int last,
		const int dim)
	{
		double minThreshold;
		double medThreshold;
		double maxThreshold;

		double ab = sqrt(theSeg->squareDiameter);
		double ab2 = theSeg->squareDiameter;

		double R = sqrt(squareDiameter);
		double R2 = squareDiameter;

		int n[5] = { 0, 0, 0, 0, 0 };
		int i;

		double mamb, am2;

		minThreshold = (R - .86602540378443864676*ab)
			* (R - .86602540378443864676*ab)
			- 0.25 * ab2;

		medThreshold = 0.5 * ab2 + 0.25 * R2
			- .43301270189221932338 * R * sqrt(4 * ab2 - R2);

		maxThreshold = 0.25 * (R2 - ab2);


		for (i = first; i <= last; i++) {

			mamb = _ScalarProduct(theList[i], theSeg->extremity1,
				theList[i], theSeg->extremity2, dim);

			if (mamb > maxThreshold) {
				n[0] ++;
				continue;
			}

			if (mamb > medThreshold) {
				n[1] ++;
				continue;
			}

			if (mamb > minThreshold) {
				am2 = _SquareDistance(theList[i], theSeg->extremity1, dim);
				if (3.0 * (am2 * ab2 - (am2 - mamb)*(am2 - mamb))
					- (R2 - ab2 - mamb)*(R2 - ab2 - mamb) < 0) {
					n[2] ++;
					continue;
				}
				n[3] ++;
				continue;
			}

			n[4] ++;
		}

		printf(" diametre courant = %g  -  double normale = %g\n",
			R, ab);
		printf(" %8d points dont\n", last - first + 1);
		printf(" - %6d : candidats extremites       R=%g\n", n[0], sqrt(maxThreshold + 0.25*ab2));
		printf(" - %6d : rien a faire               R=%g\n", n[1], sqrt(medThreshold + 0.25*ab2));
		printf(" - %6d : a tester                   \n", n[2] + n[3]);
		printf("   + %6d : a eliminer\n", n[2]);
		printf("   + %6d : a conserver\n", n[3]);
		printf(" - %6d : elimines directement       R=%g\n", n[4], sqrt(minThreshold + 0.25*ab2));
		printf("----------\n");
		printf(" %8d\n", n[0] + n[1] + n[2] + n[3] + n[4]);

	}


















	/* Find the farthest point from a sphere

	returned value :
	- its index if there are some points outside the sphere
	- else #(first index) - 1

	As the sphere diameter is an approximation of the set diameter,
	points verifying MA.MB <= (1.5 -sqrt(3)) AB^2
	can be removed

	*/


	int _FarthestPointFromSphere(typeSegment *theSeg,
		double **theList,
		const int first,
		int *last,
		const int dim,
		const int _reduction_mode_)
	{
		int i, l = (*last);
		int index = first - 1;
		double diff, maxdiff = 0.0;
		/* threshold = 1.5 - sqrt(3)
		*/
		double threshold = -0.23205080756887729352 * theSeg->squareDiameter;

		if (l < first) return(index);


		switch (_reduction_mode_) {

		case 0:

			/* NO REDUCTION CASE
			*/

			if (dim == 2) {
				for (i = first; i <= l; i++) {
					diff = _ScalarProduct2D(theList[i], theSeg->extremity1,
						theList[i], theSeg->extremity2);
					if (maxdiff < diff) {
						index = i;
						maxdiff = diff;
					}
				}
				return(index);
			}

			if (dim == 3) {
				for (i = first; i <= l; i++) {
					diff = _ScalarProduct3D(theList[i], theSeg->extremity1,
						theList[i], theSeg->extremity2);
					if (maxdiff < diff) {
						index = i;
						maxdiff = diff;
					}
				}
				return(index);
			}

			for (i = first; i <= l; i++) {
				diff = _ScalarProduct(theList[i], theSeg->extremity1,
					theList[i], theSeg->extremity2, dim);
				if (maxdiff < diff) {
					index = i;
					maxdiff = diff;
				}
			}
			return(index);

			/* END
			NO REDUCTION CASE
			*/
			break;

		default:
		case 1:

			/* REDUCTION INSIDE THE SMALLEST SPHERE
			*/

			/* AB = diameter extremities
			MA.MB = (MC+CA).(MC+CB) = MC^2 + CA.CB + MC ( CB+CA )
			= MC^2 - R^2   + 0
			*/
			if (dim == 2) {
				for (i = first; i <= l; i++) {
					diff = _ScalarProduct2D(theList[i], theSeg->extremity1,
						theList[i], theSeg->extremity2);
					if (diff > maxdiff) {
						index = i;
						maxdiff = diff;
					}
					else if (diff <= threshold) {
						_SwapPoints(theList, i, l);
						i--;   l--;
						continue;
					}
				}
				*last = l;
				return(index);
			}


			if (dim == 3) {
				for (i = first; i <= l; i++) {
					diff = _ScalarProduct3D(theList[i], theSeg->extremity1,
						theList[i], theSeg->extremity2);
					if (maxdiff < diff) {
						index = i;
						maxdiff = diff;
					}
					else if (diff <= threshold) {
						_SwapPoints(theList, i, l);
						i--;   l--;
						continue;
					}
				}
				*last = l;
				return(index);
			}


			for (i = first; i <= l; i++) {
				diff = _ScalarProduct(theList[i], theSeg->extremity1,
					theList[i], theSeg->extremity2, dim);
				if (maxdiff < diff) {
					index = i;
					maxdiff = diff;
				}
				else if (diff <= threshold) {
					_SwapPoints(theList, i, l);
					i--;   l--;
					continue;
				}
			}
			*last = l;
			return(index);
			/* END
			REDUCTION INSIDE THE SMALLEST SPHERE
			*/
		}

		return(-1);
	}













	double _MaximalSegmentInTwoLists(typeSegment *theSeg,
		const int index1,
		double **theList1,
		int *first1,
		int *last1,
		double **theList2,
		int *first2,
		int *last2,
		int dim)
	{

		int f1 = *first1;
		int l1 = *last1;

		int i1 = index1;

		int f2 = *first2;
		int l2 = *last2;

		int i2;

		double *ref;

		double d, dprevious;


		theSeg->extremity1 = (double*)NULL;
		theSeg->extremity2 = (double*)NULL;
		theSeg->squareDiameter = std::numeric_limits<double>::lowest();


		if (*first1 < 0 || *last1 < 0) return(-1.0);
		if (*first1 > *last1) return(0.0);
		if (*first2 < 0 || *last2 < 0) return(-1.0);
		if (*first2 > *last2) return(0.0);
		if (*first2 > *last2) {
			l2 = *first2;
			f2 = *last2;
		}
		if (index1 < f1 || index1 > l1) return(-1.0);


		do {

			dprevious = theSeg->squareDiameter;

			/* reference point in list #1
			*/
			ref = theList1[i1];
			/* point #i1 will be compared against all other points
			in list #2
			reject it at the beginning of the list
			do not consider it in the future
			*/
			_SwapPoints(theList1, f1, i1);
			f1++;

			/* find the furthest point from 'ref'
			*/
			d = _MaximalDistanceFromPoint(&i2, ref, theList2, f2, l2, dim);

			/* better estimate
			*/
			if (d > theSeg->squareDiameter) {

				theSeg->extremity1 = ref;
				theSeg->extremity2 = theList2[i2];
				theSeg->squareDiameter = d;

				if (f1 <= l1) {

					dprevious = theSeg->squareDiameter;

					/* reference point in list #1
					*/
					ref = theList2[i2];
					/* point #i2 will be compared against all other points
					in list #1
					reject it at the beginning of the list
					do not consider it in the future
					*/
					_SwapPoints(theList2, f2, i2);
					f2++;

					/* find the furthest point from 'ref'
					*/
					d = _MaximalDistanceFromPoint(&i1, ref, theList1, f1, l1, dim);

					/* better estimate
					*/
					if (d > theSeg->squareDiameter) {

						theSeg->extremity1 = theList1[i1];
						theSeg->extremity2 = ref;
						theSeg->squareDiameter = d;

					}

				}

			}


		} while (theSeg->squareDiameter > dprevious && f1 <= l1 && f2 <= l2);


		*first1 = f1;
		*last1 = l1;
		*first2 = f2;
		*last2 = l2;
		return(theSeg->squareDiameter);

	}












	double _MaximalSegmentInOneList(typeSegment *theSeg,
		const int index,
		double **theList,
		int *first,
		int *last,
		const int dim)
	{
		int f = *first;
		int l = *last;

		int i = index;
		double *ref;

		double d, dprevious;

		theSeg->extremity1 = (double*)NULL;
		theSeg->extremity2 = (double*)NULL;
		theSeg->squareDiameter = std::numeric_limits<double>::lowest();



		if (*first < 0 || *last < 0) return(-1.0);
		if (*first > *last) return(0.0);
		if (index < f || index > l) return(-1.0);
		if (f == l) {
			theSeg->extremity1 = theList[i];
			theSeg->extremity2 = theList[i];
			return(0.0);
		}



		do {

			dprevious = theSeg->squareDiameter;


			ref = theList[i];
			/* point #i will be compared against all other points
			reject it at the beginning of the list
			do not consider it in the future
			*/
			_SwapPoints(theList, f, i);
			f++;


			/* find the furthest point from 'ref'
			*/
			d = _MaximalDistanceFromPoint(&i, ref, theList, f, l, dim);
			if (d > theSeg->squareDiameter) {
				theSeg->extremity1 = ref;
				theSeg->extremity2 = theList[i];
				theSeg->squareDiameter = d;
			}


		} while (theSeg->squareDiameter > dprevious && f <= l);

		*first = f;
		*last = l;
		return(theSeg->squareDiameter);
	}








	double _MaximalDistanceFromPoint(int *index,
		const double *ref,
		double **theList,
		const int first,
		const int last,
		const int dim)
	{
		int f = first;
		int l = last;
		int i;
		double dmax, d;

		*index = -1;
		if (first < 0 || last < 0) return(-1.0);
		if (first > last) return(0.0);



		if (dim == 2) {

			dmax = _SquareDistance2D(theList[f], ref);
			*index = f;

			for (i = f + 1; i <= l; i++) {
				d = _SquareDistance2D(theList[i], ref);
				if (d > dmax) {
					dmax = d;
					*index = i;
				}
			}
			return(dmax);

		}

		if (dim == 3) {

			dmax = _SquareDistance3D(theList[f], ref);
			*index = f;

			for (i = f + 1; i <= l; i++) {
				d = _SquareDistance3D(theList[i], ref);
				if (d > dmax) {
					dmax = d;
					*index = i;
				}
			}
			return(dmax);

		}




		dmax = _SquareDistance(theList[f], ref, dim);
		*index = f;

		for (i = f + 1; i <= l; i++) {
			d = _SquareDistance(theList[i], ref, dim);
			if (d > dmax) {
				dmax = d;
				*index = i;
			}
		}
		return(dmax);

	}








	/* Swap two points
	*/

	void _SwapPoints(double **theList, const int i, const int j)
	{
		double *tmp;
		tmp = theList[i];
		theList[i] = theList[j];
		theList[j] = tmp;
	}





	static int _base_ = 100000000;

	void _InitCounter(typeCounter *c)
	{
		c->c1 = c->c2 = 0;
	}
	void _AddToCounter(typeCounter *c, const int i)
	{
		c->c2 += i;
		while (c->c2 >= _base_) {
			c->c2 -= _base_;
			c->c1++;
		}
		while (c->c2 < 0) {
			c->c2 += _base_;
			c->c1--;
		}
	}
	double _GetCounterAverage(typeCounter *c, const int i)
	{
		return((_base_ / (double)i) * c->c1 + c->c2 / (double)i);
	}


#ifdef _STATS_
	typeCounter scalarProducts;

	void _InitScalarProductCounter()
	{
		_InitCounter(&scalarProducts);
	}
	static void _IncScalarProductCounter()
	{
		_AddToCounter(&scalarProducts, 1);
	}
	double _GetScalarProductAverage(int n)
	{
		return(_GetCounterAverage(&scalarProducts, n));
	}
#endif




	/* square distance
	*/
	double _SquareDistance(const double *a, const double *b, const int dim)
	{
		register int i;
		register double d = 0.0;
		register double ba;
		for (i = 0; i<dim; i++) {
			ba = b[i] - a[i];
			d += ba * ba;
		}

#ifdef _STATS_
		_IncScalarProductCounter();
#endif

		return(d);
	}





	double _SquareDistance3D(const double *a, const double *b)
	{

#ifdef _STATS_
		_IncScalarProductCounter();
#endif

		return((b[0] - a[0])*(b[0] - a[0]) +
			(b[1] - a[1])*(b[1] - a[1]) +
			(b[2] - a[2])*(b[2] - a[2]));
	}


	double _SquareDistance2D(const double *a, const double *b)
	{

#ifdef _STATS_
		_IncScalarProductCounter();
#endif

		return((b[0] - a[0])*(b[0] - a[0]) +
			(b[1] - a[1])*(b[1] - a[1]));
	}





	/* dot product
	*/
	double _ScalarProduct(const double *a, const double *b,
		const double *c, const double *d, const int dim)
	{
		register int i;
		register double scalar = 0.0;
		register double ab, cd;
		for (i = 0; i<dim; i++) {
			ab = b[i] - a[i];
			cd = d[i] - c[i];
			scalar += ab * cd;
		}

#ifdef _STATS_
		_IncScalarProductCounter();
#endif

		return(scalar);
	}



	double _ScalarProduct3D(const double *a, const double *b,
		const double *c, const double *d)
	{

#ifdef _STATS_
		_IncScalarProductCounter();
#endif

		return((b[0] - a[0])*(d[0] - c[0]) +
			(b[1] - a[1])*(d[1] - c[1]) +
			(b[2] - a[2])*(d[2] - c[2]));
	}

	double _ScalarProduct2D(const double *a, const double *b,
		const double *c, const double *d)
	{

#ifdef _STATS_
		_IncScalarProductCounter();
#endif

		return((b[0] - a[0])*(d[0] - c[0]) +
			(b[1] - a[1])*(d[1] - c[1]));
	}






	int _FindPointInList(double **theList,
		const int first,
		const int last,
		double x0,
		double x1)
	{
		int i, j = -1;
		double e = 1e-5;
		for (i = first; i <= last; i++) {
			if (theList[i][0] - e < x0 && theList[i][0] + e > x0 &&
				theList[i][1] - e < x1 && theList[i][1] + e > x1) {
				if (j == -1) {
					j = i;
				}
				else {
					fprintf(stderr, "found again at #%d\n", i);
				}
			}
		}
		return(j);
	}






















	double _QuadraticDiameterInOneList(typeSegment *theDiam,
		double **theList,
		const int first,
		const int last,
		const int dim)
	{
		int i, j;
		double d;

		theDiam->extremity1 = (double*)NULL;
		theDiam->extremity2 = (double*)NULL;
		theDiam->squareDiameter = std::numeric_limits<double>::lowest();

		d = 0.0;

		if (dim == 2) {

			for (i = first; i <= last - 1; i++)
				for (j = i + 1; j <= last; j++) {
				d = _SquareDistance2D(theList[i], theList[j]);
				if (d > theDiam->squareDiameter) {
					theDiam->squareDiameter = d;
					theDiam->extremity1 = theList[i];
					theDiam->extremity2 = theList[j];
				}
				}
			return(theDiam->squareDiameter);

		}



		if (dim == 3) {

			for (i = first; i <= last - 1; i++)
				for (j = i + 1; j <= last; j++) {
				d = _SquareDistance3D(theList[i], theList[j]);
				if (d > theDiam->squareDiameter) {
					theDiam->squareDiameter = d;
					theDiam->extremity1 = theList[i];
					theDiam->extremity2 = theList[j];
				}
				}
			return(theDiam->squareDiameter);

		}



		for (i = first; i <= last - 1; i++)
			for (j = i + 1; j <= last; j++) {
			d = _SquareDistance(theList[i], theList[j], dim);
			if (d > theDiam->squareDiameter) {
				theDiam->squareDiameter = d;
				theDiam->extremity1 = theList[i];
				theDiam->extremity2 = theList[j];
			}
			}

		return(theDiam->squareDiameter);
	}








	double _QuadraticDiameterInTwoLists(typeSegment *theDiam,
		int   *index1,
		int   *index2,
		double **theList1,
		const int first1,
		const int last1,
		double **theList2,
		const int first2,
		const int last2,
		const int dim)
	{
		int i, j;
		double d;

		/*
		theDiam->extremity1 = (double*)NULL;
		theDiam->extremity2 = (double*)NULL;
		theDiam->squareDiameter = std::numeric_limits<double>::lowest();
		*/

		if (index1 != NULL && index2 != NULL) {
			*index1 = first1 - 1;
			*index2 = first2 - 1;
		}

		d = 0.0;

		if (dim == 2) {

			for (i = first1; i <= last1; i++)
				for (j = first2; j <= last2; j++) {
				d = _SquareDistance2D(theList1[i], theList2[j]);
				if (d > theDiam->squareDiameter) {
					theDiam->squareDiameter = d;
					theDiam->extremity1 = theList1[i];
					theDiam->extremity2 = theList2[j];
					if (index1 != NULL && index2 != NULL) {
						*index1 = i;
						*index2 = j;
					}
				}
				}
			return(theDiam->squareDiameter);

		}



		if (dim == 3) {

			for (i = first1; i <= last1; i++)
				for (j = first2; j <= last2; j++) {
				d = _SquareDistance3D(theList1[i], theList2[j]);
				if (d > theDiam->squareDiameter) {
					theDiam->squareDiameter = d;
					theDiam->extremity1 = theList1[i];
					theDiam->extremity2 = theList2[j];
					if (index1 != NULL && index2 != NULL) {
						*index1 = i;
						*index2 = j;
					}
				}
				}
			return(theDiam->squareDiameter);

		}



		for (i = first1; i <= last1; i++)
			for (j = first2; j <= last2; j++) {
			d = _SquareDistance(theList1[i], theList2[j], dim);
			if (d > theDiam->squareDiameter) {
				theDiam->squareDiameter = d;
				theDiam->extremity1 = theList1[i];
				theDiam->extremity2 = theList2[j];
				if (index1 != NULL && index2 != NULL) {
					*index1 = i;
					*index2 = j;
				}
			}
			}

		return(theDiam->squareDiameter);
	}

	/* === */

	static long int _random_calls_ = 0;
	static long int _random_seed_ = 0;

	long int _GetRandomCalls()
	{
		return(_random_calls_);
	}

	long int _GetRandomSeed()
	{
		return(_random_seed_);
	}



#ifdef WIN32
	void _SetRandomSeed(unsigned int seed)
	{
		srand(seed);
		_random_seed_ = seed;
		_random_calls_ = 0;
	}
#else
	void _SetRandomSeed(long int seed)
	{
		srand48(seed);
		_random_seed_ = seed;
		_random_calls_ = 0;
	}
#endif



#ifdef WIN32
	double _GetRandomDoubleNb()
	{
		_random_calls_++;
		return((double)rand() / (double)RAND_MAX);
	}
#else
	double _GetRandomDoubleNb()
	{
		_random_calls_++;
		return(drand48());
	}
#endif


	int _GetRandomIntNb(int min, int max)
	{
		if (min <= max)
			return((int)(floor(min + _GetRandomDoubleNb()*(double)(max - min + 1.0))));
		return((int)(floor(max + _GetRandomDoubleNb()*(double)(min - max + 1.0))));
	}


	double estimateDiameterInOneList(typeSegment *theDiam,
		double **theList,
		const int first,
		const int last,
		const int dim,
		double _epsilon_)
	{
		int index;

		int f = first;
		int l = last;

		int newEstimateIsSmallerThanCurrentEstimate;
		typeSegment theSeg;

		typeListOfSegments theDoubleNormals;

		double newEstimate;

		int newlast;


		int verboseWhenReducing = _GetVerboseWhenReducing();
		int _reduction_mode_in_iterative_ = _GetReductionModeInIterative();
		int tryToReduceQ = _GetTryToReduceQ();
		int _reduction_mode_of_diameter_ = _GetReductionModeOfDiameter();
		int _reduction_mode_of_dbleNorm_ = _GetReductionModeOfDbleNorm();
		int _Q_scan_ = _GetQscan();
		int _tight_bounds_ = _GetTightBounds();

		int i, j, k, n;
		int index1, index2;

		int suspicion_of_convex_hull = 0;
		int fdn, ldn, idn;

		double epsilon = _epsilon_;
		double bound, upperBound = 0.0;

		double upperSquareDiameter = 0.0;


		theDoubleNormals.n = 0;
		theDoubleNormals.nalloc = 0;
		theDoubleNormals.seg = NULL;




		theDiam->extremity1 = (double*)NULL;
		theDiam->extremity2 = (double*)NULL;
		theDiam->squareDiameter = std::numeric_limits<double>::lowest();


		if (first < 0 || last < 0) return(-1.0);
		if (first > last) {
			l = first;
			f = last;
		}
		if (f == l) {
			theDiam->extremity1 = theList[f];
			theDiam->extremity2 = theList[l];
			return(0.0);
		}



		index = _GetRandomIntNb(first, last);


		do {

			/* end conditions
			*/
			newEstimateIsSmallerThanCurrentEstimate = 0;

			/* find a double normal
			*/
			newEstimate = _MaximalSegmentInOneList(&theSeg, index, theList,
				&f, &l, dim);

			/* if we get a better estimation
			*/
			if (newEstimate > theDiam->squareDiameter) {

				/* update variables
				*/
				*theDiam = theSeg;

				/* keep the maximal segment in list
				*/
				if (_AddSegmentToList(&theSeg, &theDoubleNormals) != 1) {
					if (theDoubleNormals.nalloc > 0) free(theDoubleNormals.seg);
					return(-1.0);
				}


				/* find the farthest point outside the sphere
				*/
				newlast = l;
				index = _FarthestPointFromSphere(&theSeg, theList,
					f, &newlast, dim,
					_reduction_mode_in_iterative_);
				if (_reduction_mode_in_iterative_ == 1) {
					if (verboseWhenReducing)
						//fprintf( stdout, "...processing frth: remove %d points\n", l-newlast );
						if (newlast == l) {
						suspicion_of_convex_hull = 1;
						_reduction_mode_of_diameter_ = 0;
						_reduction_mode_of_dbleNorm_ = 0;
						}
					l = newlast;
				}




				/* stopping condition
				no point outside the sphere
				*/
				if (index < f) {
					if (theDoubleNormals.nalloc > 0) free(theDoubleNormals.seg);
					return(theDiam->squareDiameter);
				}



				/* other stopping condition

				the farthest point M outside the sphere
				is not that far away with a ball of diameter AB (and center C)

				we have MA.MB = MC^2 - |AB|^2 / 4
				leading to 4 * MC^2 = 4 * MA.MB + |AB|^2

				an upper bound of the diameter is then 2 MC
				thus (4 * MA.MB + |AB|^2) is a squared upper bound of the diameter
				*/

				bound = 4.0 * _ScalarProduct(theList[index], theDiam->extremity1,
					theList[index], theDiam->extremity2, dim) +
					theDiam->squareDiameter;

				/* stopping condition

				we want 2*MC < (1+epsilon) * d_estimate
				d_estimate = |AB|

				4 * MC^2                       < (1+epsilon)^2 * (d_estimate)^2
				4 * MA.MB + (d_estimate)^2     < (1+epsilon)^2 * (d_estimate)^2
				1 + 4 * MA.MB / (d_estimate)^2 < (1+epsilon)^2

				*/
				if (1.0 + 4.0 * _ScalarProduct(theList[index], theDiam->extremity1,
					theList[index], theDiam->extremity2, dim) /
					theDiam->squareDiameter < (1.0 + epsilon) * (1.0 + epsilon)) {
					if (theDoubleNormals.nalloc > 0) free(theDoubleNormals.seg);
					return(bound);
				}

			}
			else {

				newEstimateIsSmallerThanCurrentEstimate = 1;

				/*  I add the found segment to the list
				there is no evidence that it is a maximal segment for
				the initial set P (but it is for the current one
				ie the initial minus the points which have been removed),
				it may somehow help in case of reduction of the set
				of potential extremities


				The list of double normals is sorted (from 0 to n)
				by increasing square diameter:
				- the best diameter is added (again) at the end
				- we search the right place for this new element

				*/

				if (_AddSegmentToList(theDiam, &theDoubleNormals) != 1) {
					if (theDoubleNormals.nalloc > 0) free(theDoubleNormals.seg);
					return(-1.0);
				}

				for (n = theDoubleNormals.n - 2; n >= 0; n--) {
					if (n == 0) {
						theDoubleNormals.seg[n] = theSeg;
					}
					else {
						if (theSeg.squareDiameter <= theDoubleNormals.seg[n].squareDiameter &&
							theSeg.squareDiameter >  theDoubleNormals.seg[n - 1].squareDiameter) {
							theDoubleNormals.seg[n] = theSeg;
							n = -1;
						}
						else {
							theDoubleNormals.seg[n] = theDoubleNormals.seg[n - 1];
						}
					}
				}

			}

		} while (newEstimateIsSmallerThanCurrentEstimate == 0);










		/* last processing with the found diameter
		- points inside the smallest sphere of the
		diameter may have been already removed
		*/
		if (_reduction_mode_in_iterative_ > 0 && _reduction_mode_of_diameter_ == 1)
			_reduction_mode_of_diameter_ = 0;



		newlast = l;
		index = _LastPointOutsideSphereWithDiameter(theDiam, theDiam->squareDiameter,
			theList, f, &newlast, dim,
			_reduction_mode_of_diameter_);
		if (_reduction_mode_of_diameter_ == 1 ||
			_reduction_mode_of_diameter_ == 2) {
			if (verboseWhenReducing)
				//fprintf( stdout, "...processing diam: remove %d points\n", l-newlast );
				if (newlast == l) {
				suspicion_of_convex_hull = 1;
				_reduction_mode_of_dbleNorm_ = 0;
				}
			l = newlast;
		}


		/* in some (rare) case, the remaining points outside the largest
		sphere are removed while searching for a better diameter
		thus it is still an avantageous case.
		if any, we have
		#f       -> #index : points outside the sphere
		#index+1 -> #l     : points inside the sphere
		*/
		if (index < f) {
			if (theDoubleNormals.nalloc > 0) free(theDoubleNormals.seg);
			return(theDiam->squareDiameter);
		}









		/* do you have enough precision?
		*/
		index2 = index;
		upperSquareDiameter = theDiam->squareDiameter * (1 + epsilon) * (1 + epsilon);


		index1 = _LastPointOutsideSphereWithDiameter(theDiam, upperSquareDiameter,
			theList, f, &index2, dim, 0);
		/* there is no points outside the sphere of diameter d * (1 + epsilon)
		find the farthest one to get a better upper bound of the diameter
		*/
		if (index1 < f) {
			if (theDoubleNormals.nalloc > 0) free(theDoubleNormals.seg);

			upperBound = 4.0 * _ScalarProduct(theList[f], theDiam->extremity1,
				theList[f], theDiam->extremity2, dim) +
				theDiam->squareDiameter;

			for (k = f + 1; k <= index2; k++) {
				bound = 4.0 * _ScalarProduct(theList[k], theDiam->extremity1,
					theList[k], theDiam->extremity2, dim) +
					theDiam->squareDiameter;
				if (upperBound < bound) upperBound = bound;
			}
			return(upperBound);
		}




		/* get an upper bound of the diameter with points in [#index1+1 -> #index2]
		if there are points in this interval
		else the upper bound is simply d

		we have
		#f        -> #index1 : points outside the sphere of diameter d (1+epsilon)
		#index1+1 -> #index2 : points outside the sphere but inside the previous one
		#index2+1 -> #l      : points inside the sphere
		*/
		if (_tight_bounds_) {
			upperBound = theDiam->squareDiameter;
			if (index1 < index2) {
				for (k = index1 + 1; k <= index2; k++) {
					bound = 4.0 * _ScalarProduct(theList[k], theDiam->extremity1,
						theList[k], theDiam->extremity2, dim) +
						theDiam->squareDiameter;
					if (upperBound < bound) upperBound = bound;
				}
			}
		}
		else {
			upperBound = upperSquareDiameter;
		}









		/* to get some information on
		the points
		*/
		if (0) {
			for (n = theDoubleNormals.n - 1; n >= 0; n--) {
				_CountPointsInSpheres(&theDoubleNormals.seg[n], theDiam->squareDiameter,
					theList, f, l, dim);
			}
		}






		/* here we will reduce the set of potential extremities
		for the diameter,
		ie the set of points which are to be compared against all
		the other points

		right now, we have
		#f        -> #index1 : points outside the sphere of diameter d+epsilon
		#index1+1 -> #index2 : points outside the sphere of diameter d
		#index2+1 -> #l      : points inside the sphere

		we have a set of maximal segments
		theDoubleNormals.seg[ #i ] for #i from 0 to theDoubleNormals.n-1
		with
		theDoubleNormals.seg[ theDoubleNormals.n-1 ] == theDiam

		*/

		index = index1;
		/* right now, we have
		#f       -> #index : points outside the sphere of diameter d*(1+epsilon)
		#index+1 -> #l     : points inside the sphere of diameter d*(1+epsilon)
		*/

		if (tryToReduceQ && theDoubleNormals.n > 1) {



			for (k = 0; k < theDoubleNormals.n; k++)
				theDoubleNormals.seg[k].reduction_mode = _reduction_mode_of_dbleNorm_;



			switch (_Q_scan_) {
			default:
			case 0:
				/* backward
				*/
				ldn = 0;   fdn = theDoubleNormals.n - 2;   idn = -1;
				break;
			case 1:
				/* forward
				*/
				fdn = 0;   ldn = theDoubleNormals.n - 2;   idn = +1;
				break;
			}



			for (n = fdn; n != (ldn + idn) && index >= f; n += idn) {



				/* in [ #f #index ] find the points outside the sphere
				theDoubleNormals.seg[ n ]

				as a result
				#f   -> #i     are to be compared with all other points
				#i+1 -> #index are to be compared with a subset
				if this subset is empty, continue

				*/
				i = _LastPointOutsideSphereWithDiameter(&theDoubleNormals.seg[n],
					upperSquareDiameter,
					theList, f, &index, dim, 0);
				if (i >= index) continue;


				/* remise a jour de l'upper bound,
				a partir des points de [#i+1     -> #index]
				par rapport a la double normale courante

				Ce sont des points qui etaient candidats parce qu'au
				dela de la precision demandee, mais qui sont maintenant en
				deca de celle-ci pour la double normale courante,
				donc on reste en deca de la precision, par contre je ne sais
				s'il faut toujours mettre a jour la borne sup ...

				*/
				if (_tight_bounds_) {
					for (k = i + 1; k <= index; k++) {
						bound = 4.0 * _ScalarProduct(theList[k], theDoubleNormals.seg[n].extremity1,
							theList[k], theDoubleNormals.seg[n].extremity2, dim) +
							theDoubleNormals.seg[n].squareDiameter;
						if (upperBound < bound) upperBound = bound;
					}
				}



				/* in [ #index+1 #l ] find the points outside the sphere
				theDoubleNormals.seg[ n ]

				as a result
				#index+1 -> #j   are to be compared with the previous subset
				*/

				newlast = l;

				if (_tight_bounds_) {
					j = _LastPointOutsideSphereAndBoundWithDiameter(&theDoubleNormals.seg[n],
						upperSquareDiameter,
						theList, index + 1, &newlast, dim,
						theDoubleNormals.seg[n].reduction_mode,
						&bound);
					if (upperBound < bound) upperBound = bound;
				}
				else {
					j = _LastPointOutsideSphereWithDiameter(&theDoubleNormals.seg[n],
						upperSquareDiameter,
						theList, index + 1, &newlast, dim,
						theDoubleNormals.seg[n].reduction_mode);
				}

				if (theDoubleNormals.seg[n].reduction_mode == 1 ||
					theDoubleNormals.seg[n].reduction_mode == 2) {

					if (verboseWhenReducing)
						//fprintf( stdout, "...processing dbNR: remove %d points\n", l-newlast );
						if (newlast == l) {
						suspicion_of_convex_hull = 1;
						for (k = 0; k < theDoubleNormals.n; k++)
							theDoubleNormals.seg[k].reduction_mode = 0;
						}
					l = newlast;
				}

				if (j <= index) {
					index = i;
					continue;
				}



				/* right now
				#f       -> #i     : points to be compared with all other points
				#i+1     -> #index : points to be compared with the below set
				#index+1 -> #j     : points to be compared with the above set
				#j+1     -> #l     : remaining points

				*/


				theSeg.extremity1 = (double*)NULL;
				theSeg.extremity2 = (double*)NULL;
				theSeg.squareDiameter = 0.0;


				newEstimate = _QuadraticDiameterInTwoLists(&theSeg, NULL, NULL,
					theList, i + 1, index,
					theList, index + 1, j,
					dim);




				if (newEstimate > theDiam->squareDiameter) {
					/* update variables
					*/
					*theDiam = theSeg;
					if (upperBound < newEstimate) upperBound = newEstimate;
					if (upperSquareDiameter < newEstimate) upperSquareDiameter = newEstimate;
					/* we find a better estimate
					it is perhaps not the diameter, according that one
					diameter extremity can be in [ #f #i ]
					The question are :
					1. have we to look for a maximal segment with these two points
					or not ?
					-> Seems not necessary ...
					2. have we to consider this segment with the others ?
					-> yes
					but we can not reduce the whole set with it !!!
					*/
					theDoubleNormals.seg[n] = theSeg;
					theDoubleNormals.seg[n].reduction_mode = 0;
					n -= idn;
				}


				/* Q is now reduced
				*/
				index = i;

			}
		}










		if (theDoubleNormals.nalloc > 0) free(theDoubleNormals.seg);

		/* exhautive search

		comparison of points from #f to #index
		against all others points
		*/

		if (dim == 2) {
			for (i = f; i <= index; i++)
				for (j = i + 1; j <= l; j++) {
				newEstimate = _SquareDistance2D(theList[i], theList[j]);
				if (newEstimate > theDiam->squareDiameter) {
					theDiam->extremity1 = theList[i];
					theDiam->extremity2 = theList[j];
					theDiam->squareDiameter = newEstimate;
					if (newEstimate > upperBound) upperBound = newEstimate;
				}
				}
			return(upperBound);
		}


		if (dim == 3) {
			for (i = f; i <= index; i++)
				for (j = i + 1; j <= l; j++) {
				newEstimate = _SquareDistance3D(theList[i], theList[j]);
				if (newEstimate > theDiam->squareDiameter) {
					theDiam->extremity1 = theList[i];
					theDiam->extremity2 = theList[j];
					theDiam->squareDiameter = newEstimate;
					if (newEstimate > upperBound) upperBound = newEstimate;
				}
				}
			return(upperBound);
		}


		for (i = f; i <= index; i++)
			for (j = i + 1; j <= l; j++) {
			newEstimate = _SquareDistance(theList[i], theList[j], dim);
			if (newEstimate > theDiam->squareDiameter) {
				theDiam->extremity1 = theList[i];
				theDiam->extremity2 = theList[j];
				theDiam->squareDiameter = newEstimate;
				if (newEstimate > upperBound) upperBound = newEstimate;
			}
			}
		return(upperBound);
	}


	double estimateDiameter(typeSegment *theDiam,
		double **theList,
		const int first,
		const int last,
		const int dim,
		double epsilon)
	{

		_SetReductionModeInIterative(0);
		_SetReductionModeOfDiameter(0);
		_SetReductionModeOfDbleNorm(0);
		_DoNotTryToReduceQ();

		return estimateDiameterInOneList(theDiam, theList, first, last, dim, epsilon);
	}
};


namespace GeometryPredicates{
	/*****************************************************************************/
	/*                                                                           */
	/*  Routines for Arbitrary Precision Floating-point Arithmetic               */
	/*  and Fast Robust Geometric Predicates                                     */
	/*  (predicates.c)                                                           */
	/*                                                                           */
	/*  May 18, 1996                                                             */
	/*                                                                           */
	/*  Placed in the public domain by                                           */
	/*  Jonathan Richard Shewchuk                                                */
	/*  School of Computer Science                                               */
	/*  Carnegie Mellon University                                               */
	/*  5000 Forbes Avenue                                                       */
	/*  Pittsburgh, Pennsylvania  15213-3891                                     */
	/*  jrs@cs.cmu.edu                                                           */
	/*                                                                           */
	/*  This file contains C implementation of algorithms for exact addition     */
	/*    and multiplication of floating-point numbers, and predicates for       */
	/*    robustly performing the orientation and incircle tests used in         */
	/*    computational geometry.  The algorithms and underlying theory are      */
	/*    described in Jonathan Richard Shewchuk.  "Adaptive Precision Floating- */
	/*    Point Arithmetic and Fast Robust Geometric Predicates."  Technical     */
	/*    Report CMU-CS-96-140, School of Computer Science, Carnegie Mellon      */
	/*    University, Pittsburgh, Pennsylvania, May 1996.  (Submitted to         */
	/*    Discrete & Computational Geometry.)                                    */
	/*                                                                           */
	/*  This file, the paper listed above, and other information are available   */
	/*    from the Web page http://www.cs.cmu.edu/~quake/robust.html .           */
	/*                                                                           */
	/*****************************************************************************/

	/*****************************************************************************/
	/*                                                                           */
	/*  Using this code:                                                         */
	/*                                                                           */
	/*  First, read the short or long version of the paper (from the Web page    */
	/*    above).                                                                */
	/*                                                                           */
	/*  Be sure to call exactinit() once, before calling any of the arithmetic   */
	/*    functions or geometric predicates.  Also be sure to turn on the        */
	/*    optimizer when compiling this file.                                    */
	/*                                                                           */
	/*                                                                           */
	/*  Several geometric predicates are defined.  Their parameters are all      */
	/*    points.  Each point is an array of two or three floating-point         */
	/*    numbers.  The geometric predicates, described in the papers, are       */
	/*                                                                           */
	/*    orient2d(pa, pb, pc)                                                   */
	/*    orient2dfast(pa, pb, pc)                                               */
	/*    orient3d(pa, pb, pc, pd)                                               */
	/*    orient3dfast(pa, pb, pc, pd)                                           */
	/*    incircle(pa, pb, pc, pd)                                               */
	/*    incirclefast(pa, pb, pc, pd)                                           */
	/*    insphere(pa, pb, pc, pd, pe)                                           */
	/*    inspherefast(pa, pb, pc, pd, pe)                                       */
	/*                                                                           */
	/*  Those with suffix "fast" are approximate, non-robust versions.  Those    */
	/*    without the suffix are adaptive precision, robust versions.  There     */
	/*    are also versions with the suffices "exact" and "slow", which are      */
	/*    non-adaptive, exact arithmetic versions, which I use only for timings  */
	/*    in my arithmetic papers.                                               */
	/*                                                                           */
	/*                                                                           */
	/*  An expansion is represented by an array of floating-point numbers,       */
	/*    sorted from smallest to largest magnitude (possibly with interspersed  */
	/*    zeros).  The length of each expansion is stored as a separate integer, */
	/*    and each arithmetic function returns an integer which is the length    */
	/*    of the expansion it created.                                           */
	/*                                                                           */
	/*  Several arithmetic functions are defined.  Their parameters are          */
	/*                                                                           */
	/*    e, f           Input expansions                                        */
	/*    elen, flen     Lengths of input expansions (must be >= 1)              */
	/*    h              Output expansion                                        */
	/*    b              Input scalar                                            */
	/*                                                                           */
	/*  The arithmetic functions are                                             */
	/*                                                                           */
	/*    grow_expansion(elen, e, b, h)                                          */
	/*    grow_expansion_zeroelim(elen, e, b, h)                                 */
	/*    expansion_sum(elen, e, flen, f, h)                                     */
	/*    expansion_sum_zeroelim1(elen, e, flen, f, h)                           */
	/*    expansion_sum_zeroelim2(elen, e, flen, f, h)                           */
	/*    fast_expansion_sum(elen, e, flen, f, h)                                */
	/*    fast_expansion_sum_zeroelim(elen, e, flen, f, h)                       */
	/*    linear_expansion_sum(elen, e, flen, f, h)                              */
	/*    linear_expansion_sum_zeroelim(elen, e, flen, f, h)                     */
	/*    scale_expansion(elen, e, b, h)                                         */
	/*    scale_expansion_zeroelim(elen, e, b, h)                                */
	/*    compress(elen, e, h)                                                   */
	/*                                                                           */
	/*  All of these are described in the long version of the paper; some are    */
	/*    described in the short version.  All return an integer that is the     */
	/*    length of h.  Those with suffix _zeroelim perform zero elimination,    */
	/*    and are recommended over their counterparts.  The procedure            */
	/*    fast_expansion_sum_zeroelim() (or linear_expansion_sum_zeroelim() on   */
	/*    processors that do not use the round-to-even tiebreaking rule) is      */
	/*    recommended over expansion_sum_zeroelim().  Each procedure has a       */
	/*    little note next to it (in the code below) that tells you whether or   */
	/*    not the output expansion may be the same array as one of the input     */
	/*    expansions.                                                            */
	/*                                                                           */
	/*                                                                           */
	/*  If you look around below, you'll also find macros for a bunch of         */
	/*    simple unrolled arithmetic operations, and procedures for printing     */
	/*    expansions (commented out because they don't work with all C           */
	/*    compilers) and for generating random floating-point numbers whose      */
	/*    significand bits are all random.  Most of the macros have undocumented */
	/*    requirements that certain of their parameters should not be the same   */
	/*    variable; for safety, better to make sure all the parameters are       */
	/*    distinct variables.  Feel free to send email to jrs@cs.cmu.edu if you  */
	/*    have questions.                                                        */
	/*                                                                           */
	/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef CPU86
#include <float.h>
#endif /* CPU86 */
#ifdef LINUX
#include <fpu_control.h>
#endif /* LINUX */

	typedef double REAL;

#ifdef USE_CGAL_PREDICATES
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
	typedef CGAL::Exact_predicates_inexact_constructions_kernel cgalEpick;
	typedef cgalEpick::Point_3 Point;
	cgalEpick cgal_pred_obj;
#endif // #ifdef USE_CGAL_PREDICATES

	/* On some machines, the exact arithmetic routines might be defeated by the  */
	/*   use of internal extended precision floating-point registers.  Sometimes */
	/*   this problem can be fixed by defining certain values to be volatile,    */
	/*   thus forcing them to be stored to memory and rounded off.  This isn't   */
	/*   a great solution, though, as it slows the arithmetic down.              */
	/*                                                                           */
	/* To try this out, write "#define INEXACT volatile" below.  Normally,       */
	/*   however, INEXACT should be defined to be nothing.  ("#define INEXACT".) */

#define INEXACT                          /* Nothing */
	/* #define INEXACT volatile */

	/* #define REAL double */                      /* float or double */
#define REALPRINT doubleprint
#define REALRAND doublerand
#define NARROWRAND narrowdoublerand
#define UNIFORMRAND uniformdoublerand

	/* Which of the following two methods of finding the absolute values is      */
	/*   fastest is compiler-dependent.  A few compilers can inline and optimize */
	/*   the fabs() call; but most will incur the overhead of a function call,   */
	/*   which is disastrously slow.  A faster way on IEEE machines might be to  */
	/*   mask the appropriate bit, but that's difficult to do in C.              */

	//#define Absolute(a)  ((a) >= 0.0 ? (a) : -(a))
#define Absolute(a)  fabs(a)

	/* Many of the operations are broken up into two pieces, a main part that    */
	/*   performs an approximate operation, and a "tail" that computes the       */
	/*   roundoff error of that operation.                                       */
	/*                                                                           */
	/* The operations Fast_Two_Sum(), Fast_Two_Diff(), Two_Sum(), Two_Diff(),    */
	/*   Split(), and Two_Product() are all implemented as described in the      */
	/*   reference.  Each of these macros requires certain variables to be       */
	/*   defined in the calling routine.  The variables `bvirt', `c', `abig',    */
	/*   `_i', `_j', `_k', `_l', `_m', and `_n' are declared `INEXACT' because   */
	/*   they store the result of an operation that may incur roundoff error.    */
	/*   The input parameter `x' (or the highest numbered `x_' parameter) must   */
	/*   also be declared `INEXACT'.                                             */

#define Fast_Two_Sum_Tail(a, b, x, y) \
  bvirt = x - a; \
  y = b - bvirt

#define Fast_Two_Sum(a, b, x, y) \
  x = (REAL) (a + b); \
  Fast_Two_Sum_Tail(a, b, x, y)

#define Fast_Two_Diff_Tail(a, b, x, y) \
  bvirt = a - x; \
  y = bvirt - b

#define Fast_Two_Diff(a, b, x, y) \
  x = (REAL) (a - b); \
  Fast_Two_Diff_Tail(a, b, x, y)

#define Two_Sum_Tail(a, b, x, y) \
  bvirt = (REAL) (x - a); \
  avirt = x - bvirt; \
  bround = b - bvirt; \
  around = a - avirt; \
  y = around + bround

#define Two_Sum(a, b, x, y) \
  x = (REAL) (a + b); \
  Two_Sum_Tail(a, b, x, y)

#define Two_Diff_Tail(a, b, x, y) \
  bvirt = (REAL) (a - x); \
  avirt = x + bvirt; \
  bround = bvirt - b; \
  around = a - avirt; \
  y = around + bround

#define Two_Diff(a, b, x, y) \
  x = (REAL) (a - b); \
  Two_Diff_Tail(a, b, x, y)

#define Split(a, ahi, alo) \
  c = (REAL) (splitter * a); \
  abig = (REAL) (c - a); \
  ahi = c - abig; \
  alo = a - ahi

#define Two_Product_Tail(a, b, x, y) \
  Split(a, ahi, alo); \
  Split(b, bhi, blo); \
  err1 = x - (ahi * bhi); \
  err2 = err1 - (alo * bhi); \
  err3 = err2 - (ahi * blo); \
  y = (alo * blo) - err3

#define Two_Product(a, b, x, y) \
  x = (REAL) (a * b); \
  Two_Product_Tail(a, b, x, y)

	/* Two_Product_Presplit() is Two_Product() where one of the inputs has       */
	/*   already been split.  Avoids redundant splitting.                        */

#define Two_Product_Presplit(a, b, bhi, blo, x, y) \
  x = (REAL) (a * b); \
  Split(a, ahi, alo); \
  err1 = x - (ahi * bhi); \
  err2 = err1 - (alo * bhi); \
  err3 = err2 - (ahi * blo); \
  y = (alo * blo) - err3

	/* Two_Product_2Presplit() is Two_Product() where both of the inputs have    */
	/*   already been split.  Avoids redundant splitting.                        */

#define Two_Product_2Presplit(a, ahi, alo, b, bhi, blo, x, y) \
  x = (REAL) (a * b); \
  err1 = x - (ahi * bhi); \
  err2 = err1 - (alo * bhi); \
  err3 = err2 - (ahi * blo); \
  y = (alo * blo) - err3

	/* Square() can be done more quickly than Two_Product().                     */

#define Square_Tail(a, x, y) \
  Split(a, ahi, alo); \
  err1 = x - (ahi * ahi); \
  err3 = err1 - ((ahi + ahi) * alo); \
  y = (alo * alo) - err3

#define Square(a, x, y) \
  x = (REAL) (a * a); \
  Square_Tail(a, x, y)

	/* Macros for summing expansions of various fixed lengths.  These are all    */
	/*   unrolled versions of Expansion_Sum().                                   */

#define Two_One_Sum(a1, a0, b, x2, x1, x0) \
  Two_Sum(a0, b , _i, x0); \
  Two_Sum(a1, _i, x2, x1)

#define Two_One_Diff(a1, a0, b, x2, x1, x0) \
  Two_Diff(a0, b , _i, x0); \
  Two_Sum( a1, _i, x2, x1)

#define Two_Two_Sum(a1, a0, b1, b0, x3, x2, x1, x0) \
  Two_One_Sum(a1, a0, b0, _j, _0, x0); \
  Two_One_Sum(_j, _0, b1, x3, x2, x1)

#define Two_Two_Diff(a1, a0, b1, b0, x3, x2, x1, x0) \
  Two_One_Diff(a1, a0, b0, _j, _0, x0); \
  Two_One_Diff(_j, _0, b1, x3, x2, x1)

#define Four_One_Sum(a3, a2, a1, a0, b, x4, x3, x2, x1, x0) \
  Two_One_Sum(a1, a0, b , _j, x1, x0); \
  Two_One_Sum(a3, a2, _j, x4, x3, x2)

#define Four_Two_Sum(a3, a2, a1, a0, b1, b0, x5, x4, x3, x2, x1, x0) \
  Four_One_Sum(a3, a2, a1, a0, b0, _k, _2, _1, _0, x0); \
  Four_One_Sum(_k, _2, _1, _0, b1, x5, x4, x3, x2, x1)

#define Four_Four_Sum(a3, a2, a1, a0, b4, b3, b1, b0, x7, x6, x5, x4, x3, x2, \
                      x1, x0) \
  Four_Two_Sum(a3, a2, a1, a0, b1, b0, _l, _2, _1, _0, x1, x0); \
  Four_Two_Sum(_l, _2, _1, _0, b4, b3, x7, x6, x5, x4, x3, x2)

#define Eight_One_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b, x8, x7, x6, x5, x4, \
                      x3, x2, x1, x0) \
  Four_One_Sum(a3, a2, a1, a0, b , _j, x3, x2, x1, x0); \
  Four_One_Sum(a7, a6, a5, a4, _j, x8, x7, x6, x5, x4)

#define Eight_Two_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0, x9, x8, x7, \
                      x6, x5, x4, x3, x2, x1, x0) \
  Eight_One_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b0, _k, _6, _5, _4, _3, _2, \
                _1, _0, x0); \
  Eight_One_Sum(_k, _6, _5, _4, _3, _2, _1, _0, b1, x9, x8, x7, x6, x5, x4, \
                x3, x2, x1)

#define Eight_Four_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b4, b3, b1, b0, x11, \
                       x10, x9, x8, x7, x6, x5, x4, x3, x2, x1, x0) \
  Eight_Two_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0, _l, _6, _5, _4, _3, \
                _2, _1, _0, x1, x0); \
  Eight_Two_Sum(_l, _6, _5, _4, _3, _2, _1, _0, b4, b3, x11, x10, x9, x8, \
                x7, x6, x5, x4, x3, x2)

	/* Macros for multiplying expansions of various fixed lengths.               */

#define Two_One_Product(a1, a0, b, x3, x2, x1, x0) \
  Split(b, bhi, blo); \
  Two_Product_Presplit(a0, b, bhi, blo, _i, x0); \
  Two_Product_Presplit(a1, b, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _k, x1); \
  Fast_Two_Sum(_j, _k, x3, x2)

#define Four_One_Product(a3, a2, a1, a0, b, x7, x6, x5, x4, x3, x2, x1, x0) \
  Split(b, bhi, blo); \
  Two_Product_Presplit(a0, b, bhi, blo, _i, x0); \
  Two_Product_Presplit(a1, b, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _k, x1); \
  Fast_Two_Sum(_j, _k, _i, x2); \
  Two_Product_Presplit(a2, b, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _k, x3); \
  Fast_Two_Sum(_j, _k, _i, x4); \
  Two_Product_Presplit(a3, b, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _k, x5); \
  Fast_Two_Sum(_j, _k, x7, x6)

#define Two_Two_Product(a1, a0, b1, b0, x7, x6, x5, x4, x3, x2, x1, x0) \
  Split(a0, a0hi, a0lo); \
  Split(b0, bhi, blo); \
  Two_Product_2Presplit(a0, a0hi, a0lo, b0, bhi, blo, _i, x0); \
  Split(a1, a1hi, a1lo); \
  Two_Product_2Presplit(a1, a1hi, a1lo, b0, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _k, _1); \
  Fast_Two_Sum(_j, _k, _l, _2); \
  Split(b1, bhi, blo); \
  Two_Product_2Presplit(a0, a0hi, a0lo, b1, bhi, blo, _i, _0); \
  Two_Sum(_1, _0, _k, x1); \
  Two_Sum(_2, _k, _j, _1); \
  Two_Sum(_l, _j, _m, _2); \
  Two_Product_2Presplit(a1, a1hi, a1lo, b1, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _n, _0); \
  Two_Sum(_1, _0, _i, x2); \
  Two_Sum(_2, _i, _k, _1); \
  Two_Sum(_m, _k, _l, _2); \
  Two_Sum(_j, _n, _k, _0); \
  Two_Sum(_1, _0, _j, x3); \
  Two_Sum(_2, _j, _i, _1); \
  Two_Sum(_l, _i, _m, _2); \
  Two_Sum(_1, _k, _i, x4); \
  Two_Sum(_2, _i, _k, x5); \
  Two_Sum(_m, _k, x7, x6)

	/* An expansion of length two can be squared more quickly than finding the   */
	/*   product of two different expansions of length two, and the result is    */
	/*   guaranteed to have no more than six (rather than eight) components.     */

#define Two_Square(a1, a0, x5, x4, x3, x2, x1, x0) \
  Square(a0, _j, x0); \
  _0 = a0 + a0; \
  Two_Product(a1, _0, _k, _1); \
  Two_One_Sum(_k, _1, _j, _l, _2, x1); \
  Square(a1, _j, _1); \
  Two_Two_Sum(_j, _1, _l, _2, x5, x4, x3, x2)

	/* splitter = 2^ceiling(p / 2) + 1.  Used to split floats in half.           */
	static REAL splitter;
	static REAL epsilon;         /* = 2^(-p).  Used to estimate roundoff errors. */
	/* A set of coefficients used to calculate maximum roundoff errors.          */
	static REAL resulterrbound;
	static REAL ccwerrboundA, ccwerrboundB, ccwerrboundC;
	static REAL o3derrboundA, o3derrboundB, o3derrboundC;
	static REAL iccerrboundA, iccerrboundB, iccerrboundC;
	static REAL isperrboundA, isperrboundB, isperrboundC;

	// Options to choose types of geometric computtaions. 
	// Added by H. Si, 2012-08-23.
	static int  _use_inexact_arith; // -X option.
	static int  _use_static_filter; // Default option, disable it by -X1

	// Static filters for orient3d() and insphere(). 
	// They are pre-calcualted and set in exactinit().
	// Added by H. Si, 2012-08-23.
	static REAL o3dstaticfilter;
	static REAL ispstaticfilter;



	// The following codes were part of "IEEE 754 floating-point test software"
	//          http://www.math.utah.edu/~beebe/software/ieee/
	// The original program was "fpinfo2.c".

	double fppow2(int n)
	{
		double x, power;
		x = (n < 0) ? ((double)1.0 / (double)2.0) : (double)2.0;
		n = (n < 0) ? -n : n;
		power = (double)1.0;
		while (n-- > 0)
			power *= x;
		return (power);
	}

#ifdef SINGLE

	float fstore(float x)
	{
		return (x);
	}

	int test_float(int verbose)
	{
		float x;
		int pass = 1;

		//(void)printf("float:\n");

		if (verbose) {
			(void)printf("  sizeof(float) = %2u\n", (unsigned int)sizeof(float));
#ifdef CPU86  // <float.h>
			(void)printf("  FLT_MANT_DIG = %2d\n", FLT_MANT_DIG);
#endif
		}

		x = (float)1.0;
		while (fstore((float)1.0 + x / (float)2.0) != (float)1.0)
			x /= (float)2.0;
		if (verbose)
			(void)printf("  machine epsilon = %13.5e  ", x);

		if (x == (float)fppow2(-23)) {
			if (verbose)
				(void)printf("[IEEE 754 32-bit macheps]\n");
		}
		else {
			(void)printf("[not IEEE 754 conformant] !!\n");
			pass = 0;
		}

		x = (float)1.0;
		while (fstore(x / (float)2.0) != (float)0.0)
			x /= (float)2.0;
		if (verbose)
			(void)printf("  smallest positive number =  %13.5e  ", x);

		if (x == (float)fppow2(-149)) {
			if (verbose)
				(void)printf("[smallest 32-bit subnormal]\n");
		}
		else if (x == (float)fppow2(-126)) {
			if (verbose)
				(void)printf("[smallest 32-bit normal]\n");
		}
		else {
			(void)printf("[not IEEE 754 conformant] !!\n");
			pass = 0;
		}

		return pass;
	}

# else

	double dstore(double x)
	{
		return (x);
	}

	int test_double(int verbose)
	{
		double x;
		int pass = 1;

		// (void)printf("double:\n");
		if (verbose) {
			(void)printf("  sizeof(double) = %2u\n", (unsigned int)sizeof(double));
#ifdef CPU86  // <float.h>
			(void)printf("  DBL_MANT_DIG = %2d\n", DBL_MANT_DIG);
#endif
		}

		x = 1.0;
		while (dstore(1.0 + x / 2.0) != 1.0)
			x /= 2.0;
		if (verbose)
			(void)printf("  machine epsilon = %13.5le ", x);

		if (x == (double)fppow2(-52)) {
			if (verbose)
				(void)printf("[IEEE 754 64-bit macheps]\n");
		}
		else {
			(void)printf("[not IEEE 754 conformant] !!\n");
			pass = 0;
		}

		x = 1.0;
		while (dstore(x / 2.0) != 0.0)
			x /= 2.0;
		//if (verbose)
		//  (void)printf("  smallest positive number = %13.5le ", x);

		if (x == (double)fppow2(-1074)) {
			//if (verbose)
			//  (void)printf("[smallest 64-bit subnormal]\n");
		}
		else if (x == (double)fppow2(-1022)) {
			//if (verbose)
			//  (void)printf("[smallest 64-bit normal]\n");
		}
		else {
			(void)printf("[not IEEE 754 conformant] !!\n");
			pass = 0;
		}

		return pass;
	}

#endif

	/*****************************************************************************/
	/*                                                                           */
	/*  exactinit()   Initialize the variables used for exact arithmetic.        */
	/*                                                                           */
	/*  `epsilon' is the largest power of two such that 1.0 + epsilon = 1.0 in   */
	/*  floating-point arithmetic.  `epsilon' bounds the relative roundoff       */
	/*  error.  It is used for floating-point error analysis.                    */
	/*                                                                           */
	/*  `splitter' is used to split floating-point numbers into two half-        */
	/*  length significands for exact multiplication.                            */
	/*                                                                           */
	/*  I imagine that a highly optimizing compiler might be too smart for its   */
	/*  own good, and somehow cause this routine to fail, if it pretends that    */
	/*  floating-point arithmetic is too much like real arithmetic.              */
	/*                                                                           */
	/*  Don't change this routine unless you fully understand it.                */
	/*                                                                           */
	/*****************************************************************************/

	void exactinit(int verbose, int noexact, int nofilter, REAL maxx, REAL maxy,
		REAL maxz)
	{
		REAL half;
		REAL check, lastcheck;
		int every_other;
#ifdef LINUX
		int cword;
#endif /* LINUX */

#ifdef CPU86
#ifdef SINGLE
		_control87(_PC_24, _MCW_PC); /* Set FPU control word for single precision. */
#else /* not SINGLE */
		_control87(_PC_53, _MCW_PC); /* Set FPU control word for double precision. */
#endif /* not SINGLE */
#endif /* CPU86 */
#ifdef LINUX
#ifdef SINGLE
		/*  cword = 4223; */
		cword = 4210;                 /* set FPU control word for single precision */
#else /* not SINGLE */
		/*  cword = 4735; */
		cword = 4722;                 /* set FPU control word for double precision */
#endif /* not SINGLE */
		_FPU_SETCW(cword);
#endif /* LINUX */

		if (verbose) {
			printf("  Initializing robust predicates.\n");
		}

#ifdef USE_CGAL_PREDICATES
		if (cgal_pred_obj.Has_static_filters) {
			printf("  Use static filter.\n");
		}
		else {
			printf("  No static filter.\n");
		}
#endif // USE_CGAL_PREDICATES

#ifdef SINGLE
		test_float(verbose);
#else
		test_double(verbose);
#endif

		every_other = 1;
		half = 0.5;
		epsilon = 1.0;
		splitter = 1.0;
		check = 1.0;
		/* Repeatedly divide `epsilon' by two until it is too small to add to    */
		/*   one without causing roundoff.  (Also check if the sum is equal to   */
		/*   the previous sum, for machines that round up instead of using exact */
		/*   rounding.  Not that this library will work on such machines anyway. */
		do {
			lastcheck = check;
			epsilon *= half;
			if (every_other) {
				splitter *= 2.0;
			}
			every_other = !every_other;
			check = 1.0 + epsilon;
		} while ((check != 1.0) && (check != lastcheck));
		splitter += 1.0;

		/* Error bounds for orientation and incircle tests. */
		resulterrbound = (3.0 + 8.0 * epsilon) * epsilon;
		ccwerrboundA = (3.0 + 16.0 * epsilon) * epsilon;
		ccwerrboundB = (2.0 + 12.0 * epsilon) * epsilon;
		ccwerrboundC = (9.0 + 64.0 * epsilon) * epsilon * epsilon;
		o3derrboundA = (7.0 + 56.0 * epsilon) * epsilon;
		o3derrboundB = (3.0 + 28.0 * epsilon) * epsilon;
		o3derrboundC = (26.0 + 288.0 * epsilon) * epsilon * epsilon;
		iccerrboundA = (10.0 + 96.0 * epsilon) * epsilon;
		iccerrboundB = (4.0 + 48.0 * epsilon) * epsilon;
		iccerrboundC = (44.0 + 576.0 * epsilon) * epsilon * epsilon;
		isperrboundA = (16.0 + 224.0 * epsilon) * epsilon;
		isperrboundB = (5.0 + 72.0 * epsilon) * epsilon;
		isperrboundC = (71.0 + 1408.0 * epsilon) * epsilon * epsilon;

		// Set TetGen options.  Added by H. Si, 2012-08-23.
		_use_inexact_arith = noexact;
		_use_static_filter = !nofilter;

		// Calculate the two static filters for orient3d() and insphere() tests.
		// Added by H. Si, 2012-08-23.

		// Sort maxx < maxy < maxz. Re-use 'half' for swapping.
		assert(maxx > 0);
		assert(maxy > 0);
		assert(maxz > 0);

		if (maxx > maxz) {
			half = maxx; maxx = maxz; maxz = half;
		}
		if (maxy > maxz) {
			half = maxy; maxy = maxz; maxz = half;
		}
		else if (maxy < maxx) {
			half = maxy; maxy = maxx; maxx = half;
		}

		o3dstaticfilter = 5.1107127829973299e-15 * maxx * maxy * maxz;
		ispstaticfilter = 1.2466136531027298e-13 * maxx * maxy * maxz * (maxz * maxz);

	}

	/*****************************************************************************/
	/*                                                                           */
	/*  grow_expansion()   Add a scalar to an expansion.                         */
	/*                                                                           */
	/*  Sets h = e + b.  See the long version of my paper for details.           */
	/*                                                                           */
	/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
	/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    */
	/*  properties as well.  (That is, if e has one of these properties, so      */
	/*  will h.)                                                                 */
	/*                                                                           */
	/*****************************************************************************/

	int grow_expansion(int elen, REAL *e, REAL b, REAL *h)
		/* e and h can be the same. */
	{
		REAL Q;
		INEXACT REAL Qnew;
		int eindex;
		REAL enow;
		INEXACT REAL bvirt;
		REAL avirt, bround, around;

		Q = b;
		for (eindex = 0; eindex < elen; eindex++) {
			enow = e[eindex];
			Two_Sum(Q, enow, Qnew, h[eindex]);
			Q = Qnew;
		}
		h[eindex] = Q;
		return eindex + 1;
	}

	/*****************************************************************************/
	/*                                                                           */
	/*  grow_expansion_zeroelim()   Add a scalar to an expansion, eliminating    */
	/*                              zero components from the output expansion.   */
	/*                                                                           */
	/*  Sets h = e + b.  See the long version of my paper for details.           */
	/*                                                                           */
	/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
	/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    */
	/*  properties as well.  (That is, if e has one of these properties, so      */
	/*  will h.)                                                                 */
	/*                                                                           */
	/*****************************************************************************/

	int grow_expansion_zeroelim(int elen, REAL *e, REAL b, REAL *h)
		/* e and h can be the same. */
	{
		REAL Q, hh;
		INEXACT REAL Qnew;
		int eindex, hindex;
		REAL enow;
		INEXACT REAL bvirt;
		REAL avirt, bround, around;

		hindex = 0;
		Q = b;
		for (eindex = 0; eindex < elen; eindex++) {
			enow = e[eindex];
			Two_Sum(Q, enow, Qnew, hh);
			Q = Qnew;
			if (hh != 0.0) {
				h[hindex++] = hh;
			}
		}
		if ((Q != 0.0) || (hindex == 0)) {
			h[hindex++] = Q;
		}
		return hindex;
	}

	/*****************************************************************************/
	/*                                                                           */
	/*  expansion_sum()   Sum two expansions.                                    */
	/*                                                                           */
	/*  Sets h = e + f.  See the long version of my paper for details.           */
	/*                                                                           */
	/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
	/*  with IEEE 754), maintains the nonadjacent property as well.  (That is,   */
	/*  if e has one of these properties, so will h.)  Does NOT maintain the     */
	/*  strongly nonoverlapping property.                                        */
	/*                                                                           */
	/*****************************************************************************/

	int expansion_sum(int elen, REAL *e, int flen, REAL *f, REAL *h)
		/* e and h can be the same, but f and h cannot. */
	{
		REAL Q;
		INEXACT REAL Qnew;
		int findex, hindex, hlast;
		REAL hnow;
		INEXACT REAL bvirt;
		REAL avirt, bround, around;

		Q = f[0];
		for (hindex = 0; hindex < elen; hindex++) {
			hnow = e[hindex];
			Two_Sum(Q, hnow, Qnew, h[hindex]);
			Q = Qnew;
		}
		h[hindex] = Q;
		hlast = hindex;
		for (findex = 1; findex < flen; findex++) {
			Q = f[findex];
			for (hindex = findex; hindex <= hlast; hindex++) {
				hnow = h[hindex];
				Two_Sum(Q, hnow, Qnew, h[hindex]);
				Q = Qnew;
			}
			h[++hlast] = Q;
		}
		return hlast + 1;
	}

	/*****************************************************************************/
	/*                                                                           */
	/*  expansion_sum_zeroelim1()   Sum two expansions, eliminating zero         */
	/*                              components from the output expansion.        */
	/*                                                                           */
	/*  Sets h = e + f.  See the long version of my paper for details.           */
	/*                                                                           */
	/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
	/*  with IEEE 754), maintains the nonadjacent property as well.  (That is,   */
	/*  if e has one of these properties, so will h.)  Does NOT maintain the     */
	/*  strongly nonoverlapping property.                                        */
	/*                                                                           */
	/*****************************************************************************/

	int expansion_sum_zeroelim1(int elen, REAL *e, int flen, REAL *f, REAL *h)
		/* e and h can be the same, but f and h cannot. */
	{
		REAL Q;
		INEXACT REAL Qnew;
		int index, findex, hindex, hlast;
		REAL hnow;
		INEXACT REAL bvirt;
		REAL avirt, bround, around;

		Q = f[0];
		for (hindex = 0; hindex < elen; hindex++) {
			hnow = e[hindex];
			Two_Sum(Q, hnow, Qnew, h[hindex]);
			Q = Qnew;
		}
		h[hindex] = Q;
		hlast = hindex;
		for (findex = 1; findex < flen; findex++) {
			Q = f[findex];
			for (hindex = findex; hindex <= hlast; hindex++) {
				hnow = h[hindex];
				Two_Sum(Q, hnow, Qnew, h[hindex]);
				Q = Qnew;
			}
			h[++hlast] = Q;
		}
		hindex = -1;
		for (index = 0; index <= hlast; index++) {
			hnow = h[index];
			if (hnow != 0.0) {
				h[++hindex] = hnow;
			}
		}
		if (hindex == -1) {
			return 1;
		}
		else {
			return hindex + 1;
		}
	}

	/*****************************************************************************/
	/*                                                                           */
	/*  expansion_sum_zeroelim2()   Sum two expansions, eliminating zero         */
	/*                              components from the output expansion.        */
	/*                                                                           */
	/*  Sets h = e + f.  See the long version of my paper for details.           */
	/*                                                                           */
	/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
	/*  with IEEE 754), maintains the nonadjacent property as well.  (That is,   */
	/*  if e has one of these properties, so will h.)  Does NOT maintain the     */
	/*  strongly nonoverlapping property.                                        */
	/*                                                                           */
	/*****************************************************************************/

	int expansion_sum_zeroelim2(int elen, REAL *e, int flen, REAL *f, REAL *h)
		/* e and h can be the same, but f and h cannot. */
	{
		REAL Q, hh;
		INEXACT REAL Qnew;
		int eindex, findex, hindex, hlast;
		REAL enow;
		INEXACT REAL bvirt;
		REAL avirt, bround, around;

		hindex = 0;
		Q = f[0];
		for (eindex = 0; eindex < elen; eindex++) {
			enow = e[eindex];
			Two_Sum(Q, enow, Qnew, hh);
			Q = Qnew;
			if (hh != 0.0) {
				h[hindex++] = hh;
			}
		}
		h[hindex] = Q;
		hlast = hindex;
		for (findex = 1; findex < flen; findex++) {
			hindex = 0;
			Q = f[findex];
			for (eindex = 0; eindex <= hlast; eindex++) {
				enow = h[eindex];
				Two_Sum(Q, enow, Qnew, hh);
				Q = Qnew;
				if (hh != 0) {
					h[hindex++] = hh;
				}
			}
			h[hindex] = Q;
			hlast = hindex;
		}
		return hlast + 1;
	}

	/*****************************************************************************/
	/*                                                                           */
	/*  fast_expansion_sum()   Sum two expansions.                               */
	/*                                                                           */
	/*  Sets h = e + f.  See the long version of my paper for details.           */
	/*                                                                           */
	/*  If round-to-even is used (as with IEEE 754), maintains the strongly      */
	/*  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h   */
	/*  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent      */
	/*  properties.                                                              */
	/*                                                                           */
	/*****************************************************************************/

	int fast_expansion_sum(int elen, REAL *e, int flen, REAL *f, REAL *h)
		/* h cannot be e or f. */
	{
		REAL Q;
		INEXACT REAL Qnew;
		INEXACT REAL bvirt;
		REAL avirt, bround, around;
		int eindex, findex, hindex;
		REAL enow, fnow;

		enow = e[0];
		fnow = f[0];
		eindex = findex = 0;
		if ((fnow > enow) == (fnow > -enow)) {
			Q = enow;
			enow = e[++eindex];
		}
		else {
			Q = fnow;
			fnow = f[++findex];
		}
		hindex = 0;
		if ((eindex < elen) && (findex < flen)) {
			if ((fnow > enow) == (fnow > -enow)) {
				Fast_Two_Sum(enow, Q, Qnew, h[0]);
				enow = e[++eindex];
			}
			else {
				Fast_Two_Sum(fnow, Q, Qnew, h[0]);
				fnow = f[++findex];
			}
			Q = Qnew;
			hindex = 1;
			while ((eindex < elen) && (findex < flen)) {
				if ((fnow > enow) == (fnow > -enow)) {
					Two_Sum(Q, enow, Qnew, h[hindex]);
					enow = e[++eindex];
				}
				else {
					Two_Sum(Q, fnow, Qnew, h[hindex]);
					fnow = f[++findex];
				}
				Q = Qnew;
				hindex++;
			}
		}
		while (eindex < elen) {
			Two_Sum(Q, enow, Qnew, h[hindex]);
			enow = e[++eindex];
			Q = Qnew;
			hindex++;
		}
		while (findex < flen) {
			Two_Sum(Q, fnow, Qnew, h[hindex]);
			fnow = f[++findex];
			Q = Qnew;
			hindex++;
		}
		h[hindex] = Q;
		return hindex + 1;
	}

	/*****************************************************************************/
	/*                                                                           */
	/*  fast_expansion_sum_zeroelim()   Sum two expansions, eliminating zero     */
	/*                                  components from the output expansion.    */
	/*                                                                           */
	/*  Sets h = e + f.  See the long version of my paper for details.           */
	/*                                                                           */
	/*  If round-to-even is used (as with IEEE 754), maintains the strongly      */
	/*  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h   */
	/*  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent      */
	/*  properties.                                                              */
	/*                                                                           */
	/*****************************************************************************/

	int fast_expansion_sum_zeroelim(int elen, REAL *e, int flen, REAL *f, REAL *h)
		/* h cannot be e or f. */
	{
		REAL Q;
		INEXACT REAL Qnew;
		INEXACT REAL hh;
		INEXACT REAL bvirt;
		REAL avirt, bround, around;
		int eindex, findex, hindex;
		REAL enow, fnow;

		enow = e[0];
		fnow = f[0];
		eindex = findex = 0;
		if ((fnow > enow) == (fnow > -enow)) {
			Q = enow;
			enow = e[++eindex];
		}
		else {
			Q = fnow;
			fnow = f[++findex];
		}
		hindex = 0;
		if ((eindex < elen) && (findex < flen)) {
			if ((fnow > enow) == (fnow > -enow)) {
				Fast_Two_Sum(enow, Q, Qnew, hh);
				enow = e[++eindex];
			}
			else {
				Fast_Two_Sum(fnow, Q, Qnew, hh);
				fnow = f[++findex];
			}
			Q = Qnew;
			if (hh != 0.0) {
				h[hindex++] = hh;
			}
			while ((eindex < elen) && (findex < flen)) {
				if ((fnow > enow) == (fnow > -enow)) {
					Two_Sum(Q, enow, Qnew, hh);
					enow = e[++eindex];
				}
				else {
					Two_Sum(Q, fnow, Qnew, hh);
					fnow = f[++findex];
				}
				Q = Qnew;
				if (hh != 0.0) {
					h[hindex++] = hh;
				}
			}
		}
		while (eindex < elen) {
			Two_Sum(Q, enow, Qnew, hh);
			enow = e[++eindex];
			Q = Qnew;
			if (hh != 0.0) {
				h[hindex++] = hh;
			}
		}
		while (findex < flen) {
			Two_Sum(Q, fnow, Qnew, hh);
			fnow = f[++findex];
			Q = Qnew;
			if (hh != 0.0) {
				h[hindex++] = hh;
			}
		}
		if ((Q != 0.0) || (hindex == 0)) {
			h[hindex++] = Q;
		}
		return hindex;
	}

	/*****************************************************************************/
	/*                                                                           */
	/*  linear_expansion_sum()   Sum two expansions.                             */
	/*                                                                           */
	/*  Sets h = e + f.  See either version of my paper for details.             */
	/*                                                                           */
	/*  Maintains the nonoverlapping property.  (That is, if e is                */
	/*  nonoverlapping, h will be also.)                                         */
	/*                                                                           */
	/*****************************************************************************/

	int linear_expansion_sum(int elen, REAL *e, int flen, REAL *f, REAL *h)
		/* h cannot be e or f. */
	{
		REAL Q, q;
		INEXACT REAL Qnew;
		INEXACT REAL R;
		INEXACT REAL bvirt;
		REAL avirt, bround, around;
		int eindex, findex, hindex;
		REAL enow, fnow;
		REAL g0;

		enow = e[0];
		fnow = f[0];
		eindex = findex = 0;
		if ((fnow > enow) == (fnow > -enow)) {
			g0 = enow;
			enow = e[++eindex];
		}
		else {
			g0 = fnow;
			fnow = f[++findex];
		}
		if ((eindex < elen) && ((findex >= flen)
			|| ((fnow > enow) == (fnow > -enow)))) {
			Fast_Two_Sum(enow, g0, Qnew, q);
			enow = e[++eindex];
		}
		else {
			Fast_Two_Sum(fnow, g0, Qnew, q);
			fnow = f[++findex];
		}
		Q = Qnew;
		for (hindex = 0; hindex < elen + flen - 2; hindex++) {
			if ((eindex < elen) && ((findex >= flen)
				|| ((fnow > enow) == (fnow > -enow)))) {
				Fast_Two_Sum(enow, q, R, h[hindex]);
				enow = e[++eindex];
			}
			else {
				Fast_Two_Sum(fnow, q, R, h[hindex]);
				fnow = f[++findex];
			}
			Two_Sum(Q, R, Qnew, q);
			Q = Qnew;
		}
		h[hindex] = q;
		h[hindex + 1] = Q;
		return hindex + 2;
	}

	/*****************************************************************************/
	/*                                                                           */
	/*  linear_expansion_sum_zeroelim()   Sum two expansions, eliminating zero   */
	/*                                    components from the output expansion.  */
	/*                                                                           */
	/*  Sets h = e + f.  See either version of my paper for details.             */
	/*                                                                           */
	/*  Maintains the nonoverlapping property.  (That is, if e is                */
	/*  nonoverlapping, h will be also.)                                         */
	/*                                                                           */
	/*****************************************************************************/

	int linear_expansion_sum_zeroelim(int elen, REAL *e, int flen, REAL *f,
		REAL *h)
		/* h cannot be e or f. */
	{
		REAL Q, q, hh;
		INEXACT REAL Qnew;
		INEXACT REAL R;
		INEXACT REAL bvirt;
		REAL avirt, bround, around;
		int eindex, findex, hindex;
		int count;
		REAL enow, fnow;
		REAL g0;

		enow = e[0];
		fnow = f[0];
		eindex = findex = 0;
		hindex = 0;
		if ((fnow > enow) == (fnow > -enow)) {
			g0 = enow;
			enow = e[++eindex];
		}
		else {
			g0 = fnow;
			fnow = f[++findex];
		}
		if ((eindex < elen) && ((findex >= flen)
			|| ((fnow > enow) == (fnow > -enow)))) {
			Fast_Two_Sum(enow, g0, Qnew, q);
			enow = e[++eindex];
		}
		else {
			Fast_Two_Sum(fnow, g0, Qnew, q);
			fnow = f[++findex];
		}
		Q = Qnew;
		for (count = 2; count < elen + flen; count++) {
			if ((eindex < elen) && ((findex >= flen)
				|| ((fnow > enow) == (fnow > -enow)))) {
				Fast_Two_Sum(enow, q, R, hh);
				enow = e[++eindex];
			}
			else {
				Fast_Two_Sum(fnow, q, R, hh);
				fnow = f[++findex];
			}
			Two_Sum(Q, R, Qnew, q);
			Q = Qnew;
			if (hh != 0) {
				h[hindex++] = hh;
			}
		}
		if (q != 0) {
			h[hindex++] = q;
		}
		if ((Q != 0.0) || (hindex == 0)) {
			h[hindex++] = Q;
		}
		return hindex;
	}

	/*****************************************************************************/
	/*                                                                           */
	/*  scale_expansion()   Multiply an expansion by a scalar.                   */
	/*                                                                           */
	/*  Sets h = be.  See either version of my paper for details.                */
	/*                                                                           */
	/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
	/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    */
	/*  properties as well.  (That is, if e has one of these properties, so      */
	/*  will h.)                                                                 */
	/*                                                                           */
	/*****************************************************************************/

	int scale_expansion(int elen, REAL *e, REAL b, REAL *h)
		/* e and h cannot be the same. */
	{
		INEXACT REAL Q;
		INEXACT REAL sum;
		INEXACT REAL product1;
		REAL product0;
		int eindex, hindex;
		REAL enow;
		INEXACT REAL bvirt;
		REAL avirt, bround, around;
		INEXACT REAL c;
		INEXACT REAL abig;
		REAL ahi, alo, bhi, blo;
		REAL err1, err2, err3;

		Split(b, bhi, blo);
		Two_Product_Presplit(e[0], b, bhi, blo, Q, h[0]);
		hindex = 1;
		for (eindex = 1; eindex < elen; eindex++) {
			enow = e[eindex];
			Two_Product_Presplit(enow, b, bhi, blo, product1, product0);
			Two_Sum(Q, product0, sum, h[hindex]);
			hindex++;
			Two_Sum(product1, sum, Q, h[hindex]);
			hindex++;
		}
		h[hindex] = Q;
		return elen + elen;
	}

	/*****************************************************************************/
	/*                                                                           */
	/*  scale_expansion_zeroelim()   Multiply an expansion by a scalar,          */
	/*                               eliminating zero components from the        */
	/*                               output expansion.                           */
	/*                                                                           */
	/*  Sets h = be.  See either version of my paper for details.                */
	/*                                                                           */
	/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
	/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    */
	/*  properties as well.  (That is, if e has one of these properties, so      */
	/*  will h.)                                                                 */
	/*                                                                           */
	/*****************************************************************************/

	int scale_expansion_zeroelim(int elen, REAL *e, REAL b, REAL *h)
		/* e and h cannot be the same. */
	{
		INEXACT REAL Q, sum;
		REAL hh;
		INEXACT REAL product1;
		REAL product0;
		int eindex, hindex;
		REAL enow;
		INEXACT REAL bvirt;
		REAL avirt, bround, around;
		INEXACT REAL c;
		INEXACT REAL abig;
		REAL ahi, alo, bhi, blo;
		REAL err1, err2, err3;

		Split(b, bhi, blo);
		Two_Product_Presplit(e[0], b, bhi, blo, Q, hh);
		hindex = 0;
		if (hh != 0) {
			h[hindex++] = hh;
		}
		for (eindex = 1; eindex < elen; eindex++) {
			enow = e[eindex];
			Two_Product_Presplit(enow, b, bhi, blo, product1, product0);
			Two_Sum(Q, product0, sum, hh);
			if (hh != 0) {
				h[hindex++] = hh;
			}
			Fast_Two_Sum(product1, sum, Q, hh);
			if (hh != 0) {
				h[hindex++] = hh;
			}
		}
		if ((Q != 0.0) || (hindex == 0)) {
			h[hindex++] = Q;
		}
		return hindex;
	}

	/*****************************************************************************/
	/*                                                                           */
	/*  compress()   Compress an expansion.                                      */
	/*                                                                           */
	/*  See the long version of my paper for details.                            */
	/*                                                                           */
	/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
	/*  with IEEE 754), then any nonoverlapping expansion is converted to a      */
	/*  nonadjacent expansion.                                                   */
	/*                                                                           */
	/*****************************************************************************/

	int compress(int elen, REAL *e, REAL *h)
		/* e and h may be the same. */
	{
		REAL Q, q;
		INEXACT REAL Qnew;
		int eindex, hindex;
		INEXACT REAL bvirt;
		REAL enow, hnow;
		int top, bottom;

		bottom = elen - 1;
		Q = e[bottom];
		for (eindex = elen - 2; eindex >= 0; eindex--) {
			enow = e[eindex];
			Fast_Two_Sum(Q, enow, Qnew, q);
			if (q != 0) {
				h[bottom--] = Qnew;
				Q = q;
			}
			else {
				Q = Qnew;
			}
		}
		top = 0;
		for (hindex = bottom + 1; hindex < elen; hindex++) {
			hnow = h[hindex];
			Fast_Two_Sum(hnow, Q, Qnew, q);
			if (q != 0) {
				h[top++] = q;
			}
			Q = Qnew;
		}
		h[top] = Q;
		return top + 1;
	}

	/*****************************************************************************/
	/*                                                                           */
	/*  estimate()   Produce a one-word estimate of an expansion's value.        */
	/*                                                                           */
	/*  See either version of my paper for details.                              */
	/*                                                                           */
	/*****************************************************************************/

	REAL estimate(int elen, REAL *e)
	{
		REAL Q;
		int eindex;

		Q = e[0];
		for (eindex = 1; eindex < elen; eindex++) {
			Q += e[eindex];
		}
		return Q;
	}

	/*****************************************************************************/
	/*                                                                           */
	/*  orient2dfast()   Approximate 2D orientation test.  Nonrobust.            */
	/*  orient2dexact()   Exact 2D orientation test.  Robust.                    */
	/*  orient2dslow()   Another exact 2D orientation test.  Robust.             */
	/*  orient2d()   Adaptive exact 2D orientation test.  Robust.                */
	/*                                                                           */
	/*               Return a positive value if the points pa, pb, and pc occur  */
	/*               in counterclockwise order; a negative value if they occur   */
	/*               in clockwise order; and zero if they are collinear.  The    */
	/*               result is also a rough approximation of twice the signed    */
	/*               area of the triangle defined by the three points.           */
	/*                                                                           */
	/*  Only the first and last routine should be used; the middle two are for   */
	/*  timings.                                                                 */
	/*                                                                           */
	/*  The last three use exact arithmetic to ensure a correct answer.  The     */
	/*  result returned is the determinant of a matrix.  In orient2d() only,     */
	/*  this determinant is computed adaptively, in the sense that exact         */
	/*  arithmetic is used only to the degree it is needed to ensure that the    */
	/*  returned value has the correct sign.  Hence, orient2d() is usually quite */
	/*  fast, but will run more slowly when the input points are collinear or    */
	/*  nearly so.                                                               */
	/*                                                                           */
	/*****************************************************************************/

	REAL orient2dfast(REAL *pa, REAL *pb, REAL *pc)
	{
		REAL acx, bcx, acy, bcy;

		acx = pa[0] - pc[0];
		bcx = pb[0] - pc[0];
		acy = pa[1] - pc[1];
		bcy = pb[1] - pc[1];
		return acx * bcy - acy * bcx;
	}

	REAL orient2dexact(REAL *pa, REAL *pb, REAL *pc)
	{
		INEXACT REAL axby1, axcy1, bxcy1, bxay1, cxay1, cxby1;
		REAL axby0, axcy0, bxcy0, bxay0, cxay0, cxby0;
		REAL aterms[4], bterms[4], cterms[4];
		INEXACT REAL aterms3, bterms3, cterms3;
		REAL v[8], w[12];
		int vlength, wlength;

		INEXACT REAL bvirt;
		REAL avirt, bround, around;
		INEXACT REAL c;
		INEXACT REAL abig;
		REAL ahi, alo, bhi, blo;
		REAL err1, err2, err3;
		INEXACT REAL _i, _j;
		REAL _0;

		Two_Product(pa[0], pb[1], axby1, axby0);
		Two_Product(pa[0], pc[1], axcy1, axcy0);
		Two_Two_Diff(axby1, axby0, axcy1, axcy0,
			aterms3, aterms[2], aterms[1], aterms[0]);
		aterms[3] = aterms3;

		Two_Product(pb[0], pc[1], bxcy1, bxcy0);
		Two_Product(pb[0], pa[1], bxay1, bxay0);
		Two_Two_Diff(bxcy1, bxcy0, bxay1, bxay0,
			bterms3, bterms[2], bterms[1], bterms[0]);
		bterms[3] = bterms3;

		Two_Product(pc[0], pa[1], cxay1, cxay0);
		Two_Product(pc[0], pb[1], cxby1, cxby0);
		Two_Two_Diff(cxay1, cxay0, cxby1, cxby0,
			cterms3, cterms[2], cterms[1], cterms[0]);
		cterms[3] = cterms3;

		vlength = fast_expansion_sum_zeroelim(4, aterms, 4, bterms, v);
		wlength = fast_expansion_sum_zeroelim(vlength, v, 4, cterms, w);

		return w[wlength - 1];
	}

	REAL orient2dslow(REAL *pa, REAL *pb, REAL *pc)
	{
		INEXACT REAL acx, acy, bcx, bcy;
		REAL acxtail, acytail;
		REAL bcxtail, bcytail;
		REAL negate, negatetail;
		REAL axby[8], bxay[8];
		INEXACT REAL axby7, bxay7;
		REAL deter[16];
		int deterlen;

		INEXACT REAL bvirt;
		REAL avirt, bround, around;
		INEXACT REAL c;
		INEXACT REAL abig;
		REAL a0hi, a0lo, a1hi, a1lo, bhi, blo;
		REAL err1, err2, err3;
		INEXACT REAL _i, _j, _k, _l, _m, _n;
		REAL _0, _1, _2;

		Two_Diff(pa[0], pc[0], acx, acxtail);
		Two_Diff(pa[1], pc[1], acy, acytail);
		Two_Diff(pb[0], pc[0], bcx, bcxtail);
		Two_Diff(pb[1], pc[1], bcy, bcytail);

		Two_Two_Product(acx, acxtail, bcy, bcytail,
			axby7, axby[6], axby[5], axby[4],
			axby[3], axby[2], axby[1], axby[0]);
		axby[7] = axby7;
		negate = -acy;
		negatetail = -acytail;
		Two_Two_Product(bcx, bcxtail, negate, negatetail,
			bxay7, bxay[6], bxay[5], bxay[4],
			bxay[3], bxay[2], bxay[1], bxay[0]);
		bxay[7] = bxay7;

		deterlen = fast_expansion_sum_zeroelim(8, axby, 8, bxay, deter);

		return deter[deterlen - 1];
	}

	REAL orient2dadapt(REAL *pa, REAL *pb, REAL *pc, REAL detsum)
	{
		INEXACT REAL acx, acy, bcx, bcy;
		REAL acxtail, acytail, bcxtail, bcytail;
		INEXACT REAL detleft, detright;
		REAL detlefttail, detrighttail;
		REAL det, errbound;
		REAL B[4], C1[8], C2[12], D[16];
		INEXACT REAL B3;
		int C1length, C2length, Dlength;
		REAL u[4];
		INEXACT REAL u3;
		INEXACT REAL s1, t1;
		REAL s0, t0;

		INEXACT REAL bvirt;
		REAL avirt, bround, around;
		INEXACT REAL c;
		INEXACT REAL abig;
		REAL ahi, alo, bhi, blo;
		REAL err1, err2, err3;
		INEXACT REAL _i, _j;
		REAL _0;

		acx = (REAL)(pa[0] - pc[0]);
		bcx = (REAL)(pb[0] - pc[0]);
		acy = (REAL)(pa[1] - pc[1]);
		bcy = (REAL)(pb[1] - pc[1]);

		Two_Product(acx, bcy, detleft, detlefttail);
		Two_Product(acy, bcx, detright, detrighttail);

		Two_Two_Diff(detleft, detlefttail, detright, detrighttail,
			B3, B[2], B[1], B[0]);
		B[3] = B3;

		det = estimate(4, B);
		errbound = ccwerrboundB * detsum;
		if ((det >= errbound) || (-det >= errbound)) {
			return det;
		}

		Two_Diff_Tail(pa[0], pc[0], acx, acxtail);
		Two_Diff_Tail(pb[0], pc[0], bcx, bcxtail);
		Two_Diff_Tail(pa[1], pc[1], acy, acytail);
		Two_Diff_Tail(pb[1], pc[1], bcy, bcytail);

		if ((acxtail == 0.0) && (acytail == 0.0)
			&& (bcxtail == 0.0) && (bcytail == 0.0)) {
			return det;
		}

		errbound = ccwerrboundC * detsum + resulterrbound * Absolute(det);
		det += (acx * bcytail + bcy * acxtail)
			- (acy * bcxtail + bcx * acytail);
		if ((det >= errbound) || (-det >= errbound)) {
			return det;
		}

		Two_Product(acxtail, bcy, s1, s0);
		Two_Product(acytail, bcx, t1, t0);
		Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
		u[3] = u3;
		C1length = fast_expansion_sum_zeroelim(4, B, 4, u, C1);

		Two_Product(acx, bcytail, s1, s0);
		Two_Product(acy, bcxtail, t1, t0);
		Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
		u[3] = u3;
		C2length = fast_expansion_sum_zeroelim(C1length, C1, 4, u, C2);

		Two_Product(acxtail, bcytail, s1, s0);
		Two_Product(acytail, bcxtail, t1, t0);
		Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
		u[3] = u3;
		Dlength = fast_expansion_sum_zeroelim(C2length, C2, 4, u, D);

		return(D[Dlength - 1]);
	}

	REAL orient2d(REAL *pa, REAL *pb, REAL *pc)
	{
		REAL detleft, detright, det;
		REAL detsum, errbound;

		detleft = (pa[0] - pc[0]) * (pb[1] - pc[1]);
		detright = (pa[1] - pc[1]) * (pb[0] - pc[0]);
		det = detleft - detright;

		if (detleft > 0.0) {
			if (detright <= 0.0) {
				return det;
			}
			else {
				detsum = detleft + detright;
			}
		}
		else if (detleft < 0.0) {
			if (detright >= 0.0) {
				return det;
			}
			else {
				detsum = -detleft - detright;
			}
		}
		else {
			return det;
		}

		errbound = ccwerrboundA * detsum;
		if ((det >= errbound) || (-det >= errbound)) {
			return det;
		}

		return orient2dadapt(pa, pb, pc, detsum);
	}

	/*****************************************************************************/
	/*                                                                           */
	/*  orient3dfast()   Approximate 3D orientation test.  Nonrobust.            */
	/*  orient3dexact()   Exact 3D orientation test.  Robust.                    */
	/*  orient3dslow()   Another exact 3D orientation test.  Robust.             */
	/*  orient3d()   Adaptive exact 3D orientation test.  Robust.                */
	/*                                                                           */
	/*               Return a positive value if the point pd lies below the      */
	/*               plane passing through pa, pb, and pc; "below" is defined so */
	/*               that pa, pb, and pc appear in counterclockwise order when   */
	/*               viewed from above the plane.  Returns a negative value if   */
	/*               pd lies above the plane.  Returns zero if the points are    */
	/*               coplanar.  The result is also a rough approximation of six  */
	/*               times the signed volume of the tetrahedron defined by the   */
	/*               four points.                                                */
	/*                                                                           */
	/*  Only the first and last routine should be used; the middle two are for   */
	/*  timings.                                                                 */
	/*                                                                           */
	/*  The last three use exact arithmetic to ensure a correct answer.  The     */
	/*  result returned is the determinant of a matrix.  In orient3d() only,     */
	/*  this determinant is computed adaptively, in the sense that exact         */
	/*  arithmetic is used only to the degree it is needed to ensure that the    */
	/*  returned value has the correct sign.  Hence, orient3d() is usually quite */
	/*  fast, but will run more slowly when the input points are coplanar or     */
	/*  nearly so.                                                               */
	/*                                                                           */
	/*****************************************************************************/

	REAL orient3dfast(REAL *pa, REAL *pb, REAL *pc, REAL *pd)
	{
		REAL adx, bdx, cdx;
		REAL ady, bdy, cdy;
		REAL adz, bdz, cdz;

		adx = pa[0] - pd[0];
		bdx = pb[0] - pd[0];
		cdx = pc[0] - pd[0];
		ady = pa[1] - pd[1];
		bdy = pb[1] - pd[1];
		cdy = pc[1] - pd[1];
		adz = pa[2] - pd[2];
		bdz = pb[2] - pd[2];
		cdz = pc[2] - pd[2];

		return adx * (bdy * cdz - bdz * cdy)
			+ bdx * (cdy * adz - cdz * ady)
			+ cdx * (ady * bdz - adz * bdy);
	}

	REAL orient3dexact(REAL *pa, REAL *pb, REAL *pc, REAL *pd)
	{
		INEXACT REAL axby1, bxcy1, cxdy1, dxay1, axcy1, bxdy1;
		INEXACT REAL bxay1, cxby1, dxcy1, axdy1, cxay1, dxby1;
		REAL axby0, bxcy0, cxdy0, dxay0, axcy0, bxdy0;
		REAL bxay0, cxby0, dxcy0, axdy0, cxay0, dxby0;
		REAL ab[4], bc[4], cd[4], da[4], ac[4], bd[4];
		REAL temp8[8];
		int templen;
		REAL abc[12], bcd[12], cda[12], dab[12];
		int abclen, bcdlen, cdalen, dablen;
		REAL adet[24], bdet[24], cdet[24], ddet[24];
		int alen, blen, clen, dlen;
		REAL abdet[48], cddet[48];
		int ablen, cdlen;
		REAL deter[96];
		int deterlen;
		int i;

		INEXACT REAL bvirt;
		REAL avirt, bround, around;
		INEXACT REAL c;
		INEXACT REAL abig;
		REAL ahi, alo, bhi, blo;
		REAL err1, err2, err3;
		INEXACT REAL _i, _j;
		REAL _0;

		Two_Product(pa[0], pb[1], axby1, axby0);
		Two_Product(pb[0], pa[1], bxay1, bxay0);
		Two_Two_Diff(axby1, axby0, bxay1, bxay0, ab[3], ab[2], ab[1], ab[0]);

		Two_Product(pb[0], pc[1], bxcy1, bxcy0);
		Two_Product(pc[0], pb[1], cxby1, cxby0);
		Two_Two_Diff(bxcy1, bxcy0, cxby1, cxby0, bc[3], bc[2], bc[1], bc[0]);

		Two_Product(pc[0], pd[1], cxdy1, cxdy0);
		Two_Product(pd[0], pc[1], dxcy1, dxcy0);
		Two_Two_Diff(cxdy1, cxdy0, dxcy1, dxcy0, cd[3], cd[2], cd[1], cd[0]);

		Two_Product(pd[0], pa[1], dxay1, dxay0);
		Two_Product(pa[0], pd[1], axdy1, axdy0);
		Two_Two_Diff(dxay1, dxay0, axdy1, axdy0, da[3], da[2], da[1], da[0]);

		Two_Product(pa[0], pc[1], axcy1, axcy0);
		Two_Product(pc[0], pa[1], cxay1, cxay0);
		Two_Two_Diff(axcy1, axcy0, cxay1, cxay0, ac[3], ac[2], ac[1], ac[0]);

		Two_Product(pb[0], pd[1], bxdy1, bxdy0);
		Two_Product(pd[0], pb[1], dxby1, dxby0);
		Two_Two_Diff(bxdy1, bxdy0, dxby1, dxby0, bd[3], bd[2], bd[1], bd[0]);

		templen = fast_expansion_sum_zeroelim(4, cd, 4, da, temp8);
		cdalen = fast_expansion_sum_zeroelim(templen, temp8, 4, ac, cda);
		templen = fast_expansion_sum_zeroelim(4, da, 4, ab, temp8);
		dablen = fast_expansion_sum_zeroelim(templen, temp8, 4, bd, dab);
		for (i = 0; i < 4; i++) {
			bd[i] = -bd[i];
			ac[i] = -ac[i];
		}
		templen = fast_expansion_sum_zeroelim(4, ab, 4, bc, temp8);
		abclen = fast_expansion_sum_zeroelim(templen, temp8, 4, ac, abc);
		templen = fast_expansion_sum_zeroelim(4, bc, 4, cd, temp8);
		bcdlen = fast_expansion_sum_zeroelim(templen, temp8, 4, bd, bcd);

		alen = scale_expansion_zeroelim(bcdlen, bcd, pa[2], adet);
		blen = scale_expansion_zeroelim(cdalen, cda, -pb[2], bdet);
		clen = scale_expansion_zeroelim(dablen, dab, pc[2], cdet);
		dlen = scale_expansion_zeroelim(abclen, abc, -pd[2], ddet);

		ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
		cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
		deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, deter);

		return deter[deterlen - 1];
	}

	REAL orient3dslow(REAL *pa, REAL *pb, REAL *pc, REAL *pd)
	{
		INEXACT REAL adx, ady, adz, bdx, bdy, bdz, cdx, cdy, cdz;
		REAL adxtail, adytail, adztail;
		REAL bdxtail, bdytail, bdztail;
		REAL cdxtail, cdytail, cdztail;
		REAL negate, negatetail;
		INEXACT REAL axby7, bxcy7, axcy7, bxay7, cxby7, cxay7;
		REAL axby[8], bxcy[8], axcy[8], bxay[8], cxby[8], cxay[8];
		REAL temp16[16], temp32[32], temp32t[32];
		int temp16len, temp32len, temp32tlen;
		REAL adet[64], bdet[64], cdet[64];
		int alen, blen, clen;
		REAL abdet[128];
		int ablen;
		REAL deter[192];
		int deterlen;

		INEXACT REAL bvirt;
		REAL avirt, bround, around;
		INEXACT REAL c;
		INEXACT REAL abig;
		REAL a0hi, a0lo, a1hi, a1lo, bhi, blo;
		REAL err1, err2, err3;
		INEXACT REAL _i, _j, _k, _l, _m, _n;
		REAL _0, _1, _2;

		Two_Diff(pa[0], pd[0], adx, adxtail);
		Two_Diff(pa[1], pd[1], ady, adytail);
		Two_Diff(pa[2], pd[2], adz, adztail);
		Two_Diff(pb[0], pd[0], bdx, bdxtail);
		Two_Diff(pb[1], pd[1], bdy, bdytail);
		Two_Diff(pb[2], pd[2], bdz, bdztail);
		Two_Diff(pc[0], pd[0], cdx, cdxtail);
		Two_Diff(pc[1], pd[1], cdy, cdytail);
		Two_Diff(pc[2], pd[2], cdz, cdztail);

		Two_Two_Product(adx, adxtail, bdy, bdytail,
			axby7, axby[6], axby[5], axby[4],
			axby[3], axby[2], axby[1], axby[0]);
		axby[7] = axby7;
		negate = -ady;
		negatetail = -adytail;
		Two_Two_Product(bdx, bdxtail, negate, negatetail,
			bxay7, bxay[6], bxay[5], bxay[4],
			bxay[3], bxay[2], bxay[1], bxay[0]);
		bxay[7] = bxay7;
		Two_Two_Product(bdx, bdxtail, cdy, cdytail,
			bxcy7, bxcy[6], bxcy[5], bxcy[4],
			bxcy[3], bxcy[2], bxcy[1], bxcy[0]);
		bxcy[7] = bxcy7;
		negate = -bdy;
		negatetail = -bdytail;
		Two_Two_Product(cdx, cdxtail, negate, negatetail,
			cxby7, cxby[6], cxby[5], cxby[4],
			cxby[3], cxby[2], cxby[1], cxby[0]);
		cxby[7] = cxby7;
		Two_Two_Product(cdx, cdxtail, ady, adytail,
			cxay7, cxay[6], cxay[5], cxay[4],
			cxay[3], cxay[2], cxay[1], cxay[0]);
		cxay[7] = cxay7;
		negate = -cdy;
		negatetail = -cdytail;
		Two_Two_Product(adx, adxtail, negate, negatetail,
			axcy7, axcy[6], axcy[5], axcy[4],
			axcy[3], axcy[2], axcy[1], axcy[0]);
		axcy[7] = axcy7;

		temp16len = fast_expansion_sum_zeroelim(8, bxcy, 8, cxby, temp16);
		temp32len = scale_expansion_zeroelim(temp16len, temp16, adz, temp32);
		temp32tlen = scale_expansion_zeroelim(temp16len, temp16, adztail, temp32t);
		alen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t,
			adet);

		temp16len = fast_expansion_sum_zeroelim(8, cxay, 8, axcy, temp16);
		temp32len = scale_expansion_zeroelim(temp16len, temp16, bdz, temp32);
		temp32tlen = scale_expansion_zeroelim(temp16len, temp16, bdztail, temp32t);
		blen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t,
			bdet);

		temp16len = fast_expansion_sum_zeroelim(8, axby, 8, bxay, temp16);
		temp32len = scale_expansion_zeroelim(temp16len, temp16, cdz, temp32);
		temp32tlen = scale_expansion_zeroelim(temp16len, temp16, cdztail, temp32t);
		clen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t,
			cdet);

		ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
		deterlen = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, deter);

		return deter[deterlen - 1];
	}

	REAL orient3dadapt(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL permanent)
	{
		INEXACT REAL adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz;
		REAL det, errbound;

		INEXACT REAL bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
		REAL bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
		REAL bc[4], ca[4], ab[4];
		INEXACT REAL bc3, ca3, ab3;
		REAL adet[8], bdet[8], cdet[8];
		int alen, blen, clen;
		REAL abdet[16];
		int ablen;
		REAL *finnow, *finother, *finswap;
		REAL fin1[192], fin2[192];
		int finlength;


		REAL adxtail, bdxtail, cdxtail;
		REAL adytail, bdytail, cdytail;
		REAL adztail, bdztail, cdztail;
		INEXACT REAL at_blarge, at_clarge;
		INEXACT REAL bt_clarge, bt_alarge;
		INEXACT REAL ct_alarge, ct_blarge;
		REAL at_b[4], at_c[4], bt_c[4], bt_a[4], ct_a[4], ct_b[4];
		int at_blen, at_clen, bt_clen, bt_alen, ct_alen, ct_blen;
		INEXACT REAL bdxt_cdy1, cdxt_bdy1, cdxt_ady1;
		INEXACT REAL adxt_cdy1, adxt_bdy1, bdxt_ady1;
		REAL bdxt_cdy0, cdxt_bdy0, cdxt_ady0;
		REAL adxt_cdy0, adxt_bdy0, bdxt_ady0;
		INEXACT REAL bdyt_cdx1, cdyt_bdx1, cdyt_adx1;
		INEXACT REAL adyt_cdx1, adyt_bdx1, bdyt_adx1;
		REAL bdyt_cdx0, cdyt_bdx0, cdyt_adx0;
		REAL adyt_cdx0, adyt_bdx0, bdyt_adx0;
		REAL bct[8], cat[8], abt[8];
		int bctlen, catlen, abtlen;
		INEXACT REAL bdxt_cdyt1, cdxt_bdyt1, cdxt_adyt1;
		INEXACT REAL adxt_cdyt1, adxt_bdyt1, bdxt_adyt1;
		REAL bdxt_cdyt0, cdxt_bdyt0, cdxt_adyt0;
		REAL adxt_cdyt0, adxt_bdyt0, bdxt_adyt0;
		REAL u[4], v[12], w[16];
		INEXACT REAL u3;
		int vlength, wlength;
		REAL negate;

		INEXACT REAL bvirt;
		REAL avirt, bround, around;
		INEXACT REAL c;
		INEXACT REAL abig;
		REAL ahi, alo, bhi, blo;
		REAL err1, err2, err3;
		INEXACT REAL _i, _j, _k;
		REAL _0;


		adx = (REAL)(pa[0] - pd[0]);
		bdx = (REAL)(pb[0] - pd[0]);
		cdx = (REAL)(pc[0] - pd[0]);
		ady = (REAL)(pa[1] - pd[1]);
		bdy = (REAL)(pb[1] - pd[1]);
		cdy = (REAL)(pc[1] - pd[1]);
		adz = (REAL)(pa[2] - pd[2]);
		bdz = (REAL)(pb[2] - pd[2]);
		cdz = (REAL)(pc[2] - pd[2]);

		Two_Product(bdx, cdy, bdxcdy1, bdxcdy0);
		Two_Product(cdx, bdy, cdxbdy1, cdxbdy0);
		Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc3, bc[2], bc[1], bc[0]);
		bc[3] = bc3;
		alen = scale_expansion_zeroelim(4, bc, adz, adet);

		Two_Product(cdx, ady, cdxady1, cdxady0);
		Two_Product(adx, cdy, adxcdy1, adxcdy0);
		Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca3, ca[2], ca[1], ca[0]);
		ca[3] = ca3;
		blen = scale_expansion_zeroelim(4, ca, bdz, bdet);

		Two_Product(adx, bdy, adxbdy1, adxbdy0);
		Two_Product(bdx, ady, bdxady1, bdxady0);
		Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab3, ab[2], ab[1], ab[0]);
		ab[3] = ab3;
		clen = scale_expansion_zeroelim(4, ab, cdz, cdet);

		ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
		finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, fin1);

		det = estimate(finlength, fin1);
		errbound = o3derrboundB * permanent;
		if ((det >= errbound) || (-det >= errbound)) {
			return det;
		}

		Two_Diff_Tail(pa[0], pd[0], adx, adxtail);
		Two_Diff_Tail(pb[0], pd[0], bdx, bdxtail);
		Two_Diff_Tail(pc[0], pd[0], cdx, cdxtail);
		Two_Diff_Tail(pa[1], pd[1], ady, adytail);
		Two_Diff_Tail(pb[1], pd[1], bdy, bdytail);
		Two_Diff_Tail(pc[1], pd[1], cdy, cdytail);
		Two_Diff_Tail(pa[2], pd[2], adz, adztail);
		Two_Diff_Tail(pb[2], pd[2], bdz, bdztail);
		Two_Diff_Tail(pc[2], pd[2], cdz, cdztail);

		if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0)
			&& (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0)
			&& (adztail == 0.0) && (bdztail == 0.0) && (cdztail == 0.0)) {
			return det;
		}

		errbound = o3derrboundC * permanent + resulterrbound * Absolute(det);
		det += (adz * ((bdx * cdytail + cdy * bdxtail)
			- (bdy * cdxtail + cdx * bdytail))
			+ adztail * (bdx * cdy - bdy * cdx))
			+ (bdz * ((cdx * adytail + ady * cdxtail)
			- (cdy * adxtail + adx * cdytail))
			+ bdztail * (cdx * ady - cdy * adx))
			+ (cdz * ((adx * bdytail + bdy * adxtail)
			- (ady * bdxtail + bdx * adytail))
			+ cdztail * (adx * bdy - ady * bdx));
		if ((det >= errbound) || (-det >= errbound)) {
			return det;
		}

		finnow = fin1;
		finother = fin2;

		if (adxtail == 0.0) {
			if (adytail == 0.0) {
				at_b[0] = 0.0;
				at_blen = 1;
				at_c[0] = 0.0;
				at_clen = 1;
			}
			else {
				negate = -adytail;
				Two_Product(negate, bdx, at_blarge, at_b[0]);
				at_b[1] = at_blarge;
				at_blen = 2;
				Two_Product(adytail, cdx, at_clarge, at_c[0]);
				at_c[1] = at_clarge;
				at_clen = 2;
			}
		}
		else {
			if (adytail == 0.0) {
				Two_Product(adxtail, bdy, at_blarge, at_b[0]);
				at_b[1] = at_blarge;
				at_blen = 2;
				negate = -adxtail;
				Two_Product(negate, cdy, at_clarge, at_c[0]);
				at_c[1] = at_clarge;
				at_clen = 2;
			}
			else {
				Two_Product(adxtail, bdy, adxt_bdy1, adxt_bdy0);
				Two_Product(adytail, bdx, adyt_bdx1, adyt_bdx0);
				Two_Two_Diff(adxt_bdy1, adxt_bdy0, adyt_bdx1, adyt_bdx0,
					at_blarge, at_b[2], at_b[1], at_b[0]);
				at_b[3] = at_blarge;
				at_blen = 4;
				Two_Product(adytail, cdx, adyt_cdx1, adyt_cdx0);
				Two_Product(adxtail, cdy, adxt_cdy1, adxt_cdy0);
				Two_Two_Diff(adyt_cdx1, adyt_cdx0, adxt_cdy1, adxt_cdy0,
					at_clarge, at_c[2], at_c[1], at_c[0]);
				at_c[3] = at_clarge;
				at_clen = 4;
			}
		}
		if (bdxtail == 0.0) {
			if (bdytail == 0.0) {
				bt_c[0] = 0.0;
				bt_clen = 1;
				bt_a[0] = 0.0;
				bt_alen = 1;
			}
			else {
				negate = -bdytail;
				Two_Product(negate, cdx, bt_clarge, bt_c[0]);
				bt_c[1] = bt_clarge;
				bt_clen = 2;
				Two_Product(bdytail, adx, bt_alarge, bt_a[0]);
				bt_a[1] = bt_alarge;
				bt_alen = 2;
			}
		}
		else {
			if (bdytail == 0.0) {
				Two_Product(bdxtail, cdy, bt_clarge, bt_c[0]);
				bt_c[1] = bt_clarge;
				bt_clen = 2;
				negate = -bdxtail;
				Two_Product(negate, ady, bt_alarge, bt_a[0]);
				bt_a[1] = bt_alarge;
				bt_alen = 2;
			}
			else {
				Two_Product(bdxtail, cdy, bdxt_cdy1, bdxt_cdy0);
				Two_Product(bdytail, cdx, bdyt_cdx1, bdyt_cdx0);
				Two_Two_Diff(bdxt_cdy1, bdxt_cdy0, bdyt_cdx1, bdyt_cdx0,
					bt_clarge, bt_c[2], bt_c[1], bt_c[0]);
				bt_c[3] = bt_clarge;
				bt_clen = 4;
				Two_Product(bdytail, adx, bdyt_adx1, bdyt_adx0);
				Two_Product(bdxtail, ady, bdxt_ady1, bdxt_ady0);
				Two_Two_Diff(bdyt_adx1, bdyt_adx0, bdxt_ady1, bdxt_ady0,
					bt_alarge, bt_a[2], bt_a[1], bt_a[0]);
				bt_a[3] = bt_alarge;
				bt_alen = 4;
			}
		}
		if (cdxtail == 0.0) {
			if (cdytail == 0.0) {
				ct_a[0] = 0.0;
				ct_alen = 1;
				ct_b[0] = 0.0;
				ct_blen = 1;
			}
			else {
				negate = -cdytail;
				Two_Product(negate, adx, ct_alarge, ct_a[0]);
				ct_a[1] = ct_alarge;
				ct_alen = 2;
				Two_Product(cdytail, bdx, ct_blarge, ct_b[0]);
				ct_b[1] = ct_blarge;
				ct_blen = 2;
			}
		}
		else {
			if (cdytail == 0.0) {
				Two_Product(cdxtail, ady, ct_alarge, ct_a[0]);
				ct_a[1] = ct_alarge;
				ct_alen = 2;
				negate = -cdxtail;
				Two_Product(negate, bdy, ct_blarge, ct_b[0]);
				ct_b[1] = ct_blarge;
				ct_blen = 2;
			}
			else {
				Two_Product(cdxtail, ady, cdxt_ady1, cdxt_ady0);
				Two_Product(cdytail, adx, cdyt_adx1, cdyt_adx0);
				Two_Two_Diff(cdxt_ady1, cdxt_ady0, cdyt_adx1, cdyt_adx0,
					ct_alarge, ct_a[2], ct_a[1], ct_a[0]);
				ct_a[3] = ct_alarge;
				ct_alen = 4;
				Two_Product(cdytail, bdx, cdyt_bdx1, cdyt_bdx0);
				Two_Product(cdxtail, bdy, cdxt_bdy1, cdxt_bdy0);
				Two_Two_Diff(cdyt_bdx1, cdyt_bdx0, cdxt_bdy1, cdxt_bdy0,
					ct_blarge, ct_b[2], ct_b[1], ct_b[0]);
				ct_b[3] = ct_blarge;
				ct_blen = 4;
			}
		}

		bctlen = fast_expansion_sum_zeroelim(bt_clen, bt_c, ct_blen, ct_b, bct);
		wlength = scale_expansion_zeroelim(bctlen, bct, adz, w);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
			finother);
		finswap = finnow; finnow = finother; finother = finswap;

		catlen = fast_expansion_sum_zeroelim(ct_alen, ct_a, at_clen, at_c, cat);
		wlength = scale_expansion_zeroelim(catlen, cat, bdz, w);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
			finother);
		finswap = finnow; finnow = finother; finother = finswap;

		abtlen = fast_expansion_sum_zeroelim(at_blen, at_b, bt_alen, bt_a, abt);
		wlength = scale_expansion_zeroelim(abtlen, abt, cdz, w);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
			finother);
		finswap = finnow; finnow = finother; finother = finswap;

		if (adztail != 0.0) {
			vlength = scale_expansion_zeroelim(4, bc, adztail, v);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
				finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
		if (bdztail != 0.0) {
			vlength = scale_expansion_zeroelim(4, ca, bdztail, v);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
				finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
		if (cdztail != 0.0) {
			vlength = scale_expansion_zeroelim(4, ab, cdztail, v);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
				finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}

		if (adxtail != 0.0) {
			if (bdytail != 0.0) {
				Two_Product(adxtail, bdytail, adxt_bdyt1, adxt_bdyt0);
				Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdz, u3, u[2], u[1], u[0]);
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
					finother);
				finswap = finnow; finnow = finother; finother = finswap;
				if (cdztail != 0.0) {
					Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdztail, u3, u[2], u[1], u[0]);
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
						finother);
					finswap = finnow; finnow = finother; finother = finswap;
				}
			}
			if (cdytail != 0.0) {
				negate = -adxtail;
				Two_Product(negate, cdytail, adxt_cdyt1, adxt_cdyt0);
				Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdz, u3, u[2], u[1], u[0]);
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
					finother);
				finswap = finnow; finnow = finother; finother = finswap;
				if (bdztail != 0.0) {
					Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdztail, u3, u[2], u[1], u[0]);
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
						finother);
					finswap = finnow; finnow = finother; finother = finswap;
				}
			}
		}
		if (bdxtail != 0.0) {
			if (cdytail != 0.0) {
				Two_Product(bdxtail, cdytail, bdxt_cdyt1, bdxt_cdyt0);
				Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adz, u3, u[2], u[1], u[0]);
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
					finother);
				finswap = finnow; finnow = finother; finother = finswap;
				if (adztail != 0.0) {
					Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adztail, u3, u[2], u[1], u[0]);
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
						finother);
					finswap = finnow; finnow = finother; finother = finswap;
				}
			}
			if (adytail != 0.0) {
				negate = -bdxtail;
				Two_Product(negate, adytail, bdxt_adyt1, bdxt_adyt0);
				Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdz, u3, u[2], u[1], u[0]);
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
					finother);
				finswap = finnow; finnow = finother; finother = finswap;
				if (cdztail != 0.0) {
					Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdztail, u3, u[2], u[1], u[0]);
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
						finother);
					finswap = finnow; finnow = finother; finother = finswap;
				}
			}
		}
		if (cdxtail != 0.0) {
			if (adytail != 0.0) {
				Two_Product(cdxtail, adytail, cdxt_adyt1, cdxt_adyt0);
				Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdz, u3, u[2], u[1], u[0]);
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
					finother);
				finswap = finnow; finnow = finother; finother = finswap;
				if (bdztail != 0.0) {
					Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdztail, u3, u[2], u[1], u[0]);
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
						finother);
					finswap = finnow; finnow = finother; finother = finswap;
				}
			}
			if (bdytail != 0.0) {
				negate = -cdxtail;
				Two_Product(negate, bdytail, cdxt_bdyt1, cdxt_bdyt0);
				Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adz, u3, u[2], u[1], u[0]);
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
					finother);
				finswap = finnow; finnow = finother; finother = finswap;
				if (adztail != 0.0) {
					Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adztail, u3, u[2], u[1], u[0]);
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
						finother);
					finswap = finnow; finnow = finother; finother = finswap;
				}
			}
		}

		if (adztail != 0.0) {
			wlength = scale_expansion_zeroelim(bctlen, bct, adztail, w);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
				finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
		if (bdztail != 0.0) {
			wlength = scale_expansion_zeroelim(catlen, cat, bdztail, w);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
				finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
		if (cdztail != 0.0) {
			wlength = scale_expansion_zeroelim(abtlen, abt, cdztail, w);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
				finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}

		return finnow[finlength - 1];
	}

#ifdef USE_CGAL_PREDICATES

	REAL orient3d(REAL *pa, REAL *pb, REAL *pc, REAL *pd)
	{
		return (REAL)
			-cgal_pred_obj.orientation_3_object()
			(Point(pa[0], pa[1], pa[2]),
			Point(pb[0], pb[1], pb[2]),
			Point(pc[0], pc[1], pc[2]),
			Point(pd[0], pd[1], pd[2]));
	}

#else

	REAL orient3d(REAL *pa, REAL *pb, REAL *pc, REAL *pd)
	{
		REAL adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz;
		REAL bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady;
		REAL det;


		adx = pa[0] - pd[0];
		ady = pa[1] - pd[1];
		adz = pa[2] - pd[2];
		bdx = pb[0] - pd[0];
		bdy = pb[1] - pd[1];
		bdz = pb[2] - pd[2];
		cdx = pc[0] - pd[0];
		cdy = pc[1] - pd[1];
		cdz = pc[2] - pd[2];

		bdxcdy = bdx * cdy;
		cdxbdy = cdx * bdy;

		cdxady = cdx * ady;
		adxcdy = adx * cdy;

		adxbdy = adx * bdy;
		bdxady = bdx * ady;

		det = adz * (bdxcdy - cdxbdy)
			+ bdz * (cdxady - adxcdy)
			+ cdz * (adxbdy - bdxady);

		if (_use_inexact_arith) {
			return det;
		}

		if (_use_static_filter) {
			//if (fabs(det) > o3dstaticfilter) return det;
			if (det > o3dstaticfilter) return det;
			if (det < -o3dstaticfilter) return det;
		}


		REAL permanent, errbound;

		permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * Absolute(adz)
			+ (Absolute(cdxady) + Absolute(adxcdy)) * Absolute(bdz)
			+ (Absolute(adxbdy) + Absolute(bdxady)) * Absolute(cdz);
		errbound = o3derrboundA * permanent;
		if ((det > errbound) || (-det > errbound)) {
			return det;
		}

		return orient3dadapt(pa, pb, pc, pd, permanent);
	}

#endif // #ifdef USE_CGAL_PREDICATES

	/*****************************************************************************/
	/*                                                                           */
	/*  incirclefast()   Approximate 2D incircle test.  Nonrobust.               */
	/*  incircleexact()   Exact 2D incircle test.  Robust.                       */
	/*  incircleslow()   Another exact 2D incircle test.  Robust.                */
	/*  incircle()   Adaptive exact 2D incircle test.  Robust.                   */
	/*                                                                           */
	/*               Return a positive value if the point pd lies inside the     */
	/*               circle passing through pa, pb, and pc; a negative value if  */
	/*               it lies outside; and zero if the four points are cocircular.*/
	/*               The points pa, pb, and pc must be in counterclockwise       */
	/*               order, or the sign of the result will be reversed.          */
	/*                                                                           */
	/*  Only the first and last routine should be used; the middle two are for   */
	/*  timings.                                                                 */
	/*                                                                           */
	/*  The last three use exact arithmetic to ensure a correct answer.  The     */
	/*  result returned is the determinant of a matrix.  In incircle() only,     */
	/*  this determinant is computed adaptively, in the sense that exact         */
	/*  arithmetic is used only to the degree it is needed to ensure that the    */
	/*  returned value has the correct sign.  Hence, incircle() is usually quite */
	/*  fast, but will run more slowly when the input points are cocircular or   */
	/*  nearly so.                                                               */
	/*                                                                           */
	/*****************************************************************************/

	REAL incirclefast(REAL *pa, REAL *pb, REAL *pc, REAL *pd)
	{
		REAL adx, ady, bdx, bdy, cdx, cdy;
		REAL abdet, bcdet, cadet;
		REAL alift, blift, clift;

		adx = pa[0] - pd[0];
		ady = pa[1] - pd[1];
		bdx = pb[0] - pd[0];
		bdy = pb[1] - pd[1];
		cdx = pc[0] - pd[0];
		cdy = pc[1] - pd[1];

		abdet = adx * bdy - bdx * ady;
		bcdet = bdx * cdy - cdx * bdy;
		cadet = cdx * ady - adx * cdy;
		alift = adx * adx + ady * ady;
		blift = bdx * bdx + bdy * bdy;
		clift = cdx * cdx + cdy * cdy;

		return alift * bcdet + blift * cadet + clift * abdet;
	}

	REAL incircleexact(REAL *pa, REAL *pb, REAL *pc, REAL *pd)
	{
		INEXACT REAL axby1, bxcy1, cxdy1, dxay1, axcy1, bxdy1;
		INEXACT REAL bxay1, cxby1, dxcy1, axdy1, cxay1, dxby1;
		REAL axby0, bxcy0, cxdy0, dxay0, axcy0, bxdy0;
		REAL bxay0, cxby0, dxcy0, axdy0, cxay0, dxby0;
		REAL ab[4], bc[4], cd[4], da[4], ac[4], bd[4];
		REAL temp8[8];
		int templen;
		REAL abc[12], bcd[12], cda[12], dab[12];
		int abclen, bcdlen, cdalen, dablen;
		REAL det24x[24], det24y[24], det48x[48], det48y[48];
		int xlen, ylen;
		REAL adet[96], bdet[96], cdet[96], ddet[96];
		int alen, blen, clen, dlen;
		REAL abdet[192], cddet[192];
		int ablen, cdlen;
		REAL deter[384];
		int deterlen;
		int i;

		INEXACT REAL bvirt;
		REAL avirt, bround, around;
		INEXACT REAL c;
		INEXACT REAL abig;
		REAL ahi, alo, bhi, blo;
		REAL err1, err2, err3;
		INEXACT REAL _i, _j;
		REAL _0;

		Two_Product(pa[0], pb[1], axby1, axby0);
		Two_Product(pb[0], pa[1], bxay1, bxay0);
		Two_Two_Diff(axby1, axby0, bxay1, bxay0, ab[3], ab[2], ab[1], ab[0]);

		Two_Product(pb[0], pc[1], bxcy1, bxcy0);
		Two_Product(pc[0], pb[1], cxby1, cxby0);
		Two_Two_Diff(bxcy1, bxcy0, cxby1, cxby0, bc[3], bc[2], bc[1], bc[0]);

		Two_Product(pc[0], pd[1], cxdy1, cxdy0);
		Two_Product(pd[0], pc[1], dxcy1, dxcy0);
		Two_Two_Diff(cxdy1, cxdy0, dxcy1, dxcy0, cd[3], cd[2], cd[1], cd[0]);

		Two_Product(pd[0], pa[1], dxay1, dxay0);
		Two_Product(pa[0], pd[1], axdy1, axdy0);
		Two_Two_Diff(dxay1, dxay0, axdy1, axdy0, da[3], da[2], da[1], da[0]);

		Two_Product(pa[0], pc[1], axcy1, axcy0);
		Two_Product(pc[0], pa[1], cxay1, cxay0);
		Two_Two_Diff(axcy1, axcy0, cxay1, cxay0, ac[3], ac[2], ac[1], ac[0]);

		Two_Product(pb[0], pd[1], bxdy1, bxdy0);
		Two_Product(pd[0], pb[1], dxby1, dxby0);
		Two_Two_Diff(bxdy1, bxdy0, dxby1, dxby0, bd[3], bd[2], bd[1], bd[0]);

		templen = fast_expansion_sum_zeroelim(4, cd, 4, da, temp8);
		cdalen = fast_expansion_sum_zeroelim(templen, temp8, 4, ac, cda);
		templen = fast_expansion_sum_zeroelim(4, da, 4, ab, temp8);
		dablen = fast_expansion_sum_zeroelim(templen, temp8, 4, bd, dab);
		for (i = 0; i < 4; i++) {
			bd[i] = -bd[i];
			ac[i] = -ac[i];
		}
		templen = fast_expansion_sum_zeroelim(4, ab, 4, bc, temp8);
		abclen = fast_expansion_sum_zeroelim(templen, temp8, 4, ac, abc);
		templen = fast_expansion_sum_zeroelim(4, bc, 4, cd, temp8);
		bcdlen = fast_expansion_sum_zeroelim(templen, temp8, 4, bd, bcd);

		xlen = scale_expansion_zeroelim(bcdlen, bcd, pa[0], det24x);
		xlen = scale_expansion_zeroelim(xlen, det24x, pa[0], det48x);
		ylen = scale_expansion_zeroelim(bcdlen, bcd, pa[1], det24y);
		ylen = scale_expansion_zeroelim(ylen, det24y, pa[1], det48y);
		alen = fast_expansion_sum_zeroelim(xlen, det48x, ylen, det48y, adet);

		xlen = scale_expansion_zeroelim(cdalen, cda, pb[0], det24x);
		xlen = scale_expansion_zeroelim(xlen, det24x, -pb[0], det48x);
		ylen = scale_expansion_zeroelim(cdalen, cda, pb[1], det24y);
		ylen = scale_expansion_zeroelim(ylen, det24y, -pb[1], det48y);
		blen = fast_expansion_sum_zeroelim(xlen, det48x, ylen, det48y, bdet);

		xlen = scale_expansion_zeroelim(dablen, dab, pc[0], det24x);
		xlen = scale_expansion_zeroelim(xlen, det24x, pc[0], det48x);
		ylen = scale_expansion_zeroelim(dablen, dab, pc[1], det24y);
		ylen = scale_expansion_zeroelim(ylen, det24y, pc[1], det48y);
		clen = fast_expansion_sum_zeroelim(xlen, det48x, ylen, det48y, cdet);

		xlen = scale_expansion_zeroelim(abclen, abc, pd[0], det24x);
		xlen = scale_expansion_zeroelim(xlen, det24x, -pd[0], det48x);
		ylen = scale_expansion_zeroelim(abclen, abc, pd[1], det24y);
		ylen = scale_expansion_zeroelim(ylen, det24y, -pd[1], det48y);
		dlen = fast_expansion_sum_zeroelim(xlen, det48x, ylen, det48y, ddet);

		ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
		cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
		deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, deter);

		return deter[deterlen - 1];
	}

	REAL incircleslow(REAL *pa, REAL *pb, REAL *pc, REAL *pd)
	{
		INEXACT REAL adx, bdx, cdx, ady, bdy, cdy;
		REAL adxtail, bdxtail, cdxtail;
		REAL adytail, bdytail, cdytail;
		REAL negate, negatetail;
		INEXACT REAL axby7, bxcy7, axcy7, bxay7, cxby7, cxay7;
		REAL axby[8], bxcy[8], axcy[8], bxay[8], cxby[8], cxay[8];
		REAL temp16[16];
		int temp16len;
		REAL detx[32], detxx[64], detxt[32], detxxt[64], detxtxt[64];
		int xlen, xxlen, xtlen, xxtlen, xtxtlen;
		REAL x1[128], x2[192];
		int x1len, x2len;
		REAL dety[32], detyy[64], detyt[32], detyyt[64], detytyt[64];
		int ylen, yylen, ytlen, yytlen, ytytlen;
		REAL y1[128], y2[192];
		int y1len, y2len;
		REAL adet[384], bdet[384], cdet[384], abdet[768], deter[1152];
		int alen, blen, clen, ablen, deterlen;
		int i;

		INEXACT REAL bvirt;
		REAL avirt, bround, around;
		INEXACT REAL c;
		INEXACT REAL abig;
		REAL a0hi, a0lo, a1hi, a1lo, bhi, blo;
		REAL err1, err2, err3;
		INEXACT REAL _i, _j, _k, _l, _m, _n;
		REAL _0, _1, _2;

		Two_Diff(pa[0], pd[0], adx, adxtail);
		Two_Diff(pa[1], pd[1], ady, adytail);
		Two_Diff(pb[0], pd[0], bdx, bdxtail);
		Two_Diff(pb[1], pd[1], bdy, bdytail);
		Two_Diff(pc[0], pd[0], cdx, cdxtail);
		Two_Diff(pc[1], pd[1], cdy, cdytail);

		Two_Two_Product(adx, adxtail, bdy, bdytail,
			axby7, axby[6], axby[5], axby[4],
			axby[3], axby[2], axby[1], axby[0]);
		axby[7] = axby7;
		negate = -ady;
		negatetail = -adytail;
		Two_Two_Product(bdx, bdxtail, negate, negatetail,
			bxay7, bxay[6], bxay[5], bxay[4],
			bxay[3], bxay[2], bxay[1], bxay[0]);
		bxay[7] = bxay7;
		Two_Two_Product(bdx, bdxtail, cdy, cdytail,
			bxcy7, bxcy[6], bxcy[5], bxcy[4],
			bxcy[3], bxcy[2], bxcy[1], bxcy[0]);
		bxcy[7] = bxcy7;
		negate = -bdy;
		negatetail = -bdytail;
		Two_Two_Product(cdx, cdxtail, negate, negatetail,
			cxby7, cxby[6], cxby[5], cxby[4],
			cxby[3], cxby[2], cxby[1], cxby[0]);
		cxby[7] = cxby7;
		Two_Two_Product(cdx, cdxtail, ady, adytail,
			cxay7, cxay[6], cxay[5], cxay[4],
			cxay[3], cxay[2], cxay[1], cxay[0]);
		cxay[7] = cxay7;
		negate = -cdy;
		negatetail = -cdytail;
		Two_Two_Product(adx, adxtail, negate, negatetail,
			axcy7, axcy[6], axcy[5], axcy[4],
			axcy[3], axcy[2], axcy[1], axcy[0]);
		axcy[7] = axcy7;


		temp16len = fast_expansion_sum_zeroelim(8, bxcy, 8, cxby, temp16);

		xlen = scale_expansion_zeroelim(temp16len, temp16, adx, detx);
		xxlen = scale_expansion_zeroelim(xlen, detx, adx, detxx);
		xtlen = scale_expansion_zeroelim(temp16len, temp16, adxtail, detxt);
		xxtlen = scale_expansion_zeroelim(xtlen, detxt, adx, detxxt);
		for (i = 0; i < xxtlen; i++) {
			detxxt[i] *= 2.0;
		}
		xtxtlen = scale_expansion_zeroelim(xtlen, detxt, adxtail, detxtxt);
		x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1);
		x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2);

		ylen = scale_expansion_zeroelim(temp16len, temp16, ady, dety);
		yylen = scale_expansion_zeroelim(ylen, dety, ady, detyy);
		ytlen = scale_expansion_zeroelim(temp16len, temp16, adytail, detyt);
		yytlen = scale_expansion_zeroelim(ytlen, detyt, ady, detyyt);
		for (i = 0; i < yytlen; i++) {
			detyyt[i] *= 2.0;
		}
		ytytlen = scale_expansion_zeroelim(ytlen, detyt, adytail, detytyt);
		y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1);
		y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2);

		alen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, adet);


		temp16len = fast_expansion_sum_zeroelim(8, cxay, 8, axcy, temp16);

		xlen = scale_expansion_zeroelim(temp16len, temp16, bdx, detx);
		xxlen = scale_expansion_zeroelim(xlen, detx, bdx, detxx);
		xtlen = scale_expansion_zeroelim(temp16len, temp16, bdxtail, detxt);
		xxtlen = scale_expansion_zeroelim(xtlen, detxt, bdx, detxxt);
		for (i = 0; i < xxtlen; i++) {
			detxxt[i] *= 2.0;
		}
		xtxtlen = scale_expansion_zeroelim(xtlen, detxt, bdxtail, detxtxt);
		x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1);
		x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2);

		ylen = scale_expansion_zeroelim(temp16len, temp16, bdy, dety);
		yylen = scale_expansion_zeroelim(ylen, dety, bdy, detyy);
		ytlen = scale_expansion_zeroelim(temp16len, temp16, bdytail, detyt);
		yytlen = scale_expansion_zeroelim(ytlen, detyt, bdy, detyyt);
		for (i = 0; i < yytlen; i++) {
			detyyt[i] *= 2.0;
		}
		ytytlen = scale_expansion_zeroelim(ytlen, detyt, bdytail, detytyt);
		y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1);
		y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2);

		blen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, bdet);


		temp16len = fast_expansion_sum_zeroelim(8, axby, 8, bxay, temp16);

		xlen = scale_expansion_zeroelim(temp16len, temp16, cdx, detx);
		xxlen = scale_expansion_zeroelim(xlen, detx, cdx, detxx);
		xtlen = scale_expansion_zeroelim(temp16len, temp16, cdxtail, detxt);
		xxtlen = scale_expansion_zeroelim(xtlen, detxt, cdx, detxxt);
		for (i = 0; i < xxtlen; i++) {
			detxxt[i] *= 2.0;
		}
		xtxtlen = scale_expansion_zeroelim(xtlen, detxt, cdxtail, detxtxt);
		x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1);
		x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2);

		ylen = scale_expansion_zeroelim(temp16len, temp16, cdy, dety);
		yylen = scale_expansion_zeroelim(ylen, dety, cdy, detyy);
		ytlen = scale_expansion_zeroelim(temp16len, temp16, cdytail, detyt);
		yytlen = scale_expansion_zeroelim(ytlen, detyt, cdy, detyyt);
		for (i = 0; i < yytlen; i++) {
			detyyt[i] *= 2.0;
		}
		ytytlen = scale_expansion_zeroelim(ytlen, detyt, cdytail, detytyt);
		y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1);
		y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2);

		clen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, cdet);

		ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
		deterlen = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, deter);

		return deter[deterlen - 1];
	}

	REAL incircleadapt(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL permanent)
	{
		INEXACT REAL adx, bdx, cdx, ady, bdy, cdy;
		REAL det, errbound;

		INEXACT REAL bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
		REAL bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
		REAL bc[4], ca[4], ab[4];
		INEXACT REAL bc3, ca3, ab3;
		REAL axbc[8], axxbc[16], aybc[8], ayybc[16], adet[32];
		int axbclen, axxbclen, aybclen, ayybclen, alen;
		REAL bxca[8], bxxca[16], byca[8], byyca[16], bdet[32];
		int bxcalen, bxxcalen, bycalen, byycalen, blen;
		REAL cxab[8], cxxab[16], cyab[8], cyyab[16], cdet[32];
		int cxablen, cxxablen, cyablen, cyyablen, clen;
		REAL abdet[64];
		int ablen;
		REAL fin1[1152], fin2[1152];
		REAL *finnow, *finother, *finswap;
		int finlength;

		REAL adxtail, bdxtail, cdxtail, adytail, bdytail, cdytail;
		INEXACT REAL adxadx1, adyady1, bdxbdx1, bdybdy1, cdxcdx1, cdycdy1;
		REAL adxadx0, adyady0, bdxbdx0, bdybdy0, cdxcdx0, cdycdy0;
		REAL aa[4], bb[4], cc[4];
		INEXACT REAL aa3, bb3, cc3;
		INEXACT REAL ti1, tj1;
		REAL ti0, tj0;
		REAL u[4], v[4];
		INEXACT REAL u3, v3;
		REAL temp8[8], temp16a[16], temp16b[16], temp16c[16];
		REAL temp32a[32], temp32b[32], temp48[48], temp64[64];
		int temp8len, temp16alen, temp16blen, temp16clen;
		int temp32alen, temp32blen, temp48len, temp64len;
		REAL axtbb[8], axtcc[8], aytbb[8], aytcc[8];
		int axtbblen, axtcclen, aytbblen, aytcclen;
		REAL bxtaa[8], bxtcc[8], bytaa[8], bytcc[8];
		int bxtaalen, bxtcclen, bytaalen, bytcclen;
		REAL cxtaa[8], cxtbb[8], cytaa[8], cytbb[8];
		int cxtaalen, cxtbblen, cytaalen, cytbblen;
		REAL axtbc[8], aytbc[8], bxtca[8], bytca[8], cxtab[8], cytab[8];
		int axtbclen, aytbclen, bxtcalen, bytcalen, cxtablen, cytablen;
		REAL axtbct[16], aytbct[16], bxtcat[16], bytcat[16], cxtabt[16], cytabt[16];
		int axtbctlen, aytbctlen, bxtcatlen, bytcatlen, cxtabtlen, cytabtlen;
		REAL axtbctt[8], aytbctt[8], bxtcatt[8];
		REAL bytcatt[8], cxtabtt[8], cytabtt[8];
		int axtbcttlen, aytbcttlen, bxtcattlen, bytcattlen, cxtabttlen, cytabttlen;
		REAL abt[8], bct[8], cat[8];
		int abtlen, bctlen, catlen;
		REAL abtt[4], bctt[4], catt[4];
		int abttlen, bcttlen, cattlen;
		INEXACT REAL abtt3, bctt3, catt3;
		REAL negate;

		INEXACT REAL bvirt;
		REAL avirt, bround, around;
		INEXACT REAL c;
		INEXACT REAL abig;
		REAL ahi, alo, bhi, blo;
		REAL err1, err2, err3;
		INEXACT REAL _i, _j;
		REAL _0;

		// Avoid compiler warnings. H. Si, 2012-02-16.
		axtbclen = aytbclen = bxtcalen = bytcalen = cxtablen = cytablen = 0;

		adx = (REAL)(pa[0] - pd[0]);
		bdx = (REAL)(pb[0] - pd[0]);
		cdx = (REAL)(pc[0] - pd[0]);
		ady = (REAL)(pa[1] - pd[1]);
		bdy = (REAL)(pb[1] - pd[1]);
		cdy = (REAL)(pc[1] - pd[1]);

		Two_Product(bdx, cdy, bdxcdy1, bdxcdy0);
		Two_Product(cdx, bdy, cdxbdy1, cdxbdy0);
		Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc3, bc[2], bc[1], bc[0]);
		bc[3] = bc3;
		axbclen = scale_expansion_zeroelim(4, bc, adx, axbc);
		axxbclen = scale_expansion_zeroelim(axbclen, axbc, adx, axxbc);
		aybclen = scale_expansion_zeroelim(4, bc, ady, aybc);
		ayybclen = scale_expansion_zeroelim(aybclen, aybc, ady, ayybc);
		alen = fast_expansion_sum_zeroelim(axxbclen, axxbc, ayybclen, ayybc, adet);

		Two_Product(cdx, ady, cdxady1, cdxady0);
		Two_Product(adx, cdy, adxcdy1, adxcdy0);
		Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca3, ca[2], ca[1], ca[0]);
		ca[3] = ca3;
		bxcalen = scale_expansion_zeroelim(4, ca, bdx, bxca);
		bxxcalen = scale_expansion_zeroelim(bxcalen, bxca, bdx, bxxca);
		bycalen = scale_expansion_zeroelim(4, ca, bdy, byca);
		byycalen = scale_expansion_zeroelim(bycalen, byca, bdy, byyca);
		blen = fast_expansion_sum_zeroelim(bxxcalen, bxxca, byycalen, byyca, bdet);

		Two_Product(adx, bdy, adxbdy1, adxbdy0);
		Two_Product(bdx, ady, bdxady1, bdxady0);
		Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab3, ab[2], ab[1], ab[0]);
		ab[3] = ab3;
		cxablen = scale_expansion_zeroelim(4, ab, cdx, cxab);
		cxxablen = scale_expansion_zeroelim(cxablen, cxab, cdx, cxxab);
		cyablen = scale_expansion_zeroelim(4, ab, cdy, cyab);
		cyyablen = scale_expansion_zeroelim(cyablen, cyab, cdy, cyyab);
		clen = fast_expansion_sum_zeroelim(cxxablen, cxxab, cyyablen, cyyab, cdet);

		ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
		finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, fin1);

		det = estimate(finlength, fin1);
		errbound = iccerrboundB * permanent;
		if ((det >= errbound) || (-det >= errbound)) {
			return det;
		}

		Two_Diff_Tail(pa[0], pd[0], adx, adxtail);
		Two_Diff_Tail(pa[1], pd[1], ady, adytail);
		Two_Diff_Tail(pb[0], pd[0], bdx, bdxtail);
		Two_Diff_Tail(pb[1], pd[1], bdy, bdytail);
		Two_Diff_Tail(pc[0], pd[0], cdx, cdxtail);
		Two_Diff_Tail(pc[1], pd[1], cdy, cdytail);
		if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0)
			&& (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0)) {
			return det;
		}

		errbound = iccerrboundC * permanent + resulterrbound * Absolute(det);
		det += ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail)
			- (bdy * cdxtail + cdx * bdytail))
			+ 2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx))
			+ ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail)
			- (cdy * adxtail + adx * cdytail))
			+ 2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx))
			+ ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail)
			- (ady * bdxtail + bdx * adytail))
			+ 2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx));
		if ((det >= errbound) || (-det >= errbound)) {
			return det;
		}

		finnow = fin1;
		finother = fin2;

		if ((bdxtail != 0.0) || (bdytail != 0.0)
			|| (cdxtail != 0.0) || (cdytail != 0.0)) {
			Square(adx, adxadx1, adxadx0);
			Square(ady, adyady1, adyady0);
			Two_Two_Sum(adxadx1, adxadx0, adyady1, adyady0, aa3, aa[2], aa[1], aa[0]);
			aa[3] = aa3;
		}
		if ((cdxtail != 0.0) || (cdytail != 0.0)
			|| (adxtail != 0.0) || (adytail != 0.0)) {
			Square(bdx, bdxbdx1, bdxbdx0);
			Square(bdy, bdybdy1, bdybdy0);
			Two_Two_Sum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0, bb3, bb[2], bb[1], bb[0]);
			bb[3] = bb3;
		}
		if ((adxtail != 0.0) || (adytail != 0.0)
			|| (bdxtail != 0.0) || (bdytail != 0.0)) {
			Square(cdx, cdxcdx1, cdxcdx0);
			Square(cdy, cdycdy1, cdycdy0);
			Two_Two_Sum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0, cc3, cc[2], cc[1], cc[0]);
			cc[3] = cc3;
		}

		if (adxtail != 0.0) {
			axtbclen = scale_expansion_zeroelim(4, bc, adxtail, axtbc);
			temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, 2.0 * adx,
				temp16a);

			axtcclen = scale_expansion_zeroelim(4, cc, adxtail, axtcc);
			temp16blen = scale_expansion_zeroelim(axtcclen, axtcc, bdy, temp16b);

			axtbblen = scale_expansion_zeroelim(4, bb, adxtail, axtbb);
			temp16clen = scale_expansion_zeroelim(axtbblen, axtbb, -cdy, temp16c);

			temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp16blen, temp16b, temp32a);
			temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
				temp32alen, temp32a, temp48);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
				temp48, finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
		if (adytail != 0.0) {
			aytbclen = scale_expansion_zeroelim(4, bc, adytail, aytbc);
			temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, 2.0 * ady,
				temp16a);

			aytbblen = scale_expansion_zeroelim(4, bb, adytail, aytbb);
			temp16blen = scale_expansion_zeroelim(aytbblen, aytbb, cdx, temp16b);

			aytcclen = scale_expansion_zeroelim(4, cc, adytail, aytcc);
			temp16clen = scale_expansion_zeroelim(aytcclen, aytcc, -bdx, temp16c);

			temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp16blen, temp16b, temp32a);
			temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
				temp32alen, temp32a, temp48);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
				temp48, finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
		if (bdxtail != 0.0) {
			bxtcalen = scale_expansion_zeroelim(4, ca, bdxtail, bxtca);
			temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, 2.0 * bdx,
				temp16a);

			bxtaalen = scale_expansion_zeroelim(4, aa, bdxtail, bxtaa);
			temp16blen = scale_expansion_zeroelim(bxtaalen, bxtaa, cdy, temp16b);

			bxtcclen = scale_expansion_zeroelim(4, cc, bdxtail, bxtcc);
			temp16clen = scale_expansion_zeroelim(bxtcclen, bxtcc, -ady, temp16c);

			temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp16blen, temp16b, temp32a);
			temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
				temp32alen, temp32a, temp48);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
				temp48, finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
		if (bdytail != 0.0) {
			bytcalen = scale_expansion_zeroelim(4, ca, bdytail, bytca);
			temp16alen = scale_expansion_zeroelim(bytcalen, bytca, 2.0 * bdy,
				temp16a);

			bytcclen = scale_expansion_zeroelim(4, cc, bdytail, bytcc);
			temp16blen = scale_expansion_zeroelim(bytcclen, bytcc, adx, temp16b);

			bytaalen = scale_expansion_zeroelim(4, aa, bdytail, bytaa);
			temp16clen = scale_expansion_zeroelim(bytaalen, bytaa, -cdx, temp16c);

			temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp16blen, temp16b, temp32a);
			temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
				temp32alen, temp32a, temp48);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
				temp48, finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
		if (cdxtail != 0.0) {
			cxtablen = scale_expansion_zeroelim(4, ab, cdxtail, cxtab);
			temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, 2.0 * cdx,
				temp16a);

			cxtbblen = scale_expansion_zeroelim(4, bb, cdxtail, cxtbb);
			temp16blen = scale_expansion_zeroelim(cxtbblen, cxtbb, ady, temp16b);

			cxtaalen = scale_expansion_zeroelim(4, aa, cdxtail, cxtaa);
			temp16clen = scale_expansion_zeroelim(cxtaalen, cxtaa, -bdy, temp16c);

			temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp16blen, temp16b, temp32a);
			temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
				temp32alen, temp32a, temp48);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
				temp48, finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
		if (cdytail != 0.0) {
			cytablen = scale_expansion_zeroelim(4, ab, cdytail, cytab);
			temp16alen = scale_expansion_zeroelim(cytablen, cytab, 2.0 * cdy,
				temp16a);

			cytaalen = scale_expansion_zeroelim(4, aa, cdytail, cytaa);
			temp16blen = scale_expansion_zeroelim(cytaalen, cytaa, bdx, temp16b);

			cytbblen = scale_expansion_zeroelim(4, bb, cdytail, cytbb);
			temp16clen = scale_expansion_zeroelim(cytbblen, cytbb, -adx, temp16c);

			temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp16blen, temp16b, temp32a);
			temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
				temp32alen, temp32a, temp48);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
				temp48, finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}

		if ((adxtail != 0.0) || (adytail != 0.0)) {
			if ((bdxtail != 0.0) || (bdytail != 0.0)
				|| (cdxtail != 0.0) || (cdytail != 0.0)) {
				Two_Product(bdxtail, cdy, ti1, ti0);
				Two_Product(bdx, cdytail, tj1, tj0);
				Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
				u[3] = u3;
				negate = -bdy;
				Two_Product(cdxtail, negate, ti1, ti0);
				negate = -bdytail;
				Two_Product(cdx, negate, tj1, tj0);
				Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
				v[3] = v3;
				bctlen = fast_expansion_sum_zeroelim(4, u, 4, v, bct);

				Two_Product(bdxtail, cdytail, ti1, ti0);
				Two_Product(cdxtail, bdytail, tj1, tj0);
				Two_Two_Diff(ti1, ti0, tj1, tj0, bctt3, bctt[2], bctt[1], bctt[0]);
				bctt[3] = bctt3;
				bcttlen = 4;
			}
			else {
				bct[0] = 0.0;
				bctlen = 1;
				bctt[0] = 0.0;
				bcttlen = 1;
			}

			if (adxtail != 0.0) {
				temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, adxtail, temp16a);
				axtbctlen = scale_expansion_zeroelim(bctlen, bct, adxtail, axtbct);
				temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, 2.0 * adx,
					temp32a);
				temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
					temp32alen, temp32a, temp48);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
					temp48, finother);
				finswap = finnow; finnow = finother; finother = finswap;
				if (bdytail != 0.0) {
					temp8len = scale_expansion_zeroelim(4, cc, adxtail, temp8);
					temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail,
						temp16a);
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
						temp16a, finother);
					finswap = finnow; finnow = finother; finother = finswap;
				}
				if (cdytail != 0.0) {
					temp8len = scale_expansion_zeroelim(4, bb, -adxtail, temp8);
					temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail,
						temp16a);
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
						temp16a, finother);
					finswap = finnow; finnow = finother; finother = finswap;
				}

				temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, adxtail,
					temp32a);
				axtbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adxtail, axtbctt);
				temp16alen = scale_expansion_zeroelim(axtbcttlen, axtbctt, 2.0 * adx,
					temp16a);
				temp16blen = scale_expansion_zeroelim(axtbcttlen, axtbctt, adxtail,
					temp16b);
				temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
					temp16blen, temp16b, temp32b);
				temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
					temp32blen, temp32b, temp64);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
					temp64, finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
			if (adytail != 0.0) {
				temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, adytail, temp16a);
				aytbctlen = scale_expansion_zeroelim(bctlen, bct, adytail, aytbct);
				temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, 2.0 * ady,
					temp32a);
				temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
					temp32alen, temp32a, temp48);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
					temp48, finother);
				finswap = finnow; finnow = finother; finother = finswap;


				temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, adytail,
					temp32a);
				aytbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adytail, aytbctt);
				temp16alen = scale_expansion_zeroelim(aytbcttlen, aytbctt, 2.0 * ady,
					temp16a);
				temp16blen = scale_expansion_zeroelim(aytbcttlen, aytbctt, adytail,
					temp16b);
				temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
					temp16blen, temp16b, temp32b);
				temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
					temp32blen, temp32b, temp64);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
					temp64, finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
		}
		if ((bdxtail != 0.0) || (bdytail != 0.0)) {
			if ((cdxtail != 0.0) || (cdytail != 0.0)
				|| (adxtail != 0.0) || (adytail != 0.0)) {
				Two_Product(cdxtail, ady, ti1, ti0);
				Two_Product(cdx, adytail, tj1, tj0);
				Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
				u[3] = u3;
				negate = -cdy;
				Two_Product(adxtail, negate, ti1, ti0);
				negate = -cdytail;
				Two_Product(adx, negate, tj1, tj0);
				Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
				v[3] = v3;
				catlen = fast_expansion_sum_zeroelim(4, u, 4, v, cat);

				Two_Product(cdxtail, adytail, ti1, ti0);
				Two_Product(adxtail, cdytail, tj1, tj0);
				Two_Two_Diff(ti1, ti0, tj1, tj0, catt3, catt[2], catt[1], catt[0]);
				catt[3] = catt3;
				cattlen = 4;
			}
			else {
				cat[0] = 0.0;
				catlen = 1;
				catt[0] = 0.0;
				cattlen = 1;
			}

			if (bdxtail != 0.0) {
				temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, bdxtail, temp16a);
				bxtcatlen = scale_expansion_zeroelim(catlen, cat, bdxtail, bxtcat);
				temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, 2.0 * bdx,
					temp32a);
				temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
					temp32alen, temp32a, temp48);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
					temp48, finother);
				finswap = finnow; finnow = finother; finother = finswap;
				if (cdytail != 0.0) {
					temp8len = scale_expansion_zeroelim(4, aa, bdxtail, temp8);
					temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail,
						temp16a);
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
						temp16a, finother);
					finswap = finnow; finnow = finother; finother = finswap;
				}
				if (adytail != 0.0) {
					temp8len = scale_expansion_zeroelim(4, cc, -bdxtail, temp8);
					temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail,
						temp16a);
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
						temp16a, finother);
					finswap = finnow; finnow = finother; finother = finswap;
				}

				temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, bdxtail,
					temp32a);
				bxtcattlen = scale_expansion_zeroelim(cattlen, catt, bdxtail, bxtcatt);
				temp16alen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, 2.0 * bdx,
					temp16a);
				temp16blen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, bdxtail,
					temp16b);
				temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
					temp16blen, temp16b, temp32b);
				temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
					temp32blen, temp32b, temp64);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
					temp64, finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
			if (bdytail != 0.0) {
				temp16alen = scale_expansion_zeroelim(bytcalen, bytca, bdytail, temp16a);
				bytcatlen = scale_expansion_zeroelim(catlen, cat, bdytail, bytcat);
				temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, 2.0 * bdy,
					temp32a);
				temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
					temp32alen, temp32a, temp48);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
					temp48, finother);
				finswap = finnow; finnow = finother; finother = finswap;


				temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, bdytail,
					temp32a);
				bytcattlen = scale_expansion_zeroelim(cattlen, catt, bdytail, bytcatt);
				temp16alen = scale_expansion_zeroelim(bytcattlen, bytcatt, 2.0 * bdy,
					temp16a);
				temp16blen = scale_expansion_zeroelim(bytcattlen, bytcatt, bdytail,
					temp16b);
				temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
					temp16blen, temp16b, temp32b);
				temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
					temp32blen, temp32b, temp64);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
					temp64, finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
		}
		if ((cdxtail != 0.0) || (cdytail != 0.0)) {
			if ((adxtail != 0.0) || (adytail != 0.0)
				|| (bdxtail != 0.0) || (bdytail != 0.0)) {
				Two_Product(adxtail, bdy, ti1, ti0);
				Two_Product(adx, bdytail, tj1, tj0);
				Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
				u[3] = u3;
				negate = -ady;
				Two_Product(bdxtail, negate, ti1, ti0);
				negate = -adytail;
				Two_Product(bdx, negate, tj1, tj0);
				Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
				v[3] = v3;
				abtlen = fast_expansion_sum_zeroelim(4, u, 4, v, abt);

				Two_Product(adxtail, bdytail, ti1, ti0);
				Two_Product(bdxtail, adytail, tj1, tj0);
				Two_Two_Diff(ti1, ti0, tj1, tj0, abtt3, abtt[2], abtt[1], abtt[0]);
				abtt[3] = abtt3;
				abttlen = 4;
			}
			else {
				abt[0] = 0.0;
				abtlen = 1;
				abtt[0] = 0.0;
				abttlen = 1;
			}

			if (cdxtail != 0.0) {
				temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, cdxtail, temp16a);
				cxtabtlen = scale_expansion_zeroelim(abtlen, abt, cdxtail, cxtabt);
				temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, 2.0 * cdx,
					temp32a);
				temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
					temp32alen, temp32a, temp48);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
					temp48, finother);
				finswap = finnow; finnow = finother; finother = finswap;
				if (adytail != 0.0) {
					temp8len = scale_expansion_zeroelim(4, bb, cdxtail, temp8);
					temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail,
						temp16a);
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
						temp16a, finother);
					finswap = finnow; finnow = finother; finother = finswap;
				}
				if (bdytail != 0.0) {
					temp8len = scale_expansion_zeroelim(4, aa, -cdxtail, temp8);
					temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail,
						temp16a);
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
						temp16a, finother);
					finswap = finnow; finnow = finother; finother = finswap;
				}

				temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, cdxtail,
					temp32a);
				cxtabttlen = scale_expansion_zeroelim(abttlen, abtt, cdxtail, cxtabtt);
				temp16alen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, 2.0 * cdx,
					temp16a);
				temp16blen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, cdxtail,
					temp16b);
				temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
					temp16blen, temp16b, temp32b);
				temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
					temp32blen, temp32b, temp64);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
					temp64, finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
			if (cdytail != 0.0) {
				temp16alen = scale_expansion_zeroelim(cytablen, cytab, cdytail, temp16a);
				cytabtlen = scale_expansion_zeroelim(abtlen, abt, cdytail, cytabt);
				temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, 2.0 * cdy,
					temp32a);
				temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
					temp32alen, temp32a, temp48);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
					temp48, finother);
				finswap = finnow; finnow = finother; finother = finswap;


				temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, cdytail,
					temp32a);
				cytabttlen = scale_expansion_zeroelim(abttlen, abtt, cdytail, cytabtt);
				temp16alen = scale_expansion_zeroelim(cytabttlen, cytabtt, 2.0 * cdy,
					temp16a);
				temp16blen = scale_expansion_zeroelim(cytabttlen, cytabtt, cdytail,
					temp16b);
				temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
					temp16blen, temp16b, temp32b);
				temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
					temp32blen, temp32b, temp64);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
					temp64, finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
		}

		return finnow[finlength - 1];
	}

	REAL incircle(REAL *pa, REAL *pb, REAL *pc, REAL *pd)
	{
		REAL adx, bdx, cdx, ady, bdy, cdy;
		REAL bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady;
		REAL alift, blift, clift;
		REAL det;
		REAL permanent, errbound;

		adx = pa[0] - pd[0];
		bdx = pb[0] - pd[0];
		cdx = pc[0] - pd[0];
		ady = pa[1] - pd[1];
		bdy = pb[1] - pd[1];
		cdy = pc[1] - pd[1];

		bdxcdy = bdx * cdy;
		cdxbdy = cdx * bdy;
		alift = adx * adx + ady * ady;

		cdxady = cdx * ady;
		adxcdy = adx * cdy;
		blift = bdx * bdx + bdy * bdy;

		adxbdy = adx * bdy;
		bdxady = bdx * ady;
		clift = cdx * cdx + cdy * cdy;

		det = alift * (bdxcdy - cdxbdy)
			+ blift * (cdxady - adxcdy)
			+ clift * (adxbdy - bdxady);

		permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * alift
			+ (Absolute(cdxady) + Absolute(adxcdy)) * blift
			+ (Absolute(adxbdy) + Absolute(bdxady)) * clift;
		errbound = iccerrboundA * permanent;
		if ((det > errbound) || (-det > errbound)) {
			return det;
		}

		return incircleadapt(pa, pb, pc, pd, permanent);
	}

	/*****************************************************************************/
	/*                                                                           */
	/*  inspherefast()   Approximate 3D insphere test.  Nonrobust.               */
	/*  insphereexact()   Exact 3D insphere test.  Robust.                       */
	/*  insphereslow()   Another exact 3D insphere test.  Robust.                */
	/*  insphere()   Adaptive exact 3D insphere test.  Robust.                   */
	/*                                                                           */
	/*               Return a positive value if the point pe lies inside the     */
	/*               sphere passing through pa, pb, pc, and pd; a negative value */
	/*               if it lies outside; and zero if the five points are         */
	/*               cospherical.  The points pa, pb, pc, and pd must be ordered */
	/*               so that they have a positive orientation (as defined by     */
	/*               orient3d()), or the sign of the result will be reversed.    */
	/*                                                                           */
	/*  Only the first and last routine should be used; the middle two are for   */
	/*  timings.                                                                 */
	/*                                                                           */
	/*  The last three use exact arithmetic to ensure a correct answer.  The     */
	/*  result returned is the determinant of a matrix.  In insphere() only,     */
	/*  this determinant is computed adaptively, in the sense that exact         */
	/*  arithmetic is used only to the degree it is needed to ensure that the    */
	/*  returned value has the correct sign.  Hence, insphere() is usually quite */
	/*  fast, but will run more slowly when the input points are cospherical or  */
	/*  nearly so.                                                               */
	/*                                                                           */
	/*****************************************************************************/

	REAL inspherefast(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe)
	{
		REAL aex, bex, cex, dex;
		REAL aey, bey, cey, dey;
		REAL aez, bez, cez, dez;
		REAL alift, blift, clift, dlift;
		REAL ab, bc, cd, da, ac, bd;
		REAL abc, bcd, cda, dab;

		aex = pa[0] - pe[0];
		bex = pb[0] - pe[0];
		cex = pc[0] - pe[0];
		dex = pd[0] - pe[0];
		aey = pa[1] - pe[1];
		bey = pb[1] - pe[1];
		cey = pc[1] - pe[1];
		dey = pd[1] - pe[1];
		aez = pa[2] - pe[2];
		bez = pb[2] - pe[2];
		cez = pc[2] - pe[2];
		dez = pd[2] - pe[2];

		ab = aex * bey - bex * aey;
		bc = bex * cey - cex * bey;
		cd = cex * dey - dex * cey;
		da = dex * aey - aex * dey;

		ac = aex * cey - cex * aey;
		bd = bex * dey - dex * bey;

		abc = aez * bc - bez * ac + cez * ab;
		bcd = bez * cd - cez * bd + dez * bc;
		cda = cez * da + dez * ac + aez * cd;
		dab = dez * ab + aez * bd + bez * da;

		alift = aex * aex + aey * aey + aez * aez;
		blift = bex * bex + bey * bey + bez * bez;
		clift = cex * cex + cey * cey + cez * cez;
		dlift = dex * dex + dey * dey + dez * dez;

		return (dlift * abc - clift * dab) + (blift * cda - alift * bcd);
	}

	REAL insphereexact(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe)
	{
		INEXACT REAL axby1, bxcy1, cxdy1, dxey1, exay1;
		INEXACT REAL bxay1, cxby1, dxcy1, exdy1, axey1;
		INEXACT REAL axcy1, bxdy1, cxey1, dxay1, exby1;
		INEXACT REAL cxay1, dxby1, excy1, axdy1, bxey1;
		REAL axby0, bxcy0, cxdy0, dxey0, exay0;
		REAL bxay0, cxby0, dxcy0, exdy0, axey0;
		REAL axcy0, bxdy0, cxey0, dxay0, exby0;
		REAL cxay0, dxby0, excy0, axdy0, bxey0;
		REAL ab[4], bc[4], cd[4], de[4], ea[4];
		REAL ac[4], bd[4], ce[4], da[4], eb[4];
		REAL temp8a[8], temp8b[8], temp16[16];
		int temp8alen, temp8blen, temp16len;
		REAL abc[24], bcd[24], cde[24], dea[24], eab[24];
		REAL abd[24], bce[24], cda[24], deb[24], eac[24];
		int abclen, bcdlen, cdelen, dealen, eablen;
		int abdlen, bcelen, cdalen, deblen, eaclen;
		REAL temp48a[48], temp48b[48];
		int temp48alen, temp48blen;
		REAL abcd[96], bcde[96], cdea[96], deab[96], eabc[96];
		int abcdlen, bcdelen, cdealen, deablen, eabclen;
		REAL temp192[192];
		REAL det384x[384], det384y[384], det384z[384];
		int xlen, ylen, zlen;
		REAL detxy[768];
		int xylen;
		REAL adet[1152], bdet[1152], cdet[1152], ddet[1152], edet[1152];
		int alen, blen, clen, dlen, elen;
		REAL abdet[2304], cddet[2304], cdedet[3456];
		int ablen, cdlen;
		REAL deter[5760];
		int deterlen;
		int i;

		INEXACT REAL bvirt;
		REAL avirt, bround, around;
		INEXACT REAL c;
		INEXACT REAL abig;
		REAL ahi, alo, bhi, blo;
		REAL err1, err2, err3;
		INEXACT REAL _i, _j;
		REAL _0;


		Two_Product(pa[0], pb[1], axby1, axby0);
		Two_Product(pb[0], pa[1], bxay1, bxay0);
		Two_Two_Diff(axby1, axby0, bxay1, bxay0, ab[3], ab[2], ab[1], ab[0]);

		Two_Product(pb[0], pc[1], bxcy1, bxcy0);
		Two_Product(pc[0], pb[1], cxby1, cxby0);
		Two_Two_Diff(bxcy1, bxcy0, cxby1, cxby0, bc[3], bc[2], bc[1], bc[0]);

		Two_Product(pc[0], pd[1], cxdy1, cxdy0);
		Two_Product(pd[0], pc[1], dxcy1, dxcy0);
		Two_Two_Diff(cxdy1, cxdy0, dxcy1, dxcy0, cd[3], cd[2], cd[1], cd[0]);

		Two_Product(pd[0], pe[1], dxey1, dxey0);
		Two_Product(pe[0], pd[1], exdy1, exdy0);
		Two_Two_Diff(dxey1, dxey0, exdy1, exdy0, de[3], de[2], de[1], de[0]);

		Two_Product(pe[0], pa[1], exay1, exay0);
		Two_Product(pa[0], pe[1], axey1, axey0);
		Two_Two_Diff(exay1, exay0, axey1, axey0, ea[3], ea[2], ea[1], ea[0]);

		Two_Product(pa[0], pc[1], axcy1, axcy0);
		Two_Product(pc[0], pa[1], cxay1, cxay0);
		Two_Two_Diff(axcy1, axcy0, cxay1, cxay0, ac[3], ac[2], ac[1], ac[0]);

		Two_Product(pb[0], pd[1], bxdy1, bxdy0);
		Two_Product(pd[0], pb[1], dxby1, dxby0);
		Two_Two_Diff(bxdy1, bxdy0, dxby1, dxby0, bd[3], bd[2], bd[1], bd[0]);

		Two_Product(pc[0], pe[1], cxey1, cxey0);
		Two_Product(pe[0], pc[1], excy1, excy0);
		Two_Two_Diff(cxey1, cxey0, excy1, excy0, ce[3], ce[2], ce[1], ce[0]);

		Two_Product(pd[0], pa[1], dxay1, dxay0);
		Two_Product(pa[0], pd[1], axdy1, axdy0);
		Two_Two_Diff(dxay1, dxay0, axdy1, axdy0, da[3], da[2], da[1], da[0]);

		Two_Product(pe[0], pb[1], exby1, exby0);
		Two_Product(pb[0], pe[1], bxey1, bxey0);
		Two_Two_Diff(exby1, exby0, bxey1, bxey0, eb[3], eb[2], eb[1], eb[0]);

		temp8alen = scale_expansion_zeroelim(4, bc, pa[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, ac, -pb[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
			temp16);
		temp8alen = scale_expansion_zeroelim(4, ab, pc[2], temp8a);
		abclen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
			abc);

		temp8alen = scale_expansion_zeroelim(4, cd, pb[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, bd, -pc[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
			temp16);
		temp8alen = scale_expansion_zeroelim(4, bc, pd[2], temp8a);
		bcdlen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
			bcd);

		temp8alen = scale_expansion_zeroelim(4, de, pc[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, ce, -pd[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
			temp16);
		temp8alen = scale_expansion_zeroelim(4, cd, pe[2], temp8a);
		cdelen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
			cde);

		temp8alen = scale_expansion_zeroelim(4, ea, pd[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, da, -pe[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
			temp16);
		temp8alen = scale_expansion_zeroelim(4, de, pa[2], temp8a);
		dealen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
			dea);

		temp8alen = scale_expansion_zeroelim(4, ab, pe[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, eb, -pa[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
			temp16);
		temp8alen = scale_expansion_zeroelim(4, ea, pb[2], temp8a);
		eablen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
			eab);

		temp8alen = scale_expansion_zeroelim(4, bd, pa[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, da, pb[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
			temp16);
		temp8alen = scale_expansion_zeroelim(4, ab, pd[2], temp8a);
		abdlen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
			abd);

		temp8alen = scale_expansion_zeroelim(4, ce, pb[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, eb, pc[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
			temp16);
		temp8alen = scale_expansion_zeroelim(4, bc, pe[2], temp8a);
		bcelen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
			bce);

		temp8alen = scale_expansion_zeroelim(4, da, pc[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, ac, pd[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
			temp16);
		temp8alen = scale_expansion_zeroelim(4, cd, pa[2], temp8a);
		cdalen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
			cda);

		temp8alen = scale_expansion_zeroelim(4, eb, pd[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, bd, pe[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
			temp16);
		temp8alen = scale_expansion_zeroelim(4, de, pb[2], temp8a);
		deblen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
			deb);

		temp8alen = scale_expansion_zeroelim(4, ac, pe[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, ce, pa[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
			temp16);
		temp8alen = scale_expansion_zeroelim(4, ea, pc[2], temp8a);
		eaclen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
			eac);

		temp48alen = fast_expansion_sum_zeroelim(cdelen, cde, bcelen, bce, temp48a);
		temp48blen = fast_expansion_sum_zeroelim(deblen, deb, bcdlen, bcd, temp48b);
		for (i = 0; i < temp48blen; i++) {
			temp48b[i] = -temp48b[i];
		}
		bcdelen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
			temp48blen, temp48b, bcde);
		xlen = scale_expansion_zeroelim(bcdelen, bcde, pa[0], temp192);
		xlen = scale_expansion_zeroelim(xlen, temp192, pa[0], det384x);
		ylen = scale_expansion_zeroelim(bcdelen, bcde, pa[1], temp192);
		ylen = scale_expansion_zeroelim(ylen, temp192, pa[1], det384y);
		zlen = scale_expansion_zeroelim(bcdelen, bcde, pa[2], temp192);
		zlen = scale_expansion_zeroelim(zlen, temp192, pa[2], det384z);
		xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
		alen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, adet);

		temp48alen = fast_expansion_sum_zeroelim(dealen, dea, cdalen, cda, temp48a);
		temp48blen = fast_expansion_sum_zeroelim(eaclen, eac, cdelen, cde, temp48b);
		for (i = 0; i < temp48blen; i++) {
			temp48b[i] = -temp48b[i];
		}
		cdealen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
			temp48blen, temp48b, cdea);
		xlen = scale_expansion_zeroelim(cdealen, cdea, pb[0], temp192);
		xlen = scale_expansion_zeroelim(xlen, temp192, pb[0], det384x);
		ylen = scale_expansion_zeroelim(cdealen, cdea, pb[1], temp192);
		ylen = scale_expansion_zeroelim(ylen, temp192, pb[1], det384y);
		zlen = scale_expansion_zeroelim(cdealen, cdea, pb[2], temp192);
		zlen = scale_expansion_zeroelim(zlen, temp192, pb[2], det384z);
		xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
		blen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, bdet);

		temp48alen = fast_expansion_sum_zeroelim(eablen, eab, deblen, deb, temp48a);
		temp48blen = fast_expansion_sum_zeroelim(abdlen, abd, dealen, dea, temp48b);
		for (i = 0; i < temp48blen; i++) {
			temp48b[i] = -temp48b[i];
		}
		deablen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
			temp48blen, temp48b, deab);
		xlen = scale_expansion_zeroelim(deablen, deab, pc[0], temp192);
		xlen = scale_expansion_zeroelim(xlen, temp192, pc[0], det384x);
		ylen = scale_expansion_zeroelim(deablen, deab, pc[1], temp192);
		ylen = scale_expansion_zeroelim(ylen, temp192, pc[1], det384y);
		zlen = scale_expansion_zeroelim(deablen, deab, pc[2], temp192);
		zlen = scale_expansion_zeroelim(zlen, temp192, pc[2], det384z);
		xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
		clen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, cdet);

		temp48alen = fast_expansion_sum_zeroelim(abclen, abc, eaclen, eac, temp48a);
		temp48blen = fast_expansion_sum_zeroelim(bcelen, bce, eablen, eab, temp48b);
		for (i = 0; i < temp48blen; i++) {
			temp48b[i] = -temp48b[i];
		}
		eabclen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
			temp48blen, temp48b, eabc);
		xlen = scale_expansion_zeroelim(eabclen, eabc, pd[0], temp192);
		xlen = scale_expansion_zeroelim(xlen, temp192, pd[0], det384x);
		ylen = scale_expansion_zeroelim(eabclen, eabc, pd[1], temp192);
		ylen = scale_expansion_zeroelim(ylen, temp192, pd[1], det384y);
		zlen = scale_expansion_zeroelim(eabclen, eabc, pd[2], temp192);
		zlen = scale_expansion_zeroelim(zlen, temp192, pd[2], det384z);
		xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
		dlen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, ddet);

		temp48alen = fast_expansion_sum_zeroelim(bcdlen, bcd, abdlen, abd, temp48a);
		temp48blen = fast_expansion_sum_zeroelim(cdalen, cda, abclen, abc, temp48b);
		for (i = 0; i < temp48blen; i++) {
			temp48b[i] = -temp48b[i];
		}
		abcdlen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
			temp48blen, temp48b, abcd);
		xlen = scale_expansion_zeroelim(abcdlen, abcd, pe[0], temp192);
		xlen = scale_expansion_zeroelim(xlen, temp192, pe[0], det384x);
		ylen = scale_expansion_zeroelim(abcdlen, abcd, pe[1], temp192);
		ylen = scale_expansion_zeroelim(ylen, temp192, pe[1], det384y);
		zlen = scale_expansion_zeroelim(abcdlen, abcd, pe[2], temp192);
		zlen = scale_expansion_zeroelim(zlen, temp192, pe[2], det384z);
		xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
		elen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, edet);

		ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
		cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
		cdelen = fast_expansion_sum_zeroelim(cdlen, cddet, elen, edet, cdedet);
		deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdelen, cdedet, deter);

		return deter[deterlen - 1];
	}

	REAL insphereslow(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe)
	{
		INEXACT REAL aex, bex, cex, dex, aey, bey, cey, dey, aez, bez, cez, dez;
		REAL aextail, bextail, cextail, dextail;
		REAL aeytail, beytail, ceytail, deytail;
		REAL aeztail, beztail, ceztail, deztail;
		REAL negate, negatetail;
		INEXACT REAL axby7, bxcy7, cxdy7, dxay7, axcy7, bxdy7;
		INEXACT REAL bxay7, cxby7, dxcy7, axdy7, cxay7, dxby7;
		REAL axby[8], bxcy[8], cxdy[8], dxay[8], axcy[8], bxdy[8];
		REAL bxay[8], cxby[8], dxcy[8], axdy[8], cxay[8], dxby[8];
		REAL ab[16], bc[16], cd[16], da[16], ac[16], bd[16];
		int ablen, bclen, cdlen, dalen, aclen, bdlen;
		REAL temp32a[32], temp32b[32], temp64a[64], temp64b[64], temp64c[64];
		int temp32alen, temp32blen, temp64alen, temp64blen, temp64clen;
		REAL temp128[128], temp192[192];
		int temp128len, temp192len;
		REAL detx[384], detxx[768], detxt[384], detxxt[768], detxtxt[768];
		int xlen, xxlen, xtlen, xxtlen, xtxtlen;
		REAL x1[1536], x2[2304];
		int x1len, x2len;
		REAL dety[384], detyy[768], detyt[384], detyyt[768], detytyt[768];
		int ylen, yylen, ytlen, yytlen, ytytlen;
		REAL y1[1536], y2[2304];
		int y1len, y2len;
		REAL detz[384], detzz[768], detzt[384], detzzt[768], detztzt[768];
		int zlen, zzlen, ztlen, zztlen, ztztlen;
		REAL z1[1536], z2[2304];
		int z1len, z2len;
		REAL detxy[4608];
		int xylen;
		REAL adet[6912], bdet[6912], cdet[6912], ddet[6912];
		int alen, blen, clen, dlen;
		REAL abdet[13824], cddet[13824], deter[27648];
		int deterlen;
		int i;

		INEXACT REAL bvirt;
		REAL avirt, bround, around;
		INEXACT REAL c;
		INEXACT REAL abig;
		REAL a0hi, a0lo, a1hi, a1lo, bhi, blo;
		REAL err1, err2, err3;
		INEXACT REAL _i, _j, _k, _l, _m, _n;
		REAL _0, _1, _2;

		Two_Diff(pa[0], pe[0], aex, aextail);
		Two_Diff(pa[1], pe[1], aey, aeytail);
		Two_Diff(pa[2], pe[2], aez, aeztail);
		Two_Diff(pb[0], pe[0], bex, bextail);
		Two_Diff(pb[1], pe[1], bey, beytail);
		Two_Diff(pb[2], pe[2], bez, beztail);
		Two_Diff(pc[0], pe[0], cex, cextail);
		Two_Diff(pc[1], pe[1], cey, ceytail);
		Two_Diff(pc[2], pe[2], cez, ceztail);
		Two_Diff(pd[0], pe[0], dex, dextail);
		Two_Diff(pd[1], pe[1], dey, deytail);
		Two_Diff(pd[2], pe[2], dez, deztail);

		Two_Two_Product(aex, aextail, bey, beytail,
			axby7, axby[6], axby[5], axby[4],
			axby[3], axby[2], axby[1], axby[0]);
		axby[7] = axby7;
		negate = -aey;
		negatetail = -aeytail;
		Two_Two_Product(bex, bextail, negate, negatetail,
			bxay7, bxay[6], bxay[5], bxay[4],
			bxay[3], bxay[2], bxay[1], bxay[0]);
		bxay[7] = bxay7;
		ablen = fast_expansion_sum_zeroelim(8, axby, 8, bxay, ab);
		Two_Two_Product(bex, bextail, cey, ceytail,
			bxcy7, bxcy[6], bxcy[5], bxcy[4],
			bxcy[3], bxcy[2], bxcy[1], bxcy[0]);
		bxcy[7] = bxcy7;
		negate = -bey;
		negatetail = -beytail;
		Two_Two_Product(cex, cextail, negate, negatetail,
			cxby7, cxby[6], cxby[5], cxby[4],
			cxby[3], cxby[2], cxby[1], cxby[0]);
		cxby[7] = cxby7;
		bclen = fast_expansion_sum_zeroelim(8, bxcy, 8, cxby, bc);
		Two_Two_Product(cex, cextail, dey, deytail,
			cxdy7, cxdy[6], cxdy[5], cxdy[4],
			cxdy[3], cxdy[2], cxdy[1], cxdy[0]);
		cxdy[7] = cxdy7;
		negate = -cey;
		negatetail = -ceytail;
		Two_Two_Product(dex, dextail, negate, negatetail,
			dxcy7, dxcy[6], dxcy[5], dxcy[4],
			dxcy[3], dxcy[2], dxcy[1], dxcy[0]);
		dxcy[7] = dxcy7;
		cdlen = fast_expansion_sum_zeroelim(8, cxdy, 8, dxcy, cd);
		Two_Two_Product(dex, dextail, aey, aeytail,
			dxay7, dxay[6], dxay[5], dxay[4],
			dxay[3], dxay[2], dxay[1], dxay[0]);
		dxay[7] = dxay7;
		negate = -dey;
		negatetail = -deytail;
		Two_Two_Product(aex, aextail, negate, negatetail,
			axdy7, axdy[6], axdy[5], axdy[4],
			axdy[3], axdy[2], axdy[1], axdy[0]);
		axdy[7] = axdy7;
		dalen = fast_expansion_sum_zeroelim(8, dxay, 8, axdy, da);
		Two_Two_Product(aex, aextail, cey, ceytail,
			axcy7, axcy[6], axcy[5], axcy[4],
			axcy[3], axcy[2], axcy[1], axcy[0]);
		axcy[7] = axcy7;
		negate = -aey;
		negatetail = -aeytail;
		Two_Two_Product(cex, cextail, negate, negatetail,
			cxay7, cxay[6], cxay[5], cxay[4],
			cxay[3], cxay[2], cxay[1], cxay[0]);
		cxay[7] = cxay7;
		aclen = fast_expansion_sum_zeroelim(8, axcy, 8, cxay, ac);
		Two_Two_Product(bex, bextail, dey, deytail,
			bxdy7, bxdy[6], bxdy[5], bxdy[4],
			bxdy[3], bxdy[2], bxdy[1], bxdy[0]);
		bxdy[7] = bxdy7;
		negate = -bey;
		negatetail = -beytail;
		Two_Two_Product(dex, dextail, negate, negatetail,
			dxby7, dxby[6], dxby[5], dxby[4],
			dxby[3], dxby[2], dxby[1], dxby[0]);
		dxby[7] = dxby7;
		bdlen = fast_expansion_sum_zeroelim(8, bxdy, 8, dxby, bd);

		temp32alen = scale_expansion_zeroelim(cdlen, cd, -bez, temp32a);
		temp32blen = scale_expansion_zeroelim(cdlen, cd, -beztail, temp32b);
		temp64alen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
			temp32blen, temp32b, temp64a);
		temp32alen = scale_expansion_zeroelim(bdlen, bd, cez, temp32a);
		temp32blen = scale_expansion_zeroelim(bdlen, bd, ceztail, temp32b);
		temp64blen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
			temp32blen, temp32b, temp64b);
		temp32alen = scale_expansion_zeroelim(bclen, bc, -dez, temp32a);
		temp32blen = scale_expansion_zeroelim(bclen, bc, -deztail, temp32b);
		temp64clen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
			temp32blen, temp32b, temp64c);
		temp128len = fast_expansion_sum_zeroelim(temp64alen, temp64a,
			temp64blen, temp64b, temp128);
		temp192len = fast_expansion_sum_zeroelim(temp64clen, temp64c,
			temp128len, temp128, temp192);
		xlen = scale_expansion_zeroelim(temp192len, temp192, aex, detx);
		xxlen = scale_expansion_zeroelim(xlen, detx, aex, detxx);
		xtlen = scale_expansion_zeroelim(temp192len, temp192, aextail, detxt);
		xxtlen = scale_expansion_zeroelim(xtlen, detxt, aex, detxxt);
		for (i = 0; i < xxtlen; i++) {
			detxxt[i] *= 2.0;
		}
		xtxtlen = scale_expansion_zeroelim(xtlen, detxt, aextail, detxtxt);
		x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1);
		x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2);
		ylen = scale_expansion_zeroelim(temp192len, temp192, aey, dety);
		yylen = scale_expansion_zeroelim(ylen, dety, aey, detyy);
		ytlen = scale_expansion_zeroelim(temp192len, temp192, aeytail, detyt);
		yytlen = scale_expansion_zeroelim(ytlen, detyt, aey, detyyt);
		for (i = 0; i < yytlen; i++) {
			detyyt[i] *= 2.0;
		}
		ytytlen = scale_expansion_zeroelim(ytlen, detyt, aeytail, detytyt);
		y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1);
		y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2);
		zlen = scale_expansion_zeroelim(temp192len, temp192, aez, detz);
		zzlen = scale_expansion_zeroelim(zlen, detz, aez, detzz);
		ztlen = scale_expansion_zeroelim(temp192len, temp192, aeztail, detzt);
		zztlen = scale_expansion_zeroelim(ztlen, detzt, aez, detzzt);
		for (i = 0; i < zztlen; i++) {
			detzzt[i] *= 2.0;
		}
		ztztlen = scale_expansion_zeroelim(ztlen, detzt, aeztail, detztzt);
		z1len = fast_expansion_sum_zeroelim(zzlen, detzz, zztlen, detzzt, z1);
		z2len = fast_expansion_sum_zeroelim(z1len, z1, ztztlen, detztzt, z2);
		xylen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, detxy);
		alen = fast_expansion_sum_zeroelim(z2len, z2, xylen, detxy, adet);

		temp32alen = scale_expansion_zeroelim(dalen, da, cez, temp32a);
		temp32blen = scale_expansion_zeroelim(dalen, da, ceztail, temp32b);
		temp64alen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
			temp32blen, temp32b, temp64a);
		temp32alen = scale_expansion_zeroelim(aclen, ac, dez, temp32a);
		temp32blen = scale_expansion_zeroelim(aclen, ac, deztail, temp32b);
		temp64blen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
			temp32blen, temp32b, temp64b);
		temp32alen = scale_expansion_zeroelim(cdlen, cd, aez, temp32a);
		temp32blen = scale_expansion_zeroelim(cdlen, cd, aeztail, temp32b);
		temp64clen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
			temp32blen, temp32b, temp64c);
		temp128len = fast_expansion_sum_zeroelim(temp64alen, temp64a,
			temp64blen, temp64b, temp128);
		temp192len = fast_expansion_sum_zeroelim(temp64clen, temp64c,
			temp128len, temp128, temp192);
		xlen = scale_expansion_zeroelim(temp192len, temp192, bex, detx);
		xxlen = scale_expansion_zeroelim(xlen, detx, bex, detxx);
		xtlen = scale_expansion_zeroelim(temp192len, temp192, bextail, detxt);
		xxtlen = scale_expansion_zeroelim(xtlen, detxt, bex, detxxt);
		for (i = 0; i < xxtlen; i++) {
			detxxt[i] *= 2.0;
		}
		xtxtlen = scale_expansion_zeroelim(xtlen, detxt, bextail, detxtxt);
		x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1);
		x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2);
		ylen = scale_expansion_zeroelim(temp192len, temp192, bey, dety);
		yylen = scale_expansion_zeroelim(ylen, dety, bey, detyy);
		ytlen = scale_expansion_zeroelim(temp192len, temp192, beytail, detyt);
		yytlen = scale_expansion_zeroelim(ytlen, detyt, bey, detyyt);
		for (i = 0; i < yytlen; i++) {
			detyyt[i] *= 2.0;
		}
		ytytlen = scale_expansion_zeroelim(ytlen, detyt, beytail, detytyt);
		y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1);
		y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2);
		zlen = scale_expansion_zeroelim(temp192len, temp192, bez, detz);
		zzlen = scale_expansion_zeroelim(zlen, detz, bez, detzz);
		ztlen = scale_expansion_zeroelim(temp192len, temp192, beztail, detzt);
		zztlen = scale_expansion_zeroelim(ztlen, detzt, bez, detzzt);
		for (i = 0; i < zztlen; i++) {
			detzzt[i] *= 2.0;
		}
		ztztlen = scale_expansion_zeroelim(ztlen, detzt, beztail, detztzt);
		z1len = fast_expansion_sum_zeroelim(zzlen, detzz, zztlen, detzzt, z1);
		z2len = fast_expansion_sum_zeroelim(z1len, z1, ztztlen, detztzt, z2);
		xylen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, detxy);
		blen = fast_expansion_sum_zeroelim(z2len, z2, xylen, detxy, bdet);

		temp32alen = scale_expansion_zeroelim(ablen, ab, -dez, temp32a);
		temp32blen = scale_expansion_zeroelim(ablen, ab, -deztail, temp32b);
		temp64alen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
			temp32blen, temp32b, temp64a);
		temp32alen = scale_expansion_zeroelim(bdlen, bd, -aez, temp32a);
		temp32blen = scale_expansion_zeroelim(bdlen, bd, -aeztail, temp32b);
		temp64blen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
			temp32blen, temp32b, temp64b);
		temp32alen = scale_expansion_zeroelim(dalen, da, -bez, temp32a);
		temp32blen = scale_expansion_zeroelim(dalen, da, -beztail, temp32b);
		temp64clen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
			temp32blen, temp32b, temp64c);
		temp128len = fast_expansion_sum_zeroelim(temp64alen, temp64a,
			temp64blen, temp64b, temp128);
		temp192len = fast_expansion_sum_zeroelim(temp64clen, temp64c,
			temp128len, temp128, temp192);
		xlen = scale_expansion_zeroelim(temp192len, temp192, cex, detx);
		xxlen = scale_expansion_zeroelim(xlen, detx, cex, detxx);
		xtlen = scale_expansion_zeroelim(temp192len, temp192, cextail, detxt);
		xxtlen = scale_expansion_zeroelim(xtlen, detxt, cex, detxxt);
		for (i = 0; i < xxtlen; i++) {
			detxxt[i] *= 2.0;
		}
		xtxtlen = scale_expansion_zeroelim(xtlen, detxt, cextail, detxtxt);
		x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1);
		x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2);
		ylen = scale_expansion_zeroelim(temp192len, temp192, cey, dety);
		yylen = scale_expansion_zeroelim(ylen, dety, cey, detyy);
		ytlen = scale_expansion_zeroelim(temp192len, temp192, ceytail, detyt);
		yytlen = scale_expansion_zeroelim(ytlen, detyt, cey, detyyt);
		for (i = 0; i < yytlen; i++) {
			detyyt[i] *= 2.0;
		}
		ytytlen = scale_expansion_zeroelim(ytlen, detyt, ceytail, detytyt);
		y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1);
		y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2);
		zlen = scale_expansion_zeroelim(temp192len, temp192, cez, detz);
		zzlen = scale_expansion_zeroelim(zlen, detz, cez, detzz);
		ztlen = scale_expansion_zeroelim(temp192len, temp192, ceztail, detzt);
		zztlen = scale_expansion_zeroelim(ztlen, detzt, cez, detzzt);
		for (i = 0; i < zztlen; i++) {
			detzzt[i] *= 2.0;
		}
		ztztlen = scale_expansion_zeroelim(ztlen, detzt, ceztail, detztzt);
		z1len = fast_expansion_sum_zeroelim(zzlen, detzz, zztlen, detzzt, z1);
		z2len = fast_expansion_sum_zeroelim(z1len, z1, ztztlen, detztzt, z2);
		xylen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, detxy);
		clen = fast_expansion_sum_zeroelim(z2len, z2, xylen, detxy, cdet);

		temp32alen = scale_expansion_zeroelim(bclen, bc, aez, temp32a);
		temp32blen = scale_expansion_zeroelim(bclen, bc, aeztail, temp32b);
		temp64alen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
			temp32blen, temp32b, temp64a);
		temp32alen = scale_expansion_zeroelim(aclen, ac, -bez, temp32a);
		temp32blen = scale_expansion_zeroelim(aclen, ac, -beztail, temp32b);
		temp64blen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
			temp32blen, temp32b, temp64b);
		temp32alen = scale_expansion_zeroelim(ablen, ab, cez, temp32a);
		temp32blen = scale_expansion_zeroelim(ablen, ab, ceztail, temp32b);
		temp64clen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
			temp32blen, temp32b, temp64c);
		temp128len = fast_expansion_sum_zeroelim(temp64alen, temp64a,
			temp64blen, temp64b, temp128);
		temp192len = fast_expansion_sum_zeroelim(temp64clen, temp64c,
			temp128len, temp128, temp192);
		xlen = scale_expansion_zeroelim(temp192len, temp192, dex, detx);
		xxlen = scale_expansion_zeroelim(xlen, detx, dex, detxx);
		xtlen = scale_expansion_zeroelim(temp192len, temp192, dextail, detxt);
		xxtlen = scale_expansion_zeroelim(xtlen, detxt, dex, detxxt);
		for (i = 0; i < xxtlen; i++) {
			detxxt[i] *= 2.0;
		}
		xtxtlen = scale_expansion_zeroelim(xtlen, detxt, dextail, detxtxt);
		x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1);
		x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2);
		ylen = scale_expansion_zeroelim(temp192len, temp192, dey, dety);
		yylen = scale_expansion_zeroelim(ylen, dety, dey, detyy);
		ytlen = scale_expansion_zeroelim(temp192len, temp192, deytail, detyt);
		yytlen = scale_expansion_zeroelim(ytlen, detyt, dey, detyyt);
		for (i = 0; i < yytlen; i++) {
			detyyt[i] *= 2.0;
		}
		ytytlen = scale_expansion_zeroelim(ytlen, detyt, deytail, detytyt);
		y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1);
		y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2);
		zlen = scale_expansion_zeroelim(temp192len, temp192, dez, detz);
		zzlen = scale_expansion_zeroelim(zlen, detz, dez, detzz);
		ztlen = scale_expansion_zeroelim(temp192len, temp192, deztail, detzt);
		zztlen = scale_expansion_zeroelim(ztlen, detzt, dez, detzzt);
		for (i = 0; i < zztlen; i++) {
			detzzt[i] *= 2.0;
		}
		ztztlen = scale_expansion_zeroelim(ztlen, detzt, deztail, detztzt);
		z1len = fast_expansion_sum_zeroelim(zzlen, detzz, zztlen, detzzt, z1);
		z2len = fast_expansion_sum_zeroelim(z1len, z1, ztztlen, detztzt, z2);
		xylen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, detxy);
		dlen = fast_expansion_sum_zeroelim(z2len, z2, xylen, detxy, ddet);

		ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
		cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
		deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, deter);

		return deter[deterlen - 1];
	}

	REAL insphereadapt(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe,
		REAL permanent)
	{
		INEXACT REAL aex, bex, cex, dex, aey, bey, cey, dey, aez, bez, cez, dez;
		REAL det, errbound;

		INEXACT REAL aexbey1, bexaey1, bexcey1, cexbey1;
		INEXACT REAL cexdey1, dexcey1, dexaey1, aexdey1;
		INEXACT REAL aexcey1, cexaey1, bexdey1, dexbey1;
		REAL aexbey0, bexaey0, bexcey0, cexbey0;
		REAL cexdey0, dexcey0, dexaey0, aexdey0;
		REAL aexcey0, cexaey0, bexdey0, dexbey0;
		REAL ab[4], bc[4], cd[4], da[4], ac[4], bd[4];
		INEXACT REAL ab3, bc3, cd3, da3, ac3, bd3;
		REAL abeps, bceps, cdeps, daeps, aceps, bdeps;
		REAL temp8a[8], temp8b[8], temp8c[8], temp16[16], temp24[24], temp48[48];
		int temp8alen, temp8blen, temp8clen, temp16len, temp24len, temp48len;
		REAL xdet[96], ydet[96], zdet[96], xydet[192];
		int xlen, ylen, zlen, xylen;
		REAL adet[288], bdet[288], cdet[288], ddet[288];
		int alen, blen, clen, dlen;
		REAL abdet[576], cddet[576];
		int ablen, cdlen;
		REAL fin1[1152];
		int finlength;

		REAL aextail, bextail, cextail, dextail;
		REAL aeytail, beytail, ceytail, deytail;
		REAL aeztail, beztail, ceztail, deztail;

		INEXACT REAL bvirt;
		REAL avirt, bround, around;
		INEXACT REAL c;
		INEXACT REAL abig;
		REAL ahi, alo, bhi, blo;
		REAL err1, err2, err3;
		INEXACT REAL _i, _j;
		REAL _0;


		aex = (REAL)(pa[0] - pe[0]);
		bex = (REAL)(pb[0] - pe[0]);
		cex = (REAL)(pc[0] - pe[0]);
		dex = (REAL)(pd[0] - pe[0]);
		aey = (REAL)(pa[1] - pe[1]);
		bey = (REAL)(pb[1] - pe[1]);
		cey = (REAL)(pc[1] - pe[1]);
		dey = (REAL)(pd[1] - pe[1]);
		aez = (REAL)(pa[2] - pe[2]);
		bez = (REAL)(pb[2] - pe[2]);
		cez = (REAL)(pc[2] - pe[2]);
		dez = (REAL)(pd[2] - pe[2]);

		Two_Product(aex, bey, aexbey1, aexbey0);
		Two_Product(bex, aey, bexaey1, bexaey0);
		Two_Two_Diff(aexbey1, aexbey0, bexaey1, bexaey0, ab3, ab[2], ab[1], ab[0]);
		ab[3] = ab3;

		Two_Product(bex, cey, bexcey1, bexcey0);
		Two_Product(cex, bey, cexbey1, cexbey0);
		Two_Two_Diff(bexcey1, bexcey0, cexbey1, cexbey0, bc3, bc[2], bc[1], bc[0]);
		bc[3] = bc3;

		Two_Product(cex, dey, cexdey1, cexdey0);
		Two_Product(dex, cey, dexcey1, dexcey0);
		Two_Two_Diff(cexdey1, cexdey0, dexcey1, dexcey0, cd3, cd[2], cd[1], cd[0]);
		cd[3] = cd3;

		Two_Product(dex, aey, dexaey1, dexaey0);
		Two_Product(aex, dey, aexdey1, aexdey0);
		Two_Two_Diff(dexaey1, dexaey0, aexdey1, aexdey0, da3, da[2], da[1], da[0]);
		da[3] = da3;

		Two_Product(aex, cey, aexcey1, aexcey0);
		Two_Product(cex, aey, cexaey1, cexaey0);
		Two_Two_Diff(aexcey1, aexcey0, cexaey1, cexaey0, ac3, ac[2], ac[1], ac[0]);
		ac[3] = ac3;

		Two_Product(bex, dey, bexdey1, bexdey0);
		Two_Product(dex, bey, dexbey1, dexbey0);
		Two_Two_Diff(bexdey1, bexdey0, dexbey1, dexbey0, bd3, bd[2], bd[1], bd[0]);
		bd[3] = bd3;

		temp8alen = scale_expansion_zeroelim(4, cd, bez, temp8a);
		temp8blen = scale_expansion_zeroelim(4, bd, -cez, temp8b);
		temp8clen = scale_expansion_zeroelim(4, bc, dez, temp8c);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a,
			temp8blen, temp8b, temp16);
		temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c,
			temp16len, temp16, temp24);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, aex, temp48);
		xlen = scale_expansion_zeroelim(temp48len, temp48, -aex, xdet);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, aey, temp48);
		ylen = scale_expansion_zeroelim(temp48len, temp48, -aey, ydet);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, aez, temp48);
		zlen = scale_expansion_zeroelim(temp48len, temp48, -aez, zdet);
		xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet);
		alen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, adet);

		temp8alen = scale_expansion_zeroelim(4, da, cez, temp8a);
		temp8blen = scale_expansion_zeroelim(4, ac, dez, temp8b);
		temp8clen = scale_expansion_zeroelim(4, cd, aez, temp8c);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a,
			temp8blen, temp8b, temp16);
		temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c,
			temp16len, temp16, temp24);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, bex, temp48);
		xlen = scale_expansion_zeroelim(temp48len, temp48, bex, xdet);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, bey, temp48);
		ylen = scale_expansion_zeroelim(temp48len, temp48, bey, ydet);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, bez, temp48);
		zlen = scale_expansion_zeroelim(temp48len, temp48, bez, zdet);
		xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet);
		blen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, bdet);

		temp8alen = scale_expansion_zeroelim(4, ab, dez, temp8a);
		temp8blen = scale_expansion_zeroelim(4, bd, aez, temp8b);
		temp8clen = scale_expansion_zeroelim(4, da, bez, temp8c);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a,
			temp8blen, temp8b, temp16);
		temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c,
			temp16len, temp16, temp24);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, cex, temp48);
		xlen = scale_expansion_zeroelim(temp48len, temp48, -cex, xdet);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, cey, temp48);
		ylen = scale_expansion_zeroelim(temp48len, temp48, -cey, ydet);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, cez, temp48);
		zlen = scale_expansion_zeroelim(temp48len, temp48, -cez, zdet);
		xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet);
		clen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, cdet);

		temp8alen = scale_expansion_zeroelim(4, bc, aez, temp8a);
		temp8blen = scale_expansion_zeroelim(4, ac, -bez, temp8b);
		temp8clen = scale_expansion_zeroelim(4, ab, cez, temp8c);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a,
			temp8blen, temp8b, temp16);
		temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c,
			temp16len, temp16, temp24);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, dex, temp48);
		xlen = scale_expansion_zeroelim(temp48len, temp48, dex, xdet);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, dey, temp48);
		ylen = scale_expansion_zeroelim(temp48len, temp48, dey, ydet);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, dez, temp48);
		zlen = scale_expansion_zeroelim(temp48len, temp48, dez, zdet);
		xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet);
		dlen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, ddet);

		ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
		cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
		finlength = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, fin1);

		det = estimate(finlength, fin1);
		errbound = isperrboundB * permanent;
		if ((det >= errbound) || (-det >= errbound)) {
			return det;
		}

		Two_Diff_Tail(pa[0], pe[0], aex, aextail);
		Two_Diff_Tail(pa[1], pe[1], aey, aeytail);
		Two_Diff_Tail(pa[2], pe[2], aez, aeztail);
		Two_Diff_Tail(pb[0], pe[0], bex, bextail);
		Two_Diff_Tail(pb[1], pe[1], bey, beytail);
		Two_Diff_Tail(pb[2], pe[2], bez, beztail);
		Two_Diff_Tail(pc[0], pe[0], cex, cextail);
		Two_Diff_Tail(pc[1], pe[1], cey, ceytail);
		Two_Diff_Tail(pc[2], pe[2], cez, ceztail);
		Two_Diff_Tail(pd[0], pe[0], dex, dextail);
		Two_Diff_Tail(pd[1], pe[1], dey, deytail);
		Two_Diff_Tail(pd[2], pe[2], dez, deztail);
		if ((aextail == 0.0) && (aeytail == 0.0) && (aeztail == 0.0)
			&& (bextail == 0.0) && (beytail == 0.0) && (beztail == 0.0)
			&& (cextail == 0.0) && (ceytail == 0.0) && (ceztail == 0.0)
			&& (dextail == 0.0) && (deytail == 0.0) && (deztail == 0.0)) {
			return det;
		}

		errbound = isperrboundC * permanent + resulterrbound * Absolute(det);
		abeps = (aex * beytail + bey * aextail)
			- (aey * bextail + bex * aeytail);
		bceps = (bex * ceytail + cey * bextail)
			- (bey * cextail + cex * beytail);
		cdeps = (cex * deytail + dey * cextail)
			- (cey * dextail + dex * ceytail);
		daeps = (dex * aeytail + aey * dextail)
			- (dey * aextail + aex * deytail);
		aceps = (aex * ceytail + cey * aextail)
			- (aey * cextail + cex * aeytail);
		bdeps = (bex * deytail + dey * bextail)
			- (bey * dextail + dex * beytail);
		det += (((bex * bex + bey * bey + bez * bez)
			* ((cez * daeps + dez * aceps + aez * cdeps)
			+ (ceztail * da3 + deztail * ac3 + aeztail * cd3))
			+ (dex * dex + dey * dey + dez * dez)
			* ((aez * bceps - bez * aceps + cez * abeps)
			+ (aeztail * bc3 - beztail * ac3 + ceztail * ab3)))
			- ((aex * aex + aey * aey + aez * aez)
			* ((bez * cdeps - cez * bdeps + dez * bceps)
			+ (beztail * cd3 - ceztail * bd3 + deztail * bc3))
			+ (cex * cex + cey * cey + cez * cez)
			* ((dez * abeps + aez * bdeps + bez * daeps)
			+ (deztail * ab3 + aeztail * bd3 + beztail * da3))))
			+ 2.0 * (((bex * bextail + bey * beytail + bez * beztail)
			* (cez * da3 + dez * ac3 + aez * cd3)
			+ (dex * dextail + dey * deytail + dez * deztail)
			* (aez * bc3 - bez * ac3 + cez * ab3))
			- ((aex * aextail + aey * aeytail + aez * aeztail)
			* (bez * cd3 - cez * bd3 + dez * bc3)
			+ (cex * cextail + cey * ceytail + cez * ceztail)
			* (dez * ab3 + aez * bd3 + bez * da3)));
		if ((det >= errbound) || (-det >= errbound)) {
			return det;
		}

		return insphereexact(pa, pb, pc, pd, pe);
	}

#ifdef USE_CGAL_PREDICATES

	REAL insphere(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe)
	{
		return (REAL)
			-cgal_pred_obj.side_of_oriented_sphere_3_object()
			(Point(pa[0], pa[1], pa[2]),
			Point(pb[0], pb[1], pb[2]),
			Point(pc[0], pc[1], pc[2]),
			Point(pd[0], pd[1], pd[2]),
			Point(pe[0], pe[1], pe[2]));
	}

#else

	REAL insphere(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe)
	{
		REAL aex, bex, cex, dex;
		REAL aey, bey, cey, dey;
		REAL aez, bez, cez, dez;
		REAL aexbey, bexaey, bexcey, cexbey, cexdey, dexcey, dexaey, aexdey;
		REAL aexcey, cexaey, bexdey, dexbey;
		REAL alift, blift, clift, dlift;
		REAL ab, bc, cd, da, ac, bd;
		REAL abc, bcd, cda, dab;
		REAL det;


		aex = pa[0] - pe[0];
		bex = pb[0] - pe[0];
		cex = pc[0] - pe[0];
		dex = pd[0] - pe[0];
		aey = pa[1] - pe[1];
		bey = pb[1] - pe[1];
		cey = pc[1] - pe[1];
		dey = pd[1] - pe[1];
		aez = pa[2] - pe[2];
		bez = pb[2] - pe[2];
		cez = pc[2] - pe[2];
		dez = pd[2] - pe[2];

		aexbey = aex * bey;
		bexaey = bex * aey;
		ab = aexbey - bexaey;
		bexcey = bex * cey;
		cexbey = cex * bey;
		bc = bexcey - cexbey;
		cexdey = cex * dey;
		dexcey = dex * cey;
		cd = cexdey - dexcey;
		dexaey = dex * aey;
		aexdey = aex * dey;
		da = dexaey - aexdey;

		aexcey = aex * cey;
		cexaey = cex * aey;
		ac = aexcey - cexaey;
		bexdey = bex * dey;
		dexbey = dex * bey;
		bd = bexdey - dexbey;

		abc = aez * bc - bez * ac + cez * ab;
		bcd = bez * cd - cez * bd + dez * bc;
		cda = cez * da + dez * ac + aez * cd;
		dab = dez * ab + aez * bd + bez * da;

		alift = aex * aex + aey * aey + aez * aez;
		blift = bex * bex + bey * bey + bez * bez;
		clift = cex * cex + cey * cey + cez * cez;
		dlift = dex * dex + dey * dey + dez * dez;

		det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);

		if (_use_inexact_arith) {
			return det;
		}

		if (_use_static_filter) {
			if (fabs(det) > ispstaticfilter) return det;
			//if (det > ispstaticfilter) return det;
			//if (det < minus_ispstaticfilter) return det;

		}

		REAL aezplus, bezplus, cezplus, dezplus;
		REAL aexbeyplus, bexaeyplus, bexceyplus, cexbeyplus;
		REAL cexdeyplus, dexceyplus, dexaeyplus, aexdeyplus;
		REAL aexceyplus, cexaeyplus, bexdeyplus, dexbeyplus;
		REAL permanent, errbound;

		aezplus = Absolute(aez);
		bezplus = Absolute(bez);
		cezplus = Absolute(cez);
		dezplus = Absolute(dez);
		aexbeyplus = Absolute(aexbey);
		bexaeyplus = Absolute(bexaey);
		bexceyplus = Absolute(bexcey);
		cexbeyplus = Absolute(cexbey);
		cexdeyplus = Absolute(cexdey);
		dexceyplus = Absolute(dexcey);
		dexaeyplus = Absolute(dexaey);
		aexdeyplus = Absolute(aexdey);
		aexceyplus = Absolute(aexcey);
		cexaeyplus = Absolute(cexaey);
		bexdeyplus = Absolute(bexdey);
		dexbeyplus = Absolute(dexbey);
		permanent = ((cexdeyplus + dexceyplus) * bezplus
			+ (dexbeyplus + bexdeyplus) * cezplus
			+ (bexceyplus + cexbeyplus) * dezplus)
			* alift
			+ ((dexaeyplus + aexdeyplus) * cezplus
			+ (aexceyplus + cexaeyplus) * dezplus
			+ (cexdeyplus + dexceyplus) * aezplus)
			* blift
			+ ((aexbeyplus + bexaeyplus) * dezplus
			+ (bexdeyplus + dexbeyplus) * aezplus
			+ (dexaeyplus + aexdeyplus) * bezplus)
			* clift
			+ ((bexceyplus + cexbeyplus) * aezplus
			+ (cexaeyplus + aexceyplus) * bezplus
			+ (aexbeyplus + bexaeyplus) * cezplus)
			* dlift;
		errbound = isperrboundA * permanent;
		if ((det > errbound) || (-det > errbound)) {
			return det;
		}

		return insphereadapt(pa, pb, pc, pd, pe, permanent);
	}

#endif // #ifdef USE_CGAL_PREDICATES

	/*****************************************************************************/
	/*                                                                           */
	/*  orient4d()   Return a positive value if the point pe lies above the      */
	/*               hyperplane passing through pa, pb, pc, and pd; "above" is   */
	/*               defined in a manner best found by trial-and-error.  Returns */
	/*               a negative value if pe lies below the hyperplane.  Returns  */
	/*               zero if the points are co-hyperplanar (not affinely         */
	/*               independent).  The result is also a rough approximation of  */
	/*               24 times the signed volume of the 4-simplex defined by the  */
	/*               five points.                                                */
	/*                                                                           */
	/*  Uses exact arithmetic if necessary to ensure a correct answer.  The      */
	/*  result returned is the determinant of a matrix.  This determinant is     */
	/*  computed adaptively, in the sense that exact arithmetic is used only to  */
	/*  the degree it is needed to ensure that the returned value has the        */
	/*  correct sign.  Hence, orient4d() is usually quite fast, but will run     */
	/*  more slowly when the input points are hyper-coplanar or nearly so.       */
	/*                                                                           */
	/*  See my Robust Predicates paper for details.                              */
	/*                                                                           */
	/*****************************************************************************/

	REAL orient4dexact(REAL* pa, REAL* pb, REAL* pc, REAL* pd, REAL* pe,
		REAL aheight, REAL bheight, REAL cheight, REAL dheight,
		REAL eheight)
	{
		INEXACT REAL axby1, bxcy1, cxdy1, dxey1, exay1;
		INEXACT REAL bxay1, cxby1, dxcy1, exdy1, axey1;
		INEXACT REAL axcy1, bxdy1, cxey1, dxay1, exby1;
		INEXACT REAL cxay1, dxby1, excy1, axdy1, bxey1;
		REAL axby0, bxcy0, cxdy0, dxey0, exay0;
		REAL bxay0, cxby0, dxcy0, exdy0, axey0;
		REAL axcy0, bxdy0, cxey0, dxay0, exby0;
		REAL cxay0, dxby0, excy0, axdy0, bxey0;
		REAL ab[4], bc[4], cd[4], de[4], ea[4];
		REAL ac[4], bd[4], ce[4], da[4], eb[4];
		REAL temp8a[8], temp8b[8], temp16[16];
		int temp8alen, temp8blen, temp16len;
		REAL abc[24], bcd[24], cde[24], dea[24], eab[24];
		REAL abd[24], bce[24], cda[24], deb[24], eac[24];
		int abclen, bcdlen, cdelen, dealen, eablen;
		int abdlen, bcelen, cdalen, deblen, eaclen;
		REAL temp48a[48], temp48b[48];
		int temp48alen, temp48blen;
		REAL abcd[96], bcde[96], cdea[96], deab[96], eabc[96];
		int abcdlen, bcdelen, cdealen, deablen, eabclen;
		REAL adet[192], bdet[192], cdet[192], ddet[192], edet[192];
		int alen, blen, clen, dlen, elen;
		REAL abdet[384], cddet[384], cdedet[576];
		int ablen, cdlen;
		REAL deter[960];
		int deterlen;
		int i;

		INEXACT REAL bvirt;
		REAL avirt, bround, around;
		INEXACT REAL c;
		INEXACT REAL abig;
		REAL ahi, alo, bhi, blo;
		REAL err1, err2, err3;
		INEXACT REAL _i, _j;
		REAL _0;


		Two_Product(pa[0], pb[1], axby1, axby0);
		Two_Product(pb[0], pa[1], bxay1, bxay0);
		Two_Two_Diff(axby1, axby0, bxay1, bxay0, ab[3], ab[2], ab[1], ab[0]);

		Two_Product(pb[0], pc[1], bxcy1, bxcy0);
		Two_Product(pc[0], pb[1], cxby1, cxby0);
		Two_Two_Diff(bxcy1, bxcy0, cxby1, cxby0, bc[3], bc[2], bc[1], bc[0]);

		Two_Product(pc[0], pd[1], cxdy1, cxdy0);
		Two_Product(pd[0], pc[1], dxcy1, dxcy0);
		Two_Two_Diff(cxdy1, cxdy0, dxcy1, dxcy0, cd[3], cd[2], cd[1], cd[0]);

		Two_Product(pd[0], pe[1], dxey1, dxey0);
		Two_Product(pe[0], pd[1], exdy1, exdy0);
		Two_Two_Diff(dxey1, dxey0, exdy1, exdy0, de[3], de[2], de[1], de[0]);

		Two_Product(pe[0], pa[1], exay1, exay0);
		Two_Product(pa[0], pe[1], axey1, axey0);
		Two_Two_Diff(exay1, exay0, axey1, axey0, ea[3], ea[2], ea[1], ea[0]);

		Two_Product(pa[0], pc[1], axcy1, axcy0);
		Two_Product(pc[0], pa[1], cxay1, cxay0);
		Two_Two_Diff(axcy1, axcy0, cxay1, cxay0, ac[3], ac[2], ac[1], ac[0]);

		Two_Product(pb[0], pd[1], bxdy1, bxdy0);
		Two_Product(pd[0], pb[1], dxby1, dxby0);
		Two_Two_Diff(bxdy1, bxdy0, dxby1, dxby0, bd[3], bd[2], bd[1], bd[0]);

		Two_Product(pc[0], pe[1], cxey1, cxey0);
		Two_Product(pe[0], pc[1], excy1, excy0);
		Two_Two_Diff(cxey1, cxey0, excy1, excy0, ce[3], ce[2], ce[1], ce[0]);

		Two_Product(pd[0], pa[1], dxay1, dxay0);
		Two_Product(pa[0], pd[1], axdy1, axdy0);
		Two_Two_Diff(dxay1, dxay0, axdy1, axdy0, da[3], da[2], da[1], da[0]);

		Two_Product(pe[0], pb[1], exby1, exby0);
		Two_Product(pb[0], pe[1], bxey1, bxey0);
		Two_Two_Diff(exby1, exby0, bxey1, bxey0, eb[3], eb[2], eb[1], eb[0]);

		temp8alen = scale_expansion_zeroelim(4, bc, pa[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, ac, -pb[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
			temp16);
		temp8alen = scale_expansion_zeroelim(4, ab, pc[2], temp8a);
		abclen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
			abc);

		temp8alen = scale_expansion_zeroelim(4, cd, pb[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, bd, -pc[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
			temp16);
		temp8alen = scale_expansion_zeroelim(4, bc, pd[2], temp8a);
		bcdlen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
			bcd);

		temp8alen = scale_expansion_zeroelim(4, de, pc[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, ce, -pd[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
			temp16);
		temp8alen = scale_expansion_zeroelim(4, cd, pe[2], temp8a);
		cdelen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
			cde);

		temp8alen = scale_expansion_zeroelim(4, ea, pd[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, da, -pe[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
			temp16);
		temp8alen = scale_expansion_zeroelim(4, de, pa[2], temp8a);
		dealen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
			dea);

		temp8alen = scale_expansion_zeroelim(4, ab, pe[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, eb, -pa[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
			temp16);
		temp8alen = scale_expansion_zeroelim(4, ea, pb[2], temp8a);
		eablen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
			eab);

		temp8alen = scale_expansion_zeroelim(4, bd, pa[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, da, pb[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
			temp16);
		temp8alen = scale_expansion_zeroelim(4, ab, pd[2], temp8a);
		abdlen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
			abd);

		temp8alen = scale_expansion_zeroelim(4, ce, pb[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, eb, pc[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
			temp16);
		temp8alen = scale_expansion_zeroelim(4, bc, pe[2], temp8a);
		bcelen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
			bce);

		temp8alen = scale_expansion_zeroelim(4, da, pc[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, ac, pd[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
			temp16);
		temp8alen = scale_expansion_zeroelim(4, cd, pa[2], temp8a);
		cdalen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
			cda);

		temp8alen = scale_expansion_zeroelim(4, eb, pd[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, bd, pe[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
			temp16);
		temp8alen = scale_expansion_zeroelim(4, de, pb[2], temp8a);
		deblen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
			deb);

		temp8alen = scale_expansion_zeroelim(4, ac, pe[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, ce, pa[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
			temp16);
		temp8alen = scale_expansion_zeroelim(4, ea, pc[2], temp8a);
		eaclen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
			eac);

		temp48alen = fast_expansion_sum_zeroelim(cdelen, cde, bcelen, bce, temp48a);
		temp48blen = fast_expansion_sum_zeroelim(deblen, deb, bcdlen, bcd, temp48b);
		for (i = 0; i < temp48blen; i++) {
			temp48b[i] = -temp48b[i];
		}
		bcdelen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
			temp48blen, temp48b, bcde);
		alen = scale_expansion_zeroelim(bcdelen, bcde, aheight, adet);

		temp48alen = fast_expansion_sum_zeroelim(dealen, dea, cdalen, cda, temp48a);
		temp48blen = fast_expansion_sum_zeroelim(eaclen, eac, cdelen, cde, temp48b);
		for (i = 0; i < temp48blen; i++) {
			temp48b[i] = -temp48b[i];
		}
		cdealen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
			temp48blen, temp48b, cdea);
		blen = scale_expansion_zeroelim(cdealen, cdea, bheight, bdet);

		temp48alen = fast_expansion_sum_zeroelim(eablen, eab, deblen, deb, temp48a);
		temp48blen = fast_expansion_sum_zeroelim(abdlen, abd, dealen, dea, temp48b);
		for (i = 0; i < temp48blen; i++) {
			temp48b[i] = -temp48b[i];
		}
		deablen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
			temp48blen, temp48b, deab);
		clen = scale_expansion_zeroelim(deablen, deab, cheight, cdet);

		temp48alen = fast_expansion_sum_zeroelim(abclen, abc, eaclen, eac, temp48a);
		temp48blen = fast_expansion_sum_zeroelim(bcelen, bce, eablen, eab, temp48b);
		for (i = 0; i < temp48blen; i++) {
			temp48b[i] = -temp48b[i];
		}
		eabclen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
			temp48blen, temp48b, eabc);
		dlen = scale_expansion_zeroelim(eabclen, eabc, dheight, ddet);

		temp48alen = fast_expansion_sum_zeroelim(bcdlen, bcd, abdlen, abd, temp48a);
		temp48blen = fast_expansion_sum_zeroelim(cdalen, cda, abclen, abc, temp48b);
		for (i = 0; i < temp48blen; i++) {
			temp48b[i] = -temp48b[i];
		}
		abcdlen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
			temp48blen, temp48b, abcd);
		elen = scale_expansion_zeroelim(abcdlen, abcd, eheight, edet);

		ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
		cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
		cdelen = fast_expansion_sum_zeroelim(cdlen, cddet, elen, edet, cdedet);
		deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdelen, cdedet, deter);

		return deter[deterlen - 1];
	}

	REAL orient4dadapt(REAL* pa, REAL* pb, REAL* pc, REAL* pd, REAL* pe,
		REAL aheight, REAL bheight, REAL cheight, REAL dheight,
		REAL eheight, REAL permanent)
	{
		INEXACT REAL aex, bex, cex, dex, aey, bey, cey, dey, aez, bez, cez, dez;
		INEXACT REAL aeheight, beheight, ceheight, deheight;
		REAL det, errbound;

		INEXACT REAL aexbey1, bexaey1, bexcey1, cexbey1;
		INEXACT REAL cexdey1, dexcey1, dexaey1, aexdey1;
		INEXACT REAL aexcey1, cexaey1, bexdey1, dexbey1;
		REAL aexbey0, bexaey0, bexcey0, cexbey0;
		REAL cexdey0, dexcey0, dexaey0, aexdey0;
		REAL aexcey0, cexaey0, bexdey0, dexbey0;
		REAL ab[4], bc[4], cd[4], da[4], ac[4], bd[4];
		INEXACT REAL ab3, bc3, cd3, da3, ac3, bd3;
		REAL abeps, bceps, cdeps, daeps, aceps, bdeps;
		REAL temp8a[8], temp8b[8], temp8c[8], temp16[16], temp24[24];
		int temp8alen, temp8blen, temp8clen, temp16len, temp24len;
		REAL adet[48], bdet[48], cdet[48], ddet[48];
		int alen, blen, clen, dlen;
		REAL abdet[96], cddet[96];
		int ablen, cdlen;
		REAL fin1[192];
		int finlength;

		REAL aextail, bextail, cextail, dextail;
		REAL aeytail, beytail, ceytail, deytail;
		REAL aeztail, beztail, ceztail, deztail;
		REAL aeheighttail, beheighttail, ceheighttail, deheighttail;

		INEXACT REAL bvirt;
		REAL avirt, bround, around;
		INEXACT REAL c;
		INEXACT REAL abig;
		REAL ahi, alo, bhi, blo;
		REAL err1, err2, err3;
		INEXACT REAL _i, _j;
		REAL _0;


		aex = (REAL)(pa[0] - pe[0]);
		bex = (REAL)(pb[0] - pe[0]);
		cex = (REAL)(pc[0] - pe[0]);
		dex = (REAL)(pd[0] - pe[0]);
		aey = (REAL)(pa[1] - pe[1]);
		bey = (REAL)(pb[1] - pe[1]);
		cey = (REAL)(pc[1] - pe[1]);
		dey = (REAL)(pd[1] - pe[1]);
		aez = (REAL)(pa[2] - pe[2]);
		bez = (REAL)(pb[2] - pe[2]);
		cez = (REAL)(pc[2] - pe[2]);
		dez = (REAL)(pd[2] - pe[2]);
		aeheight = (REAL)(aheight - eheight);
		beheight = (REAL)(bheight - eheight);
		ceheight = (REAL)(cheight - eheight);
		deheight = (REAL)(dheight - eheight);

		Two_Product(aex, bey, aexbey1, aexbey0);
		Two_Product(bex, aey, bexaey1, bexaey0);
		Two_Two_Diff(aexbey1, aexbey0, bexaey1, bexaey0, ab3, ab[2], ab[1], ab[0]);
		ab[3] = ab3;

		Two_Product(bex, cey, bexcey1, bexcey0);
		Two_Product(cex, bey, cexbey1, cexbey0);
		Two_Two_Diff(bexcey1, bexcey0, cexbey1, cexbey0, bc3, bc[2], bc[1], bc[0]);
		bc[3] = bc3;

		Two_Product(cex, dey, cexdey1, cexdey0);
		Two_Product(dex, cey, dexcey1, dexcey0);
		Two_Two_Diff(cexdey1, cexdey0, dexcey1, dexcey0, cd3, cd[2], cd[1], cd[0]);
		cd[3] = cd3;

		Two_Product(dex, aey, dexaey1, dexaey0);
		Two_Product(aex, dey, aexdey1, aexdey0);
		Two_Two_Diff(dexaey1, dexaey0, aexdey1, aexdey0, da3, da[2], da[1], da[0]);
		da[3] = da3;

		Two_Product(aex, cey, aexcey1, aexcey0);
		Two_Product(cex, aey, cexaey1, cexaey0);
		Two_Two_Diff(aexcey1, aexcey0, cexaey1, cexaey0, ac3, ac[2], ac[1], ac[0]);
		ac[3] = ac3;

		Two_Product(bex, dey, bexdey1, bexdey0);
		Two_Product(dex, bey, dexbey1, dexbey0);
		Two_Two_Diff(bexdey1, bexdey0, dexbey1, dexbey0, bd3, bd[2], bd[1], bd[0]);
		bd[3] = bd3;

		temp8alen = scale_expansion_zeroelim(4, cd, bez, temp8a);
		temp8blen = scale_expansion_zeroelim(4, bd, -cez, temp8b);
		temp8clen = scale_expansion_zeroelim(4, bc, dez, temp8c);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a,
			temp8blen, temp8b, temp16);
		temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c,
			temp16len, temp16, temp24);
		alen = scale_expansion_zeroelim(temp24len, temp24, -aeheight, adet);

		temp8alen = scale_expansion_zeroelim(4, da, cez, temp8a);
		temp8blen = scale_expansion_zeroelim(4, ac, dez, temp8b);
		temp8clen = scale_expansion_zeroelim(4, cd, aez, temp8c);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a,
			temp8blen, temp8b, temp16);
		temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c,
			temp16len, temp16, temp24);
		blen = scale_expansion_zeroelim(temp24len, temp24, beheight, bdet);

		temp8alen = scale_expansion_zeroelim(4, ab, dez, temp8a);
		temp8blen = scale_expansion_zeroelim(4, bd, aez, temp8b);
		temp8clen = scale_expansion_zeroelim(4, da, bez, temp8c);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a,
			temp8blen, temp8b, temp16);
		temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c,
			temp16len, temp16, temp24);
		clen = scale_expansion_zeroelim(temp24len, temp24, -ceheight, cdet);

		temp8alen = scale_expansion_zeroelim(4, bc, aez, temp8a);
		temp8blen = scale_expansion_zeroelim(4, ac, -bez, temp8b);
		temp8clen = scale_expansion_zeroelim(4, ab, cez, temp8c);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a,
			temp8blen, temp8b, temp16);
		temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c,
			temp16len, temp16, temp24);
		dlen = scale_expansion_zeroelim(temp24len, temp24, deheight, ddet);

		ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
		cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
		finlength = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, fin1);

		det = estimate(finlength, fin1);
		errbound = isperrboundB * permanent;
		if ((det >= errbound) || (-det >= errbound)) {
			return det;
		}

		Two_Diff_Tail(pa[0], pe[0], aex, aextail);
		Two_Diff_Tail(pa[1], pe[1], aey, aeytail);
		Two_Diff_Tail(pa[2], pe[2], aez, aeztail);
		Two_Diff_Tail(aheight, eheight, aeheight, aeheighttail);
		Two_Diff_Tail(pb[0], pe[0], bex, bextail);
		Two_Diff_Tail(pb[1], pe[1], bey, beytail);
		Two_Diff_Tail(pb[2], pe[2], bez, beztail);
		Two_Diff_Tail(bheight, eheight, beheight, beheighttail);
		Two_Diff_Tail(pc[0], pe[0], cex, cextail);
		Two_Diff_Tail(pc[1], pe[1], cey, ceytail);
		Two_Diff_Tail(pc[2], pe[2], cez, ceztail);
		Two_Diff_Tail(cheight, eheight, ceheight, ceheighttail);
		Two_Diff_Tail(pd[0], pe[0], dex, dextail);
		Two_Diff_Tail(pd[1], pe[1], dey, deytail);
		Two_Diff_Tail(pd[2], pe[2], dez, deztail);
		Two_Diff_Tail(dheight, eheight, deheight, deheighttail);
		if ((aextail == 0.0) && (aeytail == 0.0) && (aeztail == 0.0)
			&& (bextail == 0.0) && (beytail == 0.0) && (beztail == 0.0)
			&& (cextail == 0.0) && (ceytail == 0.0) && (ceztail == 0.0)
			&& (dextail == 0.0) && (deytail == 0.0) && (deztail == 0.0)
			&& (aeheighttail == 0.0) && (beheighttail == 0.0)
			&& (ceheighttail == 0.0) && (deheighttail == 0.0)) {
			return det;
		}

		errbound = isperrboundC * permanent + resulterrbound * Absolute(det);
		abeps = (aex * beytail + bey * aextail)
			- (aey * bextail + bex * aeytail);
		bceps = (bex * ceytail + cey * bextail)
			- (bey * cextail + cex * beytail);
		cdeps = (cex * deytail + dey * cextail)
			- (cey * dextail + dex * ceytail);
		daeps = (dex * aeytail + aey * dextail)
			- (dey * aextail + aex * deytail);
		aceps = (aex * ceytail + cey * aextail)
			- (aey * cextail + cex * aeytail);
		bdeps = (bex * deytail + dey * bextail)
			- (bey * dextail + dex * beytail);
		det += ((beheight
			* ((cez * daeps + dez * aceps + aez * cdeps)
			+ (ceztail * da3 + deztail * ac3 + aeztail * cd3))
			+ deheight
			* ((aez * bceps - bez * aceps + cez * abeps)
			+ (aeztail * bc3 - beztail * ac3 + ceztail * ab3)))
			- (aeheight
			* ((bez * cdeps - cez * bdeps + dez * bceps)
			+ (beztail * cd3 - ceztail * bd3 + deztail * bc3))
			+ ceheight
			* ((dez * abeps + aez * bdeps + bez * daeps)
			+ (deztail * ab3 + aeztail * bd3 + beztail * da3))))
			+ ((beheighttail * (cez * da3 + dez * ac3 + aez * cd3)
			+ deheighttail * (aez * bc3 - bez * ac3 + cez * ab3))
			- (aeheighttail * (bez * cd3 - cez * bd3 + dez * bc3)
			+ ceheighttail * (dez * ab3 + aez * bd3 + bez * da3)));
		if ((det >= errbound) || (-det >= errbound)) {
			return det;
		}

		return orient4dexact(pa, pb, pc, pd, pe,
			aheight, bheight, cheight, dheight, eheight);
	}

	REAL orient4d(REAL* pa, REAL* pb, REAL* pc, REAL* pd, REAL* pe,
		REAL aheight, REAL bheight, REAL cheight, REAL dheight,
		REAL eheight)
	{
		REAL aex, bex, cex, dex;
		REAL aey, bey, cey, dey;
		REAL aez, bez, cez, dez;
		REAL aexbey, bexaey, bexcey, cexbey, cexdey, dexcey, dexaey, aexdey;
		REAL aexcey, cexaey, bexdey, dexbey;
		REAL aeheight, beheight, ceheight, deheight;
		REAL ab, bc, cd, da, ac, bd;
		REAL abc, bcd, cda, dab;
		REAL aezplus, bezplus, cezplus, dezplus;
		REAL aexbeyplus, bexaeyplus, bexceyplus, cexbeyplus;
		REAL cexdeyplus, dexceyplus, dexaeyplus, aexdeyplus;
		REAL aexceyplus, cexaeyplus, bexdeyplus, dexbeyplus;
		REAL det;
		REAL permanent, errbound;


		aex = pa[0] - pe[0];
		bex = pb[0] - pe[0];
		cex = pc[0] - pe[0];
		dex = pd[0] - pe[0];
		aey = pa[1] - pe[1];
		bey = pb[1] - pe[1];
		cey = pc[1] - pe[1];
		dey = pd[1] - pe[1];
		aez = pa[2] - pe[2];
		bez = pb[2] - pe[2];
		cez = pc[2] - pe[2];
		dez = pd[2] - pe[2];
		aeheight = aheight - eheight;
		beheight = bheight - eheight;
		ceheight = cheight - eheight;
		deheight = dheight - eheight;

		aexbey = aex * bey;
		bexaey = bex * aey;
		ab = aexbey - bexaey;
		bexcey = bex * cey;
		cexbey = cex * bey;
		bc = bexcey - cexbey;
		cexdey = cex * dey;
		dexcey = dex * cey;
		cd = cexdey - dexcey;
		dexaey = dex * aey;
		aexdey = aex * dey;
		da = dexaey - aexdey;

		aexcey = aex * cey;
		cexaey = cex * aey;
		ac = aexcey - cexaey;
		bexdey = bex * dey;
		dexbey = dex * bey;
		bd = bexdey - dexbey;

		abc = aez * bc - bez * ac + cez * ab;
		bcd = bez * cd - cez * bd + dez * bc;
		cda = cez * da + dez * ac + aez * cd;
		dab = dez * ab + aez * bd + bez * da;

		det = (deheight * abc - ceheight * dab) + (beheight * cda - aeheight * bcd);

		aezplus = Absolute(aez);
		bezplus = Absolute(bez);
		cezplus = Absolute(cez);
		dezplus = Absolute(dez);
		aexbeyplus = Absolute(aexbey);
		bexaeyplus = Absolute(bexaey);
		bexceyplus = Absolute(bexcey);
		cexbeyplus = Absolute(cexbey);
		cexdeyplus = Absolute(cexdey);
		dexceyplus = Absolute(dexcey);
		dexaeyplus = Absolute(dexaey);
		aexdeyplus = Absolute(aexdey);
		aexceyplus = Absolute(aexcey);
		cexaeyplus = Absolute(cexaey);
		bexdeyplus = Absolute(bexdey);
		dexbeyplus = Absolute(dexbey);
		permanent = ((cexdeyplus + dexceyplus) * bezplus
			+ (dexbeyplus + bexdeyplus) * cezplus
			+ (bexceyplus + cexbeyplus) * dezplus)
			* Absolute(aeheight)
			+ ((dexaeyplus + aexdeyplus) * cezplus
			+ (aexceyplus + cexaeyplus) * dezplus
			+ (cexdeyplus + dexceyplus) * aezplus)
			* Absolute(beheight)
			+ ((aexbeyplus + bexaeyplus) * dezplus
			+ (bexdeyplus + dexbeyplus) * aezplus
			+ (dexaeyplus + aexdeyplus) * bezplus)
			* Absolute(ceheight)
			+ ((bexceyplus + cexbeyplus) * aezplus
			+ (cexaeyplus + aexceyplus) * bezplus
			+ (aexbeyplus + bexaeyplus) * cezplus)
			* Absolute(deheight);
		errbound = isperrboundA * permanent;
		if ((det > errbound) || (-det > errbound)) {
			return det;
		}

		return orient4dadapt(pa, pb, pc, pd, pe,
			aheight, bheight, cheight, dheight, eheight, permanent);
	}




}

namespace ContainerFunctions {

/** Move elements in the range [b,e) to front of the container if Func returns true
* This function is especially efficient if we have little items which need to move to the front
* This function respects the order of the elements
* @param b iterator to the first item
* @param e iterator to the last item
* @return Iterator r  where the range [r,e] is the back part of the vector where Func returned true
*/
template<typename Iterator, typename Func>
Iterator moveElementsToFrontIf(Iterator b, Iterator  e, Func f) {

    Iterator w = b;  // write pointer
    while  ( b != e ) {
        if( !f(*b) ) { // if not move to front increment read pointer
            b++;
            continue;
        }

        if(w!=b) {   // copy only if not at same position!
            *w = *b; // copy  value to front (test is not true)
        }
        ++w;
        ++b;
    }
    return w;
}

/** Move elements in the range [b,e] to the back of the container if Func returns true
* This function is especially efficient if we have little items which need to move to the front
* Caution: This function does not respect the order of the items (so do not use if container should stay sorted)!
* @param b begin iterator
* @param e end iterator (no past the end iterator of a container!)
* @return Iterator r  where the range [r,e] is the back part of the vector where Func returned true
*/
template<typename Iterator,typename Func>
Iterator moveElementsToBackIf(Iterator b, Iterator  e, Func f) {

    if( b == e  ) {
        return b;
    }
    while(b != e){

        if(f(*b)){
            // Move end to next swappable item
            while(f(*e)){
                --e;
                if( e == b){
                    return e;
                };
            }

            // swap with back
            if(b!=e){
                std::swap(*b,*e);
            }
        }
        ++b;
    }
    return f(*e)? e : ++e;
}


/** Move all sequences S in [b,e) to the front where each element i
* of the subsequence S =[start,end] fullfils c(start,i).
* Example: 1 2 2 2 3 3 3 4 5 5 6  -> 1 2 3 4 5 6 (with Comp c=AlmostEqual)
* This function is especially efficient if we have little items which need to move to the front
* This function respects the order of the elements
* @return Iterator r  where the range [r,e] is the back part of the vector where  Comp c returned true
*/
template<typename Iterator, typename Comp>
Iterator moveConsecutiveToFrontIf(Iterator b, Iterator  e, Comp c) {


    if( std::distance(b,e)<2 ) {
        return e;
    }

    Iterator comp = b;
    Iterator write = ++b;
    while  ( b != e ) {
        if( !c(*comp, *b ) ) { // if we can skip element or if
            ++b;
            continue;
        }

        if(write!=b) {   // copy only if not at same position!
            *write = *b; // copy  value to front (test is not true)
        }
        //std::cout << *dest << std::endl;
        comp = write++;
        ++b;
    }

    return write;
}
/** Move all sequences S in [b,e) to the front where each element i
* of the subsequence S =[start,end] fullfils Func(start,i).
* This function skips elements for which s(i) is true
* Example: 1 2 2 2 3 3 3 4 5 5 6  -> 1 2 3 4 5 6 (with Func=AlmostEqual)
* This function is especially efficient if we have little items which need to move to the front
* This function respects the order of the elements
* @return Iterator r  where the range [r,e] is the back part of the vector where Func returned true
*/
template<typename Iterator, typename Func, typename Skip>
Iterator moveConsecutiveToFrontIf(Iterator b, Iterator  e, Func f, Skip s) {


    if( std::distance(b,e)<2 ) {
        return e;
    }
    Iterator comp;
    Iterator write = b;
    // skip all at beginning
    while(b!=e && s(*b)){
        ++b;
    }
    if(b!=e){
        // write first element to front
        if(write!=b){
            *write = *b;
        }
        comp = write++; // shift write pointer (comp is previous)

        // Start loop over all elements
        while  ( b != e ) {
            if(s(*b) || !f(*comp, *b ) ) {
                ++b;
                continue;
            }

            if(write!=b) {   // copy only if not at same position!
                *write = *b; // copy  value to front (test is not true)
            }
            //std::cout << *dest << std::endl;
            comp = write++;
            ++b;
        }
    }

    return write;
}



};

#define ApproxMVBB_ERRORMSG(x)
#define ApproxMVBB_ASSERTMSG(x)
#define ApproxMVBB_STATIC_ASSERT(x)

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

namespace AngleFunctions {

/** Maps angle to the range [-pi,pi]
    * Sawtooth function as g = f(x+a/2) - a/2 , f = x - floor(x/a)*a , a = 2*pi
    */
inline double mapToPi(double x) {
    return   x -  std::floor( x/(2.0*M_PI) + 0.5 ) * M_PI * 2.0;
}

/** Maps angle to the range [0,2*pi]
    * Sawtooth function as f = x - floor(x/a)*a , a = 2*pi
    */
inline double mapTo2Pi(double x) {
    return   x -  std::floor( x/(2.0*M_PI)) * M_PI * 2.0;
}

/** Relative angle measured from angle */
inline double relativeAnglePi(double angle, double angle2) {
    return mapToPi(angle2-angle);
}
/** Relative angle measured from angle */
inline double relativeAngle2Pi(double angle, double angle2) {
    return mapTo2Pi(angle2-angle);
}

};


namespace MathFunctions {
/** Greates common divisor */
template<bool argsPositive = false,typename T>
typename std::enable_if<std::is_integral<T>::value,T >::type
gcd2( T  a, T  b ) {
    if(!argsPositive) {
        a = std::abs( a );
        b = std::abs( b );
    }
    if( a == 0  ||  a == b ) {
        return  b;
    }
    if( b == 0 ) {
        return  a;
    }
    if( a > b ) {
        return  gcd2<true,T>( a % b, b );
    } else {
        return  gcd2<true,T>( a, b % a );
    }
}
/** Greates common divisor */
template<bool argsPositive = false, typename T>
typename std::enable_if<std::is_integral<T>::value,T>::type
gcd3( T  a, T  b, T  c ) {
    if(!argsPositive) {
        a = std::abs( a );
        b = std::abs( b );
        c = std::abs( c );
    }
    if   ( a == 0 )
        return  gcd2<true,T>( b, c );
    if   ( b == 0 )
        return  gcd2<true,T>( a, c );
    if   ( c == 0 )
        return  gcd2<true,T>( a, b );

    return  gcd2<true,T>( a, gcd2<true,T>( b, c ) );
}
}

namespace PointFunctions {
template<typename Derived>
void applyRandomRotTrans(Eigen::MatrixBase<Derived> & points) {

    EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Derived,3, Eigen::Dynamic)
            Quaternion q;
    q.coeffs().setRandom();
    q.normalize();
    Matrix33 R = q.matrix();
    Vector3 trans;
    trans.setRandom();
    points = R*points;
    points.colwise() += trans;
    std::cout << "Applied Transformation" << std::endl;
}


template<typename VecT1, typename VecT2>
inline bool almostEqualAbs(const VecT1  & a, const VecT2  & b, double eps = 1.0e-8 ) {
    return ((a-b).array().abs() <= eps).all();
}


template <size_t size>
class TypeWithSize {
public:
    // This prevents the user from using TypeWithSize<N> with incorrect
    // values of N.
    typedef void UInt;
};

// The specialization for size 4.
template <>
class TypeWithSize<4> {
public:
    // unsigned int has size 4 in both gcc and MSVC.
    //
    // As base/basictypes.h doesn't compile on Windows, we cannot use
    // uint32, uint64, and etc here.
    typedef int Int;
    typedef unsigned int UInt;
};

// The specialization for size 8.
template <>
class TypeWithSize<8> {
public:
#if GTEST_OS_WINDOWS
    typedef __int64 Int;
    typedef unsigned __int64 UInt;
#else
    typedef long long Int;  // NOLINT
    typedef unsigned long long UInt;  // NOLINT
#endif  // GTEST_OS_WINDOWS
};

template <typename RawType>
class FloatingPoint {
public:
    // Defines the unsigned integer type that has the same size as the
    // floating point number.
    typedef typename TypeWithSize<sizeof(RawType)>::UInt Bits;

    // Constants.

    // # of bits in a number.
    static const size_t kBitCount = 8*sizeof(RawType);

    // # of fraction bits in a number.
    static const size_t kFractionBitCount =
            std::numeric_limits<RawType>::digits - 1;

    // # of exponent bits in a number.
    static const size_t kExponentBitCount = kBitCount - 1 - kFractionBitCount;

    // The mask for the sign bit.
    static const Bits kSignBitMask = static_cast<Bits>(1) << (kBitCount - 1);

    // The mask for the fraction bits.
    static const Bits kFractionBitMask =
            ~static_cast<Bits>(0) >> (kExponentBitCount + 1);

    // The mask for the exponent bits.
    static const Bits kExponentBitMask = ~(kSignBitMask | kFractionBitMask);

    // How many ULP's (Units in the Last Place) we want to tolerate when
    // comparing two numbers.  The larger the value, the more error we
    // allow.  A 0 value means that two numbers must be exactly the same
    // to be considered equal.
    //
    // The maximum error of a single floating-point operation is 0.5
    // units in the last place.  On Intel CPU's, all floating-point
    // calculations are done with 80-bit doubleision, while double has 64
    // bits.  Therefore, 4 should be enough for ordinary use.
    //
    // See the following article for more details on ULP:
    // http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
    static const size_t kMaxUlps = 4;

    // Constructs a FloatingPoint from a raw floating-point number.
    //
    // On an Intel CPU, passing a non-normalized NAN (Not a Number)
    // around may change its bits, although the new value is guaranteed
    // to be also a NAN.  Therefore, don't expect this constructor to
    // preserve the bits in x when x is a NAN.
    explicit FloatingPoint(const RawType& x) { u_.value_ = x; }

    // Static methods

    // Reinterprets a bit pattern as a floating-point number.
    //
    // This function is needed to test the AlmostEquals() method.
    static RawType ReinterpretBits(const Bits bits) {
        FloatingPoint fp(0);
        fp.u_.bits_ = bits;
        return fp.u_.value_;
    }

    // Returns the floating-point number that represent positive infinity.
    static RawType Infinity() {
        return ReinterpretBits(kExponentBitMask);
    }

    // Returns the maximum representable finite floating-point number.
    static RawType Max();

    // Non-static methods

    // Returns the bits that represents this number.
    const Bits &bits() const { return u_.bits_; }

    // Returns the exponent bits of this number.
    Bits exponent_bits() const { return kExponentBitMask & u_.bits_; }

    // Returns the fraction bits of this number.
    Bits fraction_bits() const { return kFractionBitMask & u_.bits_; }

    // Returns the sign bit of this number.
    Bits sign_bit() const { return kSignBitMask & u_.bits_; }

    // Returns true iff this is NAN (not a number).
    bool is_nan() const {
        // It's a NAN if the exponent bits are all ones and the fraction
        // bits are not entirely zeros.
        return (exponent_bits() == kExponentBitMask) && (fraction_bits() != 0);
    }

    // Returns true iff this number is at most kMaxUlps ULP's away from
    // rhs.  In particular, this function:
    //
    //   - returns false if either number is (or both are) NAN.
    //   - treats really large numbers as almost equal to infinity.
    //   - thinks +0.0 and -0.0 are 0 DLP's apart.
    bool AlmostEquals(const FloatingPoint& rhs) const {
        // The IEEE standard says that any comparison operation involving
        // a NAN must return false.
        if (is_nan() || rhs.is_nan()) return false;

        return DistanceBetweenSignAndMagnitudeNumbers(u_.bits_, rhs.u_.bits_)
                <= kMaxUlps;
    }

private:
    // The data type used to store the actual floating-point number.
    union FloatingPointUnion {
        RawType value_;  // The raw floating-point number.
        Bits bits_;      // The bits that represent the number.
    };

    // Converts an integer from the sign-and-magnitude representation to
    // the biased representation.  More doubleisely, let N be 2 to the
    // power of (kBitCount - 1), an integer x is represented by the
    // unsigned number x + N.
    //
    // For instance,
    //
    //   -N + 1 (the most negative number representable using
    //          sign-and-magnitude) is represented by 1;
    //   0      is represented by N; and
    //   N - 1  (the biggest number representable using
    //          sign-and-magnitude) is represented by 2N - 1.
    //
    // Read http://en.wikipedia.org/wiki/Signed_number_representations
    // for more details on signed number representations.
    static Bits SignAndMagnitudeToBiased(const Bits &sam) {
        if (kSignBitMask & sam) {
            // sam represents a negative number.
            return ~sam + 1;
        } else {
            // sam represents a positive number.
            return kSignBitMask | sam;
        }
    }

    // Given two numbers in the sign-and-magnitude representation,
    // returns the distance between them as an unsigned number.
    static Bits DistanceBetweenSignAndMagnitudeNumbers(const Bits &sam1,
                                                       const Bits &sam2) {
        const Bits biased1 = SignAndMagnitudeToBiased(sam1);
        const Bits biased2 = SignAndMagnitudeToBiased(sam2);
        return (biased1 >= biased2) ? (biased1 - biased2) : (biased2 - biased1);
    }

    FloatingPointUnion u_;
};

template<typename VecT1, typename VecT2>
inline bool almostEqualUlp(const VecT1  & a, const VecT2  & b, double eps = 1.0e-8 ) {
    bool ret = true;
    for(unsigned int i=0;i<a.size();i++){
        ret = ret && FloatingPoint<double>(a(i)).AlmostEquals(FloatingPoint<double>(b(i)));
    }
    return ret;
}

template<typename VecT1, typename VecT2>
inline bool equal(const VecT1  & a, const VecT2  & b) {
    return (a.array() == b.array()).all();
}

/** vec1 = b-a and vec2 = c-a */
template<typename VecT1, typename VecT2, typename VecT3>
inline int orient2d(const VecT1  & a, const VecT2  & b, const VecT3  & c) {
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(VecT1,2)
            EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(VecT2,2)
            EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(VecT3,2)

            double f_A = GeometryPredicates::orient2d(const_cast<double*>(a.data()),
                                                      const_cast<double*>(b.data()),
                                                      const_cast<double*>(c.data()));

    return  ( ( f_A < 0.0 )?   -1   :   ( (f_A > 0.0)? 1 : 0) );

}

/** vec1 = b-a and vec2 = c-a */
template<typename VecT1, typename VecT2, typename VecT3>
inline bool leftTurn(const VecT1  & a, const VecT2  & b, const VecT3  & c) {
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(VecT1,2)
            EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(VecT2,2)
            EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(VecT3,2)
            return orient2d(a,b,c) > 0;
}

/** Postcondition: Points need to be collinear */
template<typename VecT1, typename VecT2, typename VecT3>
inline int collinearAreOrderedAlongLine(const VecT1  & a, const VecT2  & b, const VecT3  & c) {
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(VecT1,2)
            EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(VecT2,2)
            EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(VecT3,2)

            if (a(0) < b(0)) return !(c(0) < b(0));
    if (b(0) < a(0)) return !(b(0) < c(0));
    if (a(1) < b(1)) return !(c(1) < b(1));
    if (b(1) < a(1)) return !(b(1) < c(1));
    return true; // a==b

}



/** Get angle measures from x-Axis through point a */
template<typename VecT1, typename VecT2>
inline double getAngle(const VecT1  & a, const VecT2 & b) {
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(VecT1,2)
            EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(VecT2,2)
            Vector2 t = b - a ;
    double angle = std::atan2(t(1),t(0));
    if(angle<0.0) {
        angle += 2*M_PI;
    }
    return angle;
}

template<typename VecT1, typename VecT2>
Vector2 intersectLines(const VecT1 & p1, double ang1,
                       const VecT2 & p2, double ang2, double eps = 1e-10) {
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(VecT1,2)
            EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(VecT2,2)
            using namespace std;
    // Two lines p1 + a*t1 =  p2 + b*t2;
    // [-a, b]^-1 * (p1-p2) = [t1;t2]
    // p = (p1-p2)
    // =>
    // t1 = px by - bx py / (ay bx - ax by);

    Vector2 a;
    a(0) = cos(ang1);
    a(1) = sin(ang1);

    double bx = cos(ang2);
    double by = sin(ang2);

    double nom = (p1(0)-p2(0))*by  -  (p1(1)-p2(1)) *bx;

    double det = a(1)*bx   -  a(0)*by;

    if( det == 0.0 ) { // lines are collinear
        if(abs(nom) < eps ) { // if the two lines are almost identical! rot( p1-p2, 90 degrees) almost orthogonal to b
            return p1;
        }
        a(0) = std::numeric_limits<double>::infinity();
        a(1) = std::numeric_limits<double>::infinity();
        return a;
    }

    return  (nom / det) *a + p1;
}

/** Traverse first y-Axis then if equal check x-Axis */
template<typename Derived>
inline unsigned int minPointYX(const Eigen::MatrixBase<Derived> & points) {
    EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Derived,2, Eigen::Dynamic)
            unsigned int index = 0;
    for(unsigned int i=1; i<points.cols(); ++i) {
        if( points(1,i) < points(1,index) ) {
            index = i;
        } else if( points(1,i) == points(1,index)  && points(0,i) <  points(0,index) ) {
            index = i;
        }
    }
    return index;
}

/** Traverse first y-Axis then if equal check x-Axis */
template<typename Derived>
inline unsigned int minPointXY(const Eigen::MatrixBase<Derived> & points) {
    EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Derived,2, Eigen::Dynamic)
            unsigned int index = 0;
    for(unsigned int i=1; i<points.cols(); ++i) {
        if( points(0,i) < points(0,index) ) {
            index = i;
        } else if( points(0,i) == points(0,index)  && points(1,i) <  points(1,index) ) {
            index = i;
        }
    }
    return index;
}

template<unsigned int Dimension,typename TVector,typename Derived>
std::pair<TVector,TVector> estimateDiameter(const Eigen::MatrixBase<Derived> & points, const double epsilon) {

    ApproxMVBB_STATIC_ASSERT(Derived::RowsAtCompileTime == Dimension);

    Eigen::MatrixBase<Derived> & pp = const_cast< Eigen::MatrixBase<Derived> &>(points);

    // Construct pointer list
    auto size = pp.cols();
    double* * pList = new double*[size];
    for(decltype(size) i=0; i<size; ++i) {
        pList[i] = const_cast<double*>(pp.col(i).data());
    }

    Diameter::typeSegment pairP;
    Diameter::estimateDiameter(&pairP,pList,0,(int)(size-1),Dimension,epsilon);

    Eigen::Map<TVector> p1(pairP.extremity1);
    Eigen::Map<TVector> p2(pairP.extremity2);

    //    std::cout << "p1: " << p1.transpose() << std::endl
    //              << "p2: " << p2.transpose() << std::endl
    //              << " l: " << std::sqrt(pairP.squareDiameter) << std::endl;
    delete[] pList;
    return std::pair<TVector,TVector>(p1,p2);

}

class CompareByAngle {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using PointData = std::pair<unsigned int, bool>;

    /** Cosntructor, points is not a temporary, it accepts all sorts of matrix expressions,
            * however the construction of MatrixRef<> might create a temporary but this is stored in m_p!
            */
    template<typename Derived>
    CompareByAngle(const Eigen::MatrixBase<Derived> & points,
                   const Vector2 &base,
                   unsigned int baseIdx,
                   unsigned int & deletedPoints): m_p(points), m_base(base),m_baseIdx(baseIdx), m_deletedPoints(deletedPoints) {
        EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Derived,2, Eigen::Dynamic)
    }

    /** True if b is positively rotated from a, stricly weak ordering! */
    int operator()( const PointData & point1In,
                    const PointData & point2In )
    {
        using namespace PointFunctions;
        PointData & point1 = const_cast<PointData &>(point1In);
        PointData & point2 = const_cast<PointData &>(point2In);
        unsigned int idx1 = point1.first;
        unsigned int idx2 = point2.first;
        //
        //            if(idx1<idx2){
        //                if(almostEqualUlp(m_p.col(idx1),m_p.col(idx2))){
        //                    if(!point1.second){point1.second = true; ++m_deletedPoints;}
        //                    return false;
        //                }
        //            }else{
        //                if(almostEqualUlp(m_p.col(idx2),m_p.col(idx1))){
        //                    if(!point2.second){point2.second = true; ++m_deletedPoints;}
        //                    return false;
        //                }
        //            }

        // Compare by Area Sign (by ascending positive (z-Axis Rotation) angle in x-y Plane)
        // always  insert the smaller index first , and the larger second (as the function is not completely symmetric!
        if(idx1<idx2) {
            int sgn = orient2d( m_base, m_p.col(idx1), m_p.col(idx2) );
            if(sgn != 0){ return (sgn > 0);}
        } else {
            int sgn = orient2d( m_base, m_p.col(idx2), m_p.col(idx1) );
            if(sgn != 0) {return (sgn < 0);}

        }
        // points are collinear

        if(PointFunctions::equal(m_base, m_p.col(idx1))) return false;
        if(PointFunctions::equal(m_base, m_p.col(idx2))) return true;
        if(PointFunctions::equal(m_p.col(idx1), m_p.col(idx2))) return false;

        // if idx2 lies between mbase and idx1 then it should go after idx1
        return collinearAreOrderedAlongLine(m_base,m_p.col(idx2),m_p.col(idx1));
    }
private:
    unsigned int & m_deletedPoints;
    const Vector2 m_base; const unsigned int m_baseIdx;
    const Eigen::Ref<const Matrix2Dyn> m_p;
};

};


class AABB2d {
public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    AABB2d() {
        reset();
    };

    void reset() {
        // Violating the constraint min<max for making a completey empty box!
        m_minPoint(0) = std::numeric_limits<double>::max();
        m_maxPoint(0) = std::numeric_limits<double>::lowest();
        m_minPoint(1) = std::numeric_limits<double>::max();
        m_maxPoint(1) = std::numeric_limits<double>::lowest();
    }

    template<typename Derived>
    AABB2d& unite(const Eigen::MatrixBase<Derived> &p) {
        EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived,2);
        m_maxPoint(0) = std::max(m_maxPoint(0),p(0));
        m_maxPoint(1) = std::max(m_maxPoint(1),p(1));
        m_minPoint(0) = std::min( m_minPoint(0),p(0));
        m_minPoint(1) = std::min( m_minPoint(1),p(1));
        return *this;
    };

    inline Array2 extent() const{
        return (m_maxPoint - m_minPoint).array();
    };

    inline double maxExtent() const{
        return (m_maxPoint - m_minPoint).maxCoeff();
    };

    using AffineTrafo2d = Eigen::Transform<double,2,Eigen::TransformTraits::Affine>;

    inline bool overlaps(const AABB2d & box) const {
        bool x = (m_maxPoint(0) >= box.m_minPoint(0)) && ( m_minPoint(0) <= box.m_maxPoint(0));
        bool y = (m_maxPoint(1) >= box.m_minPoint(1)) && ( m_minPoint(1) <= box.m_maxPoint(1));
        return (x && y);
    };

    inline bool inside(const Vector2 &p) const {
        return (
                    p(0) >= m_minPoint(0) && p(0) <= m_maxPoint(0) &&
                    p(1) >= m_minPoint(1) && p(1) <= m_maxPoint(1));
    };

    inline bool isEmpty() const {
        return m_maxPoint(0) <= m_minPoint(0) || m_maxPoint(1) <= m_minPoint(1);
    }

    inline void expand(double d) {
        ApproxMVBB_ASSERTMSG(d>=0,"d>=0")
                m_minPoint -= Vector2(d,d);
        m_maxPoint += Vector2(d,d);
    };

    inline void expand(Vector2 d) {
        ApproxMVBB_ASSERTMSG(d(0)>=0 && d(1)>=0,"d>=0")
                m_minPoint -= d;
        m_maxPoint += d;
    };

    inline double area() const {
        Vector2 d = m_maxPoint- m_minPoint;
        return d(0) * d(1);
    };


    AABB2d( const Vector2 &p) {
        m_minPoint = p;
        m_maxPoint = p;
    };
    AABB2d( const Vector2 &l, const Vector2 &u) {
        m_minPoint = Vector2(std::min(l(0),u(0)),std::min(l(1),u(1)));
        m_maxPoint = Vector2(std::max(l(0),u(0)),std::max(l(1),u(1)));
    };

    AABB2d& unite(const AABB2d & box) {
        m_maxPoint(0) = std::max(m_maxPoint(0),box.m_maxPoint(0));
        m_maxPoint(1) = std::max(m_maxPoint(1),box.m_maxPoint(1));
        m_minPoint(0) = std::min( m_minPoint(0),box.m_minPoint(0));
        m_minPoint(1) = std::min( m_minPoint(1),box.m_minPoint(1));
        return *this;
    };

    AABB2d operator+ (const Vector2 &p) {
        AABB2d a;
        a.m_maxPoint(0) = std::max(m_maxPoint(0),p(0));
        a.m_maxPoint(1) = std::max(m_maxPoint(1),p(1));
        a. m_minPoint(0) = std::min( m_minPoint(0),p(0));
        a. m_minPoint(1) = std::min( m_minPoint(1),p(1));
        return a;
    };

    AABB2d operator+ (const AABB2d & box) {
        AABB2d a;
        a.m_maxPoint(0) = std::max(m_maxPoint(0),box.m_maxPoint(0));
        a.m_maxPoint(1) = std::max(m_maxPoint(1),box.m_maxPoint(1));
        a. m_minPoint(0) = std::min( m_minPoint(0),box. m_minPoint(0));
        a. m_minPoint(1) = std::min( m_minPoint(1),box. m_minPoint(1));
        return a;
    };

    AABB2d& operator+= (const AABB2d & box) {
        m_maxPoint(0) = std::max(m_maxPoint(0),box.m_maxPoint(0));
        m_maxPoint(1) = std::max(m_maxPoint(1),box.m_maxPoint(1));
        m_minPoint(0) = std::min( m_minPoint(0),box. m_minPoint(0));
        m_minPoint(1) = std::min( m_minPoint(1),box. m_minPoint(1));
        return *this;
    };


    AABB2d & transform(const AffineTrafo2d & M) {

        AABB2d ret( M*(Vector2( m_minPoint(0), m_minPoint(1))));
        ret.unite(M*(Vector2( m_maxPoint(0), m_minPoint(1))));
        ret.unite(M*(Vector2( m_minPoint(0), m_maxPoint(1))));
        ret.unite(M*(Vector2( m_maxPoint(0), m_maxPoint(1))));
        *this = ret;
        return *this;
    };

    //info about axis aligned bounding box
    Vector2 m_minPoint;
    Vector2 m_maxPoint;
};


class AABB {
public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    AABB() {
        reset();
    };

    template<typename Derived>
    AABB& unite(const Eigen::MatrixBase<Derived> &p) {
        m_maxPoint(0) = std::max(m_maxPoint(0),p(0));
        m_maxPoint(1) = std::max(m_maxPoint(1),p(1));
        m_maxPoint(2) = std::max(m_maxPoint(2),p(2));
        m_minPoint(0) = std::min( m_minPoint(0),p(0));
        m_minPoint(1) = std::min( m_minPoint(1),p(1));
        m_minPoint(2) = std::min( m_minPoint(2),p(2));
        return *this;
    };


    using AffineTrafo = Eigen::Transform<double,3,Eigen::TransformTraits::Affine>;

    inline bool overlaps(const AABB & box) const {
        bool x = (m_maxPoint(0) >= box. m_minPoint(0)) && ( m_minPoint(0) <= box.m_maxPoint(0));
        bool y = (m_maxPoint(1) >= box. m_minPoint(1)) && ( m_minPoint(1) <= box.m_maxPoint(1));
        bool z = (m_maxPoint(2) >= box. m_minPoint(2)) && ( m_minPoint(2) <= box.m_maxPoint(2));
        return (x && y && z);
    };

    inline bool inside(const Vector3 &p) const {
        return (
                    p(0) >= m_minPoint(0) && p(0) <= m_maxPoint(0) &&
                    p(1) >= m_minPoint(1) && p(1) <= m_maxPoint(1) &&
                    p(2) >= m_minPoint(2) && p(2) <= m_maxPoint(2));
    };

    inline Array3 extent() const{
        return (m_maxPoint - m_minPoint).array();
    };

    inline double maxExtent() const{
        return (m_maxPoint - m_minPoint).maxCoeff();
    };

    inline bool isEmpty() const {
        return m_maxPoint(0) <= m_minPoint(0) || m_maxPoint(1) <= m_minPoint(1) || m_maxPoint(2) <= m_minPoint(2);
    }

    inline void expand(double d) {
        ApproxMVBB_ASSERTMSG(d>=0,"d>=0")
                m_minPoint -= Vector3(d,d,d);
        m_maxPoint += Vector3(d,d,d);
    };

    inline void expand(Vector3 d) {
        ApproxMVBB_ASSERTMSG(d(0)>=0 && d(1)>=0 && d(2)>=0,"d>=0")
                m_minPoint -= d;
        m_maxPoint += d;
    };

    inline double volume() const {
        Vector3 d = m_maxPoint- m_minPoint;
        return d(0) * d(1) * d(2);
    };

    //info about axis aligned bounding box
    Vector3 m_minPoint;
    Vector3 m_maxPoint;

    void reset() {
        // Violating the constraint min<max for making a completey empty box!
        m_minPoint(0) = std::numeric_limits<double>::max();
        m_maxPoint(0) = std::numeric_limits<double>::lowest();
        m_minPoint(1) = std::numeric_limits<double>::max();
        m_maxPoint(1) = std::numeric_limits<double>::lowest();
        m_minPoint(2) = std::numeric_limits<double>::max();
        m_maxPoint(2) = std::numeric_limits<double>::lowest();
    }


    AABB( const Vector3 &p) {
        m_minPoint = p;
        m_maxPoint = p;
    };
    AABB( const Vector3 &l, const Vector3 &u) {
        m_minPoint = Vector3(std::min(l(0),u(0)),std::min(l(1),u(1)),std::min(l(2),u(2)));
        m_maxPoint = Vector3(std::max(l(0),u(0)),std::max(l(1),u(1)),std::max(l(2),u(2)));
    };


    AABB& unite(const AABB & box) {
        m_maxPoint(0) = std::max(m_maxPoint(0),box.m_maxPoint(0));
        m_maxPoint(1) = std::max(m_maxPoint(1),box.m_maxPoint(1));
        m_maxPoint(2) = std::max(m_maxPoint(2),box.m_maxPoint(2));
        m_minPoint(0) = std::min( m_minPoint(0),box.m_minPoint(0));
        m_minPoint(1) = std::min( m_minPoint(1),box.m_minPoint(1));
        m_minPoint(2) = std::min( m_minPoint(2),box.m_minPoint(2));
        return *this;
    };

    AABB operator+ (const Vector3 &p) {
        AABB a;
        a.m_maxPoint(0) = std::max(m_maxPoint(0),p(0));
        a.m_maxPoint(1) = std::max(m_maxPoint(1),p(1));
        a.m_maxPoint(2) = std::max(m_maxPoint(2),p(2));
        a. m_minPoint(0) = std::min( m_minPoint(0),p(0));
        a. m_minPoint(1) = std::min( m_minPoint(1),p(1));
        a. m_minPoint(2) = std::min( m_minPoint(2),p(2));
        return a;
    };

    AABB operator+ (const AABB & box) {
        AABB a;
        a.m_maxPoint(0) = std::max(m_maxPoint(0),box.m_maxPoint(0));
        a.m_maxPoint(1) = std::max(m_maxPoint(1),box.m_maxPoint(1));
        a.m_maxPoint(2) = std::max(m_maxPoint(2),box.m_maxPoint(2));
        a. m_minPoint(0) = std::min( m_minPoint(0),box. m_minPoint(0));
        a. m_minPoint(1) = std::min( m_minPoint(1),box. m_minPoint(1));
        a. m_minPoint(2) = std::min( m_minPoint(2),box. m_minPoint(2));
        return a;
    };

    AABB& operator+= (const AABB & box) {
        m_maxPoint(0) = std::max(m_maxPoint(0),box.m_maxPoint(0));
        m_maxPoint(1) = std::max(m_maxPoint(1),box.m_maxPoint(1));
        m_maxPoint(2) = std::max(m_maxPoint(2),box.m_maxPoint(2));
        m_minPoint(0) = std::min( m_minPoint(0),box. m_minPoint(0));
        m_minPoint(1) = std::min( m_minPoint(1),box. m_minPoint(1));
        m_minPoint(2) = std::min( m_minPoint(2),box. m_minPoint(2));
        return *this;
    };

    AABB&  transform(const AffineTrafo & M) {

        AABB ret( M*(Vector3( m_minPoint(0), m_minPoint(1), m_minPoint(2))));
        ret.unite(M*(Vector3( m_maxPoint(0), m_minPoint(1), m_minPoint(2))));
        ret.unite(M*(Vector3( m_minPoint(0), m_maxPoint(1), m_minPoint(2))));
        ret.unite(M*(Vector3( m_minPoint(0), m_minPoint(1), m_maxPoint(2))));
        ret.unite(M*(Vector3( m_minPoint(0), m_maxPoint(1), m_maxPoint(2))));
        ret.unite(M*(Vector3( m_maxPoint(0), m_maxPoint(1), m_minPoint(2))));
        ret.unite(M*(Vector3( m_maxPoint(0), m_minPoint(1), m_maxPoint(2))));
        ret.unite(M*(Vector3( m_maxPoint(0), m_maxPoint(1), m_maxPoint(2))));
        *this = ret;
        return *this;
    };
};


class ConvexHull2D {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /** Cosntructor, points is not a temporary, it accepts all sorts of matrix expressions,
        * however the construction of MatrixRef<> might create a temporary but this is stored in m_p!
        */
    template<typename Derived>
    ConvexHull2D(const Eigen::MatrixBase<Derived> & points) :m_p(points) {
        EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Derived,2, Eigen::Dynamic)
                ApproxMVBB_ASSERTMSG( m_p.data() == points.derived().data() ," You store a temporary in a Ref<> which works here, but do you really want this?")
    }

    void computeMonotonChain();

    //bool verifyHull();

    inline std::vector<unsigned int> & getIndices() {
        return m_indicesCH;
    }

    void chainHull();

    std::vector<unsigned int> m_indicesCH;
    const Eigen::Ref< const Matrix2Dyn > m_p;

    void compute()
    {
        using namespace PointFunctions;

        m_indicesCH.clear();

        // Need at least 1 point! m_p.col(0)
        if(m_p.cols()==0) {
            return;
        }

        // Compute min point p0
        unsigned int position = minPointYX(m_p);
        // std::cout << "min:" << position << std::endl;
        Vector2 base = m_p.col(position);

        // Indices into m_p
        // first = index into m_p, second =  delete flag!
        using PointData = std::pair<unsigned int, bool>;
        std::vector< PointData > indices; indices.reserve(m_p.cols());

        //Save p0 as first point in indices list (temporary)
        indices.emplace_back(position,false);
        //Add all indices of points not almostEqual to position
        for(unsigned int i = 0; i<m_p.cols(); ++i) {
            if( (i != position) && !almostEqualUlp(base,m_p.col(i))) {
                indices.emplace_back(i,false);
            }
        }

        // Convex hull consists of only 1 or 2 points!
        if(indices.size() <= 2  ) {
            for(auto & pa : indices){ m_indicesCH.emplace_back(pa.first);}  return;
        }

        // Sort by angle
        unsigned int deletedPoints = 0;
        CompareByAngle comp(m_p,base,position,deletedPoints);
        std::sort( indices.begin()+1, indices.end(), comp );

        std::vector<unsigned int> indicesT(indices.size()-deletedPoints);
        unsigned int k=0;
        for(auto & p: indices){ if(!p.second){indicesT[k++]=p.first;} }

        // Remove almost equal elements
        // move to front if not almost equal to first point
        // skip all deleted points
        auto ifFunc =  [&](const unsigned int & p1, const unsigned int & p2){
            return !almostEqualUlp(this->m_p.col(p1),this->m_p.col(p2));
        };
        //    auto skipFunc = [](const unsigned int & p1){ return false; };
        auto d2 = ContainerFunctions::moveConsecutiveToFrontIf(indicesT.begin(),indicesT.end(), ifFunc);
        indicesT.resize( std::distance(indicesT.begin(),d2) );

        // Convex hull consists of only 1 or 2 points!
        if(indicesT.size() <= 2  ) {
            for(auto & pa : indicesT){ m_indicesCH.emplace_back(pa);}  return;
        }

        //for(auto  a: indicesT){
        //    std::cout << a.first<<",";
        //}
        //std::cout << "Graham Scan points: " << indicesT.size() << std::endl;
        unsigned int nPoints  = indicesT.size();
        m_indicesCH.reserve(nPoints);
        unsigned int lastIdx  = indicesT[0];
        unsigned int firstIdx = indicesT[1];
        m_indicesCH.push_back( lastIdx  );
        m_indicesCH.push_back( firstIdx );

        unsigned int lPtIdx  = lastIdx;
        unsigned int mPtIdx  = firstIdx;
        unsigned int currIdx ;
        unsigned int i = 2; // current point;

        // skip the first non left turns in the sequence!
        //    std::cout << "lastIdx point: " <<lastIdx << ","<< m_p.col(lastIdx).transpose() << std::endl;
        //    std::cout << "firstIdx point: " <<firstIdx << ","<< m_p.col(firstIdx).transpose() << std::endl;

        while( i<nPoints){
            currIdx = indicesT[i];
            if(leftTurn(m_p.col(lastIdx), m_p.col(firstIdx), m_p.col(currIdx))){
                break;
            }
            ++i;
        };


        //    std::cout << "i:"<< i << std::endl;
        //    std::cout << "currIdx point: " <<currIdx << std::endl;
        //    std::cout << ","<< m_p.col(currIdx).transpose() << std::endl << "===="<<std::endl;
        //    std::cout << "0,5,8: :" << orient2d(m_p.col(0), m_p.col(5), m_p.col(8)) << std::endl;
        //    std::cout << "1,5,8: :" << orient2d(m_p.col(1), m_p.col(5), m_p.col(8)) << std::endl;
        //    std::cout << "5,8,0: :" << orient2d(m_p.col(5), m_p.col(8), m_p.col(0)) << std::endl;

        if ( i < nPoints )
        {
            m_indicesCH.push_back( currIdx );
            decltype(m_indicesCH.rbegin()) revIter;
            lPtIdx  = mPtIdx;
            mPtIdx  = currIdx;

            for (++i ; i < nPoints; ++i )
            {
                currIdx = indicesT[i];

                if ( leftTurn(m_p.col(mPtIdx), m_p.col(currIdx), m_p.col(lastIdx)) )
                {
                    while ( !leftTurn(m_p.col(lPtIdx), m_p.col(mPtIdx), m_p.col(currIdx))   )
                    {
                        //std::cout << "right turn: " <<lPtIdx << ","<< mPtIdx << "," << currIdx << std::endl;
                        ApproxMVBB_ASSERTMSG(m_indicesCH.size()>2,"");
                        m_indicesCH.pop_back();

                        if( m_indicesCH.size() <= 1) { // Degenerate Case if we come back to the beginning
                            mPtIdx  = lPtIdx; // such that lPtIdx stays the same , mPtIdx becomes currIdx and we go to next point!
                            break;
                        }else{
                            mPtIdx = lPtIdx;                      // middle point becomes last
                            revIter = m_indicesCH.rbegin();
                            lPtIdx = *(++revIter);                // previous of .back();
                        }
                    }
                    //std::cout << "left turn: " <<lPtIdx << ","<< mPtIdx << "," << currIdx << std::endl;
                    m_indicesCH.push_back( currIdx );
                    lPtIdx  = mPtIdx;          // last becomes middle
                    mPtIdx  = currIdx;         // middle becomes currIdx
                }/*else{
                  std::cout << "skip point: " << currIdx << std::endl;
              }*/
            }

        }

    }
};

class MinAreaRectangle {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /** Cosntructor, points is not a temporary, it accepts all sorts of matrix expressions,
        * however the construction of MatrixRef<> might create a temporary but this is stored in m_p!
        * MatrixRef<>  m_p is handed further to m_conv
        */
    template<typename Derived>
    MinAreaRectangle(const Eigen::MatrixBase<Derived> & points)
        : m_p(points), m_convh(m_p) {
        EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Derived,2, Eigen::Dynamic)
                ApproxMVBB_ASSERTMSG( m_p.data() == points.derived().data() ," You store a temporary in a Ref<> which works here, but do you really want this?")
    }

    struct Box2d {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        inline void reset() {
            m_p.setZero();
            m_u.setZero();
            m_v.setZero();
            m_area = 0.0;
            m_uL = 1.0;
            m_vL = 1.0;
        }
        Vector2 m_p;   ///< first corner x = m_p0 + m_u*m_uL * u + m_v*m_vL * v  , u,v in [0,1]
        Vector2 m_u;   ///< vector of first side (x-Axis)  (normalized)
        Vector2 m_v;   ///< vector of second side (y-Axis) (normalized)

        double m_uL = 1.0;     ///< Scalar factor direction u
        double m_vL = 1.0;     ///< Scalar factor direction v

        double m_area = 0.0;
    };

    inline const Box2d &
    getMinRectangle() {
        return m_minBox;
    }


    void compute() {
        // Compute minimum area rectangle

        // Clear minBox;
        m_minBox.reset();

        // If no points return early!
        if(m_p.cols() == 0) {
            return;
        }

        // Generate Convex Hull
        m_convh.compute();
        ApproxMVBB_ASSERTMSG(m_convh.verifyHull(), "Convex hull not ok!")

        // Compute Rectangle
        computeRectangle();
    }

    void computeRectangle() {

        // Copy Points;
        m_hullIdx = m_convh.getIndices();
        unsigned int nPoints = m_hullIdx.size();

        //std::cout << "Convex hull points:" << nPoints << std::endl;

        // The below code works for points >= 2
        // Anyway catch the cases n=1 and 2 and return early!
        if(nPoints == 0) {
            return;
        } else if( nPoints == 1) {
            m_minBox.m_p = m_p.col(m_hullIdx[0]);
            m_minBox.m_u.setZero();
            m_minBox.m_v.setZero();
            adjustRectangle();
            return;
        } else if(nPoints == 2) {
            m_minBox.m_p = m_p.col(m_hullIdx[0]);
            m_minBox.m_u = m_p.col(m_hullIdx[1]) - m_p.col(m_hullIdx[0]);
            m_minBox.m_v.setZero();
            adjustRectangle();
            return;
        }


        // Performing Rotating Calipers method
        // 1. find the vertices with the minimum and maximum x an y coordinates.
        //    These vertices (points) will be denoted by {PI, PJ, PK, PL}
        // 2. construct CALL and CALJ as the first set of calipers parallel to the x-axis,
        //    and CALI and CALK as the second set of calipers parallel to the y-axis
        // 3. create for every CALX a next caliper,
        //    NEXT_CALX, based on the next point in the convex hull and calculate the tangent of CALX and it's NEXT_CALX
        // 4. get the smallest positive tangent between all CALX and NEXT_CALX and rotate every caliper by that smallest gradient
        // 5. repeat step 3 and 4 untill the calipers have turned more then 90 degrees

        // We implement this algorithm a bit differently:
        // 1. Compute all angles of all edges angle[i] is angle of edge p[i],p[i+1].
        // 2. Rotatate the first caliper A from edge to edge.
        // 3. Determine all new angles for B,C,D such that it forms a box (rotation of pi from angle of caliper A)
        // 4. Determine for each caliper I  the location (the corresponding vertex index) where the
        //    caliper I forms a tangent (starting from the currend index).
        // 5. Determine the box by intersecting the calipers.
        // 6. Compare to area to minimum and return to 2.

        // Calcualte all edge angles

        m_angles.resize(nPoints);
        for(unsigned int i=0; i<nPoints; ++i) {
            m_angles[i] = PointFunctions::getAngle(m_p.col(m_hullIdx[i]),
                                                   m_p.col(m_hullIdx[(i+1) % nPoints])
                                                  );
            //std::cout << "angle: " << m_angles[i] << std::endl;
        }

        // 2. Construct 4 Calipers
        Caliper calipers[4];
        // Set A
        calipers[0].m_idx = 0;
        calipers[0].m_ptIdx = m_hullIdx[0];
        // Init A,B,C,D
        updateCalipers(m_angles[0],calipers);

        Box2d box;
        // Init minBox
        getBox(calipers,m_minBox);

        // Rotate calipers over all edges and take minimum box
        for(unsigned int i=1; i<nPoints; ++i) {

            updateCalipers(m_angles[i],calipers);
            // Get Box from Calipers
            getBox(calipers,box);

            if(box.m_area < m_minBox.m_area) {

                m_minBox = box;
            }
        }


        adjustRectangle();

    }

    struct Caliper {
        unsigned int m_idx = 0;   // index in m_hullIdx
        unsigned int m_ptIdx = 0; // index in m_p
        double m_currAngle = 0.0;
    };

    // determine the vertex v for which the edge angle is greater than c.m_currAngle
    // the find algroithm starts at c.m_idx
    void findVertex(Caliper & c) {

        double matchAngle = c.m_currAngle;    // between 0-2pi
        //std::cout << "Match Angle: " << matchAngle << std::endl;

        unsigned int currIdx = c.m_idx;

        double currAngle = m_angles[currIdx]; // all between 0-2pi
        double nextAngle;

        unsigned int nPoints = m_hullIdx.size();
        bool found = false;
        unsigned int i = 0;
        while(!found && i < nPoints) {

            currIdx = (currIdx + 1) % nPoints;
            // std::cout << "c: " << currIdx << std::endl;
            nextAngle = m_angles[currIdx];
            //std::cout << "diff: " << std::abs(currAngle-nextAngle) << std::endl;
            //                                         curr       next
            // if we are not at the boundary [ ..., 45 degree, 50 degree , ....]
            // and in
            if( nextAngle > currAngle ) {
                //std::cout << " greater" <<std::endl;
                found = (currAngle < matchAngle && matchAngle <= nextAngle);
            } else {
                //std::cout << " boundary" <<std::endl;
                //                                     curr       next
                // if we are at the boundary [ ..., 359 degree] [ 5 degree , .... ]
                // skip this if we have an angle difference |curr -next| < eps
                // this might happen due to floating point nonsense
                found = (std::abs(currAngle-nextAngle)>1e-10) && (currAngle < matchAngle || matchAngle <= nextAngle);
            }

            currAngle = nextAngle;
            ++i;
        }

        if(found) {
            c.m_idx = currIdx;
            c.m_ptIdx = m_hullIdx[currIdx];
            //std::cout << "caliper idx:" << currIdx << std::endl;
        } else {
            std::stringstream ss;
            for(auto & a : m_angles){ ss << a <<",";}
            ApproxMVBB_ERRORMSG("Could not find vertex with angle greater than: "
                     << matchAngle << "in angles: " << ss.str());
        }
    }

    void getBox(Caliper (&c)[4], Box2d & box) {
            using namespace PointFunctions;
            // Intersect Caliper A with B,
            box.m_u = intersectLines(   m_p.col(c[0].m_ptIdx), c[0].m_currAngle,
                                        m_p.col(c[1].m_ptIdx), c[1].m_currAngle);
            // Intersect Caliper D with C,
            box.m_v = intersectLines(   m_p.col(c[3].m_ptIdx), c[3].m_currAngle,
                                        m_p.col(c[2].m_ptIdx), c[2].m_currAngle);

            // Intersect Caliper A with D,
            box.m_p = intersectLines(  m_p.col(c[0].m_ptIdx), c[0].m_currAngle,
                                       m_p.col(c[3].m_ptIdx), c[3].m_currAngle);

            box.m_u -= box.m_p; //(p1-p0)
            box.m_v -= box.m_p; //(p2-p0)
            box.m_area = box.m_u.norm()*box.m_v.norm();

        }

    using Vector2U = Eigen::Matrix<unsigned int,2,1>;

    inline void adjustRectangle(){

        // The rectangle might be a slight parallelogram due to numerics
        //     |-l---*
        //  *  +++++++=================+
        //  |  +     |               | +
        //  h  +   | <-v  Rect     |   +
        //  |  + |               |     +
        //  -  +=================+++++++
        //     p        u

        double uNorm = m_minBox.m_u.norm();
        double vNorm = m_minBox.m_v.norm();

        // Normalize u direction (if u close to zero -> this becomes not finite)
        Vector2 uN = m_minBox.m_u.array() / uNorm;
        Vector2 vN = m_minBox.m_v.array() / vNorm;

        bool uF = uN.allFinite();
        bool vF = vN.allFinite();

        if(uF && vF){

            // make orthogonal (x axis is u)
            Vector2 uNT;
            uNT(0) = -uN(1);
            uNT(1) =  uN(0);

            double h =  uNT.dot(m_minBox.m_v);
            double l =  uN.dot (m_minBox.m_v);

            if(l>=0.0){
                m_minBox.m_uL = uNorm + l;
            }else{
                m_minBox.m_uL = uNorm - l;
                m_minBox.m_p += uN*l; // move start point back
            }

            // if v vector pointed downwards (negative h)
            if( h<0 ){
                // invert uNT (and switch u,v)
                uNT *= -1.0;
                m_minBox.m_u = uNT;
                m_minBox.m_v = uN;
                m_minBox.m_vL = m_minBox.m_uL;
                m_minBox.m_uL = -h;
            }else{
                m_minBox.m_vL = h;
                m_minBox.m_u = uN;
                m_minBox.m_v = uNT;
            }

        }else if(uF && !vF){
            // adjust v direction
            m_minBox.m_v(0) = -m_minBox.m_u(1);
            m_minBox.m_v(1) = m_minBox.m_u(0);

            m_minBox.m_uL = uNorm;
            m_minBox.m_vL = 0.0;

        }else if(!uF && vF){
            // adjust u direction
            m_minBox.m_u(0) = -m_minBox.m_v(1);
            m_minBox.m_u(1) = m_minBox.m_v(0);

            m_minBox.m_uL = 0.0;
            m_minBox.m_vL = vNorm;
        }else{
            // adjust both directions
            m_minBox.m_u(0) = 1.0; m_minBox.m_u(1) = 0.0;
            m_minBox.m_v(0) = 0.0;  m_minBox.m_v(1) = 1.0;
            m_minBox.m_uL = 0.0; m_minBox.m_vL = 0.0;
        }

    }

    inline void updateCalipers(double edgeAngle, Caliper (&c)[4]) {

        updateAngles(edgeAngle, c);
        for(unsigned char i=0; i<4; i++) {
            findVertex(c[i]);
        }

    }


    //determine caliper angles according to the box given with c[0].m_currAngle = edgeAngle;
    inline void updateAngles(double edgeAngle, Caliper (&c)[4]) {
        c[0].m_currAngle = AngleFunctions::mapTo2Pi(edgeAngle);
        c[1].m_currAngle = AngleFunctions::mapTo2Pi(c[0].m_currAngle + 0.5*M_PI);
        c[2].m_currAngle = AngleFunctions::mapTo2Pi(c[1].m_currAngle + 0.5*M_PI);
        c[3].m_currAngle = AngleFunctions::mapTo2Pi(c[2].m_currAngle + 0.5*M_PI);

        //        std::cout << "caliper 1 angle:" <<  c[0].m_currAngle << std::endl;
        //        std::cout << "caliper 2 angle:" <<  c[1].m_currAngle << std::endl;
        //        std::cout << "caliper 3 angle:" <<  c[2].m_currAngle << std::endl;
        //        std::cout << "caliper 4 angle:" <<  c[3].m_currAngle << std::endl;
    }

    std::vector<unsigned int> m_hullIdx;

    std::vector<double> m_angles;
    Box2d m_minBox;
    const Eigen::Ref<const Matrix2Dyn> m_p;

    ConvexHull2D m_convh;

};

class ProjectedPointSet {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW


    /** Computes the true MVBB by constructing the convex hull in 2d
        * and then building the minimum area rectangle around it and afterwards
        * pulling the rectangle out in direction m_dirZ again which build the MVBB in direction m_dirZ
        */
    template<typename Derived>
    OOBB computeMVBB(const Vector3 & zDir,
                     const Eigen::MatrixBase<Derived> & points ) {

        EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Derived,3,Eigen::Dynamic)

                //using namespace CoordinateSystem;
                m_zDir = zDir;
        project(points);

        // compute minimum area rectangle first
        // std::cout << "Dump points DEBUG:" << std::endl;
        //TestFunctions::dumpPointsMatrixBinary("DumpedPoints.bin",m_p);


        MinAreaRectangle mar(m_p);
        mar.compute();
        auto rect = mar.getMinRectangle();

        //std::cout << "Dump RECT DEBUG:" << std::endl;
        // Vector2List p;
        //p.push_back( rect.m_p);
        //p.push_back( rect.m_p + rect.m_u );
        //p.push_back( rect.m_p + rect.m_u + rect.m_v );
        //p.push_back( rect.m_p + rect.m_v );
        //p.push_back( rect.m_p);

        //TestFunctions::dumpPoints("./MinAreaRectangleTest13" "Out.txt",p);


        // Box coordinates are in K Frame

        Matrix22 A2_KM;

        //        std::cout << "u:" << rect.m_u.norm() << std::endl;
        //        std::cout << "v:" << rect.m_v.norm() << std::endl;

        A2_KM.col(0) = rect.m_u;
        A2_KM.col(1) = rect.m_v;

        Vector2 M_p = A2_KM.transpose()*rect.m_p;

        Vector3 M_min;
        M_min.head<2>() = M_p;
        M_min(2) = m_minZValue;

        Vector3 M_max(rect.m_uL, rect.m_vL, 0.0);
        M_max.head<2>() += M_p;
        M_max(2) = m_maxZValue;

        // Make coordinate transformation from M frame (Minimum Rectancle)
        // to K frame (Projection Plane);
        Matrix33 A_IM;
        // Make A_KM
        A_IM.setIdentity();
        A_IM.block<2,2>(0,0) =  A2_KM;
        // Make A_IM;
        A_IM = m_A_KI.transpose()*A_IM; // A_IM = A_IK * A_KM

        return OOBB(M_min,M_max, A_IM );
    }

    inline void makeCoordinateSystem(      Vector3 &v1,
                                           Vector3 &v2,
                                           Vector3 &v3)
    {

        using std::abs;
        using std::sqrt;

        v1.normalize();

        if(abs(v1(0)) > abs(v1(2))){
            double invLen = 1.0 / sqrt(v1(0)*v1(0) + v1(2) * v1(2));
            v2 = Vector3(-v1(2) * invLen, 0, v1(0) *invLen);
        }
        else{
            double invLen = 1.0 / sqrt(v1(1)*v1(1) + v1(2) * v1(2));
            v2 = Vector3(0, v1(2) * invLen, -v1(1) *invLen);
        }
        v3 = v1.cross(v2);

        v2.normalize();
        v3.normalize();

    };

    template<typename Derived>
    void project(const Eigen::MatrixBase<Derived> & points) {

        if(points.cols()==0) {
            ApproxMVBB_ERRORMSG("Point set empty!");
        }

        // Generate Orthonormal Bases
        Vector3 xDir,yDir;
        //std::cout << "dir: " <<  m_zDir << std::endl;
        makeCoordinateSystem(m_zDir,xDir,yDir);

        //Make coodinate transform from frame I to K!
        m_A_KI.col(0) = xDir;
        m_A_KI.col(1) = yDir;
        m_A_KI.col(2) = m_zDir;
        m_A_KI.transposeInPlace();
        //ApproxMVBB_ASSERTMSG(checkOrthogonality(xDir,yDir,m_zDir,1e-6),
        //          "Not orthogonal: x:"<<xDir.transpose()<< " y: "<< yDir.transpose() << " z: "  << m_zDir.transpose());

        // Project Points onto xDir,yDir Halfspace
        auto size = points.cols();
        m_p.resize(2,size);
        //m_p = m_A_KI * points;  // Project points! (below is faster)
        Vector3 p;
        m_maxZValue = std::numeric_limits<double>::lowest();
        m_minZValue = std::numeric_limits<double>::max();
        for(decltype(size) i=0; i<size; ++i) {
            p = m_A_KI*points.col(i);
            m_p.col(i) = p.head<2>();
            m_maxZValue = std::max(m_maxZValue,p(2));
            m_minZValue = std::min(m_minZValue,p(2));
        }


    }

    Matrix2Dyn  m_p; ///< Projected points in frame K

    Vector3 m_zDir;
    Matrix33 m_A_KI; ///< Transformation from I frame into the projection frame K
    double m_minZValue, m_maxZValue;


    template<typename Derived>
    OOBB computeMVBBApprox(const Vector3 & zDir,
                           const Eigen::MatrixBase<Derived> & points,
                           const double epsilon) {

        EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Derived,3,Eigen::Dynamic)

                using namespace PointFunctions;
        m_zDir = zDir;
        project(points);
        //std::cout <<"projected points" <<std::endl;

        // Estimate diameter in 2d projective plane
        std::pair<Vector2,Vector2> pp = estimateDiameter<2,Vector2>(m_p,epsilon);


        Vector2 dirX = pp.first - pp.second;
        if( ( pp.second.array() >=  pp.first.array()).all() ) {
            dirX *= -1;
        }
        //std::cout <<"estimated 2d diameter: " << dirX.transpose() << " eps: " << epsilon << std::endl;
        // Built Coordinate Trafo from frame K to frame M
        using Matrix22 = Eigen::Matrix<double, 2, 2>;
        Matrix22 A2_MK;
        dirX.normalize();

        // If normalized direction INf/NaN, use (1,0)
        if( !dirX.allFinite() ) {
            dirX.setZero();
            dirX(0)= 1;
        }

        Vector2 dirY(-dirX(1),dirX(0));  // Positive rotation of 90 degrees
        A2_MK.col(0) = dirX;
        A2_MK.col(1) = dirY;
        A2_MK.transposeInPlace();

        AABB2d aabb;
        Vector2 p;
        auto size = m_p.cols();
        for(unsigned int i=0; i < size; ++i) {
            p.noalias() = A2_MK * m_p.col(i); // Transform all points
            aabb.unite(p);
        }

        //         std::cout << "aabb_min: "<<  aabb.m_minPoint.transpose() << std::endl;
        //         std::cout << "aabb_max: "<< aabb.m_maxPoint.transpose() << std::endl;

        Vector3 M_min;
        M_min.head<2>() = aabb.m_minPoint;
        M_min(2) = m_minZValue;

        Vector3 M_max;
        M_max.head<2>() = aabb.m_maxPoint;
        M_max(2) = m_maxZValue;

        // Make coordinate transformation from M frame (Minimum Rectancle)
        // to K frame (Projection Plane);
        Matrix33 A_KM;
        A_KM.setIdentity();
        A_KM.block<2,2>(0,0) =  A2_MK.transpose();

        //         std::cout << "M_min: "<< M_min.transpose() << std::endl;
        //         std::cout << "M_max: "<< M_max.transpose() << std::endl;

        return OOBB(M_min,M_max, m_A_KI.transpose()*A_KM );
    }
};


/**
* Function to optimize oriented bounding box volume.
* This performs an exhaustive grid search over a given tighly fitted bounding box (use approximateMVBBDiam)
* to find a tighter volume.
* @param gridSize of the 3d Grid
* @param optLoops how many optimization loops are preformed
*        for the oobb computed in the given discrete sampled direction in the grid  (see optimizeMVBB)
*/
template<typename Derived>
OOBB approximateMVBBGridSearch(const Eigen::MatrixBase<Derived> & points,
                               OOBB oobb,
                               double epsilon,
                               const unsigned int gridSize = 5,
                               const unsigned int optLoops = 6
        ) {
    EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Derived,3,Eigen::Dynamic)

            // Define the volume lower bound above we accept a new volume as
            double volumeAcceptTol = oobb.volume() * 1e-6;

    //Get the direction of the input OOBB in I frame:
    Vector3 dir1 = oobb.getDirection(0);
    Vector3 dir2 = oobb.getDirection(1);
    Vector3 dir3 = oobb.getDirection(2);

    Vector3 dir;

    ProjectedPointSet proj;

    for(int x = -int(gridSize); x <= (int)gridSize; ++x ) {
        for(int  y = -int(gridSize); y <= (int)gridSize; ++y ) {
            for(int z = 0; z <= (int)gridSize; ++z ) {


                if( MathFunctions::gcd3(x,y,z)> 1 ||  ((x==0) && (y==0) &&  (z==0))  ) {
                    continue;
                }

                // Make direction
                dir = x*dir1 + y*dir2 + z*dir3;
                //std::cout << "dir: " << dir.transpose() << std::endl;

                // Compute MVBB in dirZ
                auto res = proj.computeMVBB(dir,points);

                if(optLoops){
                    res = optimizeMVBB(points,res,optLoops);
                }

                if(res.volume() < oobb.volume() && res.volume()>volumeAcceptTol ) {
                    oobb = res;
                }

            }
        }
    }
    return oobb;
}



/** We are given a point set, and (hopefully) a tight fitting
*   bounding box. We compute a sample of the given size nPoints that represents
*   the point-set. The only guarenteed is that if we use sample of size m,
*   we get an approximation of quality about 1/\sqrt{m}. Note that we pad
*   the sample if necessary to get the desired size.
*   This function changes the oobb and sets the z Axis to the greates extent!
*   @param nPoints needs to be greater or equal than 2
*/
template<typename Derived>
void samplePointsGrid(Matrix3Dyn & newPoints,
                      const Eigen::MatrixBase<Derived> & points,
                      const unsigned int nPoints,
                      OOBB & oobb) {


    if(nPoints >= points.cols() || nPoints < 2) {
        ApproxMVBB_ERRORMSG("Wrong arguements!")
    }

    newPoints.resize(3,nPoints);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<unsigned int> dis(0, points.cols()-1);

    //total points = bottomPoints=gridSize^2  + topPoints=gridSize^2
    unsigned int gridSize = std::max( static_cast<unsigned int>( std::sqrt( static_cast<double>(nPoints) / 2.0 )) , 1U );

    // Set z-Axis to longest dimension
    //std::cout << oobb.m_minPoint.transpose() << std::endl;
    oobb.setZAxisLongest();

    unsigned int halfSampleSize = gridSize*gridSize;
    std::vector< std::pair<unsigned int , double > > topPoints(halfSampleSize,    std::pair<unsigned int,double>{} );    // grid of indices of the top points (indexed from 1 )
    std::vector< std::pair<unsigned int , double > > bottomPoints(halfSampleSize, std::pair<unsigned int,double>{} ); // grid of indices of the bottom points (indexed from 1 )

    using LongInt = long long int;
    Eigen::Array<LongInt, 2, 1> idx; // Normalized P
    //std::cout << oobb.extent() << std::endl;
    //std::cout << oobb.m_minPoint.transpose() << std::endl;
    Array2 dxdyInv =  Array2(gridSize,gridSize) / oobb.extent().head<2>(); // in K Frame;
    Vector3 K_p;

    Matrix33 A_KI(oobb.m_q_KI.matrix().transpose());

    // Register points in grid
    auto size = points.cols();
    for(unsigned int i=0; i<size; ++i) {

        K_p = A_KI * points.col(i);
        // get x index in grid
        idx = (  (K_p - oobb.m_minPoint).head<2>().array() * dxdyInv ).template cast<LongInt>();
        // map to grid
        idx(0) = std::max(   std::min( LongInt(gridSize-1), idx(0)),  0LL   );
        idx(1) = std::max(   std::min( LongInt(gridSize-1), idx(1)),  0LL   );
        //std::cout << idx.transpose() << std::endl;
        unsigned int pos = idx(0) + idx(1)*gridSize;

        // Register points in grid
        // if z axis of p is > topPoints[pos]  -> set new top point at pos
        // if z axis of p is < bottom[pos]     -> set new bottom point at pos

        if( topPoints[pos].first == 0) {
            topPoints[pos].first  = bottomPoints[pos].first  = i+1;
            topPoints[pos].second = bottomPoints[pos].second = K_p(2);
        } else {
            if( topPoints[pos].second < K_p(2) ) {
                topPoints[pos].first = i+1;
                topPoints[pos].second = K_p(2);
            }
            if( bottomPoints[pos].second > K_p(2) ) {
                bottomPoints[pos].first = i+1;
                bottomPoints[pos].second = K_p(2);
            }
        }
    }

    // Copy top and bottom points
    unsigned int k=0;

    // k does not overflow -> 2* halfSampleSize = 2*gridSize*gridSize <= nPoints;
    for(unsigned int i=0; i<halfSampleSize; ++i) {
        if( topPoints[i].first != 0 ) {
            // comment in if you want the points top points of the grid
            //            Array3 a(i % gridSize,i/gridSize,oobb.m_maxPoint(2)-oobb.m_minPoint(2));
            //            a.head<2>()*=dxdyInv.inverse();
            newPoints.col(k++) =  points.col(topPoints[i].first-1);  //  A_KI.transpose()*(oobb.m_minPoint + a.matrix()).eval() ;
            if(topPoints[i].first != bottomPoints[i].first) {
                // comment in if you want the bottom points of the grid
                //                Array3 a(i % gridSize,i/gridSize,0);
                //                a.head<2>()*=dxdyInv.inverse();
                newPoints.col(k++) = points.col(bottomPoints[i].first-1); //  A_KI.transpose()*(oobb.m_minPoint + a.matrix()).eval() ;
            }
        }
    }
    // Add random points!
    while( k < nPoints) {
        newPoints.col(k++) = points.col( dis(gen) ); //= Vector3(0,0,0);//
    }
}

/**
* Function to optimize oriented bounding box volume.
* Projecting several times into the direction of the axis of the current oobb,
* constructing the mvbb and overwriting the current oobb if volume is smaller
*/
template<typename Derived>
OOBB optimizeMVBB( const Eigen::MatrixBase<Derived> & points,
                   OOBB oobb, unsigned int times = 10) {

    EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Derived,3,Eigen::Dynamic)

            if( oobb.volume() == 0.0 || times == 0) {
        return oobb;
    }

    // Define the volume lower bound above we accept a new volume as
    double volumeAcceptTol = oobb.volume() * 1e-6;

    bool sameAsCache = true;
    unsigned int cacheIdx = 0; // current write Idx into the cache
    Vector3 dirCache[3]; // save the last three directions (avoiding cycling in choosen axis)

    Vector3 dir;
    ProjectedPointSet proj;
    for(unsigned int loop = 0; loop < times; ++loop ) {

        // Determine Direction (choose x or y axis)
        //std::cout << oobb.m_q_KI.matrix() << std::endl;
        dir = oobb.getDirection(0);

        // check against all chache values
        for(unsigned char i=0; i<3 && i<loop; ++i) {
            double dotp = std::abs(dir.dot(dirCache[i])); //
            if( std::abs(dotp - 1.0) <= 1e-3 ) {
                //std::cout << "Change Dir" << std::endl;
                // direction are almost the same as in the cache, choose another one
                dir = oobb.getDirection(1);
                break;
            }
        }
        // Write to cache and shift write idx
        dirCache[cacheIdx] = dir;
        cacheIdx = (cacheIdx + 1) % 3;

        //std::cout << "Optimizing dir: " << dir << std::endl;
        OOBB o = proj.computeMVBB( dir, points);

        if( o.volume() < oobb.volume() && o.volume()>volumeAcceptTol) {
            oobb = o;
        }
    }

    return  oobb;
}


template<typename Derived>
OOBB approximateMVBB(const Eigen::MatrixBase<Derived> & points,
                     const double epsilon,
                     const unsigned int pointSamples,
                     const unsigned int gridSize,
                     const unsigned int mvbbDiamOptLoops,
                     const unsigned int gridSearchOptLoops
        ) {
    EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Derived,3,Eigen::Dynamic)

            // Approx MVBB with Diameter
            auto oobb = approximateMVBBDiam(points,epsilon,mvbbDiamOptLoops);

    if(pointSamples<points.cols()) {

        // sample points
        Matrix3Dyn sampled;
        samplePointsGrid(sampled,points,pointSamples,oobb);

        // Exhaustive grid search with sampled points
        oobb = approximateMVBBGridSearch(sampled,oobb,epsilon,gridSize,gridSearchOptLoops);

    } else {
        oobb = approximateMVBBGridSearch(points,oobb,epsilon,gridSize,gridSearchOptLoops);
    }
    return oobb;
}


/**
* Function to optimize oriented bounding box volume.
* This constructs an approximation of a tightly fitted bounding box by computing
* the diameter d in 3d and afterwards the projection of the points in the plane perpendicular to direction d
* and then the diameter f in 2d and extruding the OOBB in 2d to the final OOBB approximation in 3d.
*/
template<typename Derived>
OOBB approximateMVBBDiam(const Eigen::MatrixBase<Derived> & points,
                         const double epsilon,
                         const unsigned int optLoops = 10
        ) {
    EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Derived,3,Eigen::Dynamic)

            using namespace PointFunctions;

    auto pp = estimateDiameter<3,Vector3>(points,epsilon);

    Vector3 dirZ = pp.first - pp.second;
    if( ( pp.second.array() >=  pp.first.array()).all() ) {
        dirZ *= -1;
    }
    // If direction zero, use (1,0)
    if( (dirZ.array() == 0.0).all() ) {
        dirZ.setZero();
        dirZ(0)= 1;
    }
    //std::cout <<"estimated 3d diameter: " << dirZ.transpose() << " eps: " << epsilon << std::endl;


    // Compute MVBB in dirZ
    ProjectedPointSet proj;
    //OOBB oobb = proj.computeMVBB();
    // or faster estimate diameter in projected plane and build coordinate system
    OOBB oobb = proj.computeMVBBApprox(dirZ,points,epsilon);

    if(optLoops) {
        oobb = optimizeMVBB(points,oobb,optLoops);
    }
    return oobb;
}

}

template ApproxMVBB::OOBB ApproxMVBB::approximateMVBBDiam(const Eigen::MatrixBase< Eigen::Matrix<double,3,Eigen::Dynamic> > & points,
        const double epsilon,
        const unsigned int pointSamples,
        const unsigned int gridSize,
        const unsigned int mvbbDiamOptLoops,
        const unsigned int gridSearchOptLoops);

