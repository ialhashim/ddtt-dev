// From: https://github.com/pekrau/MolScript
#pragma once
#undef ico1
#undef ico2
#undef ico3

namespace icosahedron{
	struct vector3{
		double x, y, z;
	};

	std::vector<vector3> sample(int level)
	{
		auto v3_initialize = [&](vector3 *dest, const double x, const double y, const double z){
			dest->x = x;
			dest->y = y;
			dest->z = z;
		};

		auto v3_between = [&](vector3 *dest, const vector3 *p1, const vector3 *p2, const double fraction){
			dest->x = fraction * (p2->x - p1->x) + p1->x;
			dest->y = fraction * (p2->y - p1->y) + p1->y;
			dest->z = fraction * (p2->z - p1->z) + p1->z;
		};

		auto v3_normalize = [&](vector3 *v){
			register double ir;

			ir = 1.0 / sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
			v->x *= ir;
			v->y *= ir;
			v->z *= ir;
		};

		vector3 *icosahedron = NULL;
		vector3 *current;
		int count;
		int parts;
		double fraction;

		/*------------------------------------------------------------*/
		auto sphere_triangles = [&](vector3 *ico1, vector3 *ico2, vector3 *ico3)
		{
			int slot;
			vector3 p1, p2;

			v3_between(&p1, ico1, ico3, fraction);
			v3_normalize(&p1);
			v3_between(&p2, ico2, ico3, fraction);
			v3_normalize(&p2);
			for (slot = 1; slot <= parts; slot++) {
				v3_between(current + count, &p1, &p2, (double)slot / (double)parts);
				v3_normalize(current + count);
				count++;
			}
		};

		/*------------------------------------------------------------*/
		auto sphere_ico_point_count = [&](int level)
			/*
			Return the number of points generated for the sphere given
			the level, using the icosahedron algorithm.
			*/
		{
			int i, count;

			count = 12;			/* vertices */
			count += 30 * (level - 1);	/* edges */
			for (i = 1; i < level - 1; i++) count += 20 * i; /* faces */

			return count;
		};

		auto icosahedron_vertices = [&]()
			/*
			Return a new array containing the 12 vertices of an icosahedron.
			The vertices are on the unit sphere. The order of the
			vertices are essential for other dependent routines.
			*/
		{
			double a = 0.850650808352039932;
			double b = 0.525731112119133606;
			vector3 *v, *p;
			int i, j;

			p = v = (vector3 *)malloc(12 * sizeof(vector3));

			for (i = 1; i <= 2; i++) {
				a = -a;
				for (j = 1; j <= 2; j++) {
					b = -b;
					v3_initialize(p++, 0.0, a, b);
					v3_initialize(p++, b, 0.0, a);
					v3_initialize(p++, a, b, 0.0);
				}
			}

			return v;
		};

		/*------------------------------------------------------------*/
		auto sphere_ico_points = [&](int level)
			/*
			  Return an array of points on the unit sphere, using the
			  icosahedron algorithm. The level determines the number of points.
			  */
		{
			int vertex;

			if (icosahedron == NULL) icosahedron = icosahedron_vertices();
			
			//vector3 *new_ico;
			//new_ico = (vector3 *)malloc(sphere_ico_point_count(level) * sizeof(vector3));
			std::vector<vector3> samples(sphere_ico_point_count(level));
			vector3 *new_ico = &samples[0];

			current = new_ico;
			count = 0;

			*(current + count++) = *icosahedron;
			*(current + count++) = *(icosahedron + 9);

			for (vertex = 0; vertex < level; vertex++) {
				fraction = (double)vertex / (double)level;
				parts = level - vertex;
				sphere_triangles(icosahedron + 1, icosahedron + 4, icosahedron + 0);
				sphere_triangles(icosahedron + 4, icosahedron + 8, icosahedron + 0);
				sphere_triangles(icosahedron + 8, icosahedron + 3, icosahedron + 0);
				sphere_triangles(icosahedron + 3, icosahedron + 2, icosahedron + 0);
				sphere_triangles(icosahedron + 2, icosahedron + 1, icosahedron + 0);
				sphere_triangles(icosahedron + 5, icosahedron + 7, icosahedron + 9);
				sphere_triangles(icosahedron + 7, icosahedron + 10, icosahedron + 9);
				sphere_triangles(icosahedron + 10, icosahedron + 11, icosahedron + 9);
				sphere_triangles(icosahedron + 11, icosahedron + 6, icosahedron + 9);
				sphere_triangles(icosahedron + 6, icosahedron + 5, icosahedron + 9);
			}

			for (vertex = 1; vertex < level; vertex++) {
				fraction = (double)vertex / (double)level;
				parts = level - vertex;
				sphere_triangles(icosahedron + 2, icosahedron + 1, icosahedron + 5);
				sphere_triangles(icosahedron + 5, icosahedron + 6, icosahedron + 1);
				sphere_triangles(icosahedron + 1, icosahedron + 4, icosahedron + 6);
				sphere_triangles(icosahedron + 6, icosahedron + 11, icosahedron + 4);
				sphere_triangles(icosahedron + 4, icosahedron + 8, icosahedron + 11);
				sphere_triangles(icosahedron + 11, icosahedron + 10, icosahedron + 8);
				sphere_triangles(icosahedron + 8, icosahedron + 3, icosahedron + 10);
				sphere_triangles(icosahedron + 10, icosahedron + 7, icosahedron + 3);
				sphere_triangles(icosahedron + 3, icosahedron + 2, icosahedron + 7);
				sphere_triangles(icosahedron + 7, icosahedron + 5, icosahedron + 2);
			}

			assert(count == sphere_ico_point_count(level));

			// Clean up
			delete icosahedron;

			return samples;
		};

		return sphere_ico_points(level);
	}
}
