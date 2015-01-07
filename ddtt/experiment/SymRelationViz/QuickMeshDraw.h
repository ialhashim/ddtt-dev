#pragma once

#include <qgl.h>
#include "SurfaceMeshModel.h"

using namespace SurfaceMesh;

#define glColorQt(c) glColor4d(c.redF(), c.greenF(), c.blueF(), c.alphaF())

struct QuickMeshDraw{

	static void drawMeshSolid(SurfaceMeshModel * mesh, QColor c = QColor(255, 255, 255, 255), Surface_mesh::Vector3 translation = Surface_mesh::Vector3(0, 0, 0))
	{
		if (!mesh) return;

		if (!mesh->property("hasNormals").toBool())
		{
			mesh->update_face_normals();
			mesh->update_vertex_normals();
			mesh->setProperty("hasNormals", true);
		}

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_LIGHTING);

		glColorQt(c);

		Surface_mesh::Vertex_property<Surface_mesh::Vector3> points = mesh->vertex_property<Surface_mesh::Vector3>("v:point");
		Surface_mesh::Face_property<Surface_mesh::Vector3> fnormals = mesh->face_property<Surface_mesh::Vector3>("f:normal");

		Surface_mesh::Face_iterator fit, fend = mesh->faces_end();
		Surface_mesh::Vertex_around_face_circulator fvit, fvend;

		glPushMatrix();
		glTranslated(translation[0], translation[1], translation[2]);

		glBegin(GL_TRIANGLES);
		for (fit = mesh->faces_begin(); fit != fend; ++fit){
			glNormal3dv(fnormals[fit].data());
			fvit = fvend = mesh->vertices(fit);
			do{ glVertex3dv(points[fvit].data()); } while (++fvit != fvend);
		}
		glEnd();

		glPopMatrix();
	}

	static void drawMeshWireFrame(SurfaceMeshModel * mesh)
	{
		if (!mesh) return;

		Surface_mesh::Face_iterator fit, fend = mesh->faces_end();
		Surface_mesh::Vertex_around_face_circulator fvit, fvend;
		Surface_mesh::Vertex_property<Surface_mesh::Vector3> points = mesh->vertex_property<Surface_mesh::Vector3>("v:point");
		Surface_mesh::Face_property<Surface_mesh::Vector3> fnormals = mesh->face_property<Surface_mesh::Vector3>("f:normal");

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		glColor4d(0, 1, 1, 0.25);
		glLineWidth(1.0f);
		glEnable(GL_CULL_FACE);
		glCullFace(GL_BACK);

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		for (fit = mesh->faces_begin(); fit != fend; ++fit){
			glBegin(GL_POLYGON);
			glNormal3dv(fnormals[fit].data());
			fvit = fvend = mesh->vertices(fit);
			do{ glVertex3dv(points[fvit].data()); } while (++fvit != fvend);
			glEnd();
		}
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		glDisable(GL_CULL_FACE);
	}

	static void drawMeshName(SurfaceMeshModel * mesh, int name = 0)
	{
		glPushName(name);

		Surface_mesh::Vertex_property<Surface_mesh::Vector3> points = mesh->vertex_property<Surface_mesh::Vector3>("v:point");
		Surface_mesh::Face_iterator fit, fend = mesh->faces_end();
		Surface_mesh::Vertex_around_face_circulator fvit, fvend;

		glBegin(GL_TRIANGLES);
		for (fit = mesh->faces_begin(); fit != fend; ++fit){
			fvit = fvend = mesh->vertices(fit);
			do{ glVertex3dv(points[fvit].data()); } while (++fvit != fvend);
		}
		glEnd();

		glPopName();
	}

	static void drawAABBox(Eigen::AlignedBox3d box, QColor c = QColor())
	{
		GLfloat black[] = { 0.0f, 0.0f, 0.0f, 1.0f };
		GLfloat red[] = { 1.0, 0.3, 0.3, 1.0 };
		GLfloat green[] = { 0.5f, 1.0f, 0.5f, 1.0f };
		GLfloat axis_col[3][4] = { { 1.0f, 0.3f, 0.3f, 1.0f }, { 0.2f, 0.8f, 0.2f, 1.0f }, { 0.1f, 0.1f, 1.0f, 1.0f } };
		GLfloat yellow[] = { 0.9f, 0.9f, 0.0f, 1.0f };

		GLfloat axis[3][3] = { { 1.0, 0, 0 }, { 0, 1.0, 0 }, { 0, 0, 1.0 } };

		Surface_mesh::Vector3 vp[8] = {
			box.corner(Eigen::AlignedBox3d::TopRightCeil),
			box.corner(Eigen::AlignedBox3d::TopRightFloor),
			box.corner(Eigen::AlignedBox3d::BottomRightFloor),
			box.corner(Eigen::AlignedBox3d::BottomRightCeil),
			box.corner(Eigen::AlignedBox3d::TopLeftCeil),
			box.corner(Eigen::AlignedBox3d::TopLeftFloor),
			box.corner(Eigen::AlignedBox3d::BottomLeftFloor),
			box.corner(Eigen::AlignedBox3d::BottomLeftCeil)
		};

		glPushAttrib(GL_LIGHTING_BIT | GL_ENABLE_BIT | GL_HINT_BIT | GL_LINE_BIT | GL_CURRENT_BIT);

		// draw edge
		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glDisable(GL_LIGHTING);

		glColor3fv(black);
		glLineWidth(1.0);

		glBegin(GL_LINE_LOOP);
		glVertex3f(vp[0][0], vp[0][1], vp[0][2]);
		glVertex3f(vp[1][0], vp[1][1], vp[1][2]);
		glVertex3f(vp[2][0], vp[2][1], vp[2][2]);
		glVertex3f(vp[3][0], vp[3][1], vp[3][2]);
		glEnd();
		glBegin(GL_LINE_LOOP);
		glVertex3f(vp[4][0], vp[4][1], vp[4][2]);
		glVertex3f(vp[5][0], vp[5][1], vp[5][2]);
		glVertex3f(vp[6][0], vp[6][1], vp[6][2]);
		glVertex3f(vp[7][0], vp[7][1], vp[7][2]);
		glEnd();
		glBegin(GL_LINES);
		glVertex3f(vp[0][0], vp[0][1], vp[0][2]);
		glVertex3f(vp[4][0], vp[4][1], vp[4][2]);
		glVertex3f(vp[1][0], vp[1][1], vp[1][2]);
		glVertex3f(vp[5][0], vp[5][1], vp[5][2]);
		glVertex3f(vp[2][0], vp[2][1], vp[2][2]);
		glVertex3f(vp[6][0], vp[6][1], vp[6][2]);
		glVertex3f(vp[3][0], vp[3][1], vp[3][2]);
		glVertex3f(vp[7][0], vp[7][1], vp[7][2]);
		glEnd();

		// draw face
		glPushAttrib(GL_DEPTH_BUFFER_BIT);
		glEnable(GL_LIGHTING);
		glEnable(GL_BLEND);
		glDepthMask(GL_FALSE);
		glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
		glEnable(GL_COLOR_MATERIAL);

		if (c!=QColor())
		{
			glColor4f(c.redF(),c.greenF(),c.blueF(),0.2);
		}
		else
		{
			glColor4f(green[0], green[1], green[2], 0.2);
		}

		glBegin(GL_TRIANGLES);
		for (int i = 0; i < 12; i++) {
			glNormal3f(axis[BoxFaceNormalCoeff[i][0]][0] * BoxFaceNormalCoeff[i][1],
				axis[BoxFaceNormalCoeff[i][0]][1] * BoxFaceNormalCoeff[i][1],
				axis[BoxFaceNormalCoeff[i][0]][2] * BoxFaceNormalCoeff[i][1]);
			glVertex3f(vp[BoxFaceVertex[i][0]][0], vp[BoxFaceVertex[i][0]][1], vp[BoxFaceVertex[i][0]][2]);
			glVertex3f(vp[BoxFaceVertex[i][1]][0], vp[BoxFaceVertex[i][1]][1], vp[BoxFaceVertex[i][1]][2]);
			glVertex3f(vp[BoxFaceVertex[i][2]][0], vp[BoxFaceVertex[i][2]][1], vp[BoxFaceVertex[i][2]][2]);
		}
		glEnd();
		glPopAttrib();

		glPopAttrib();
	}

};