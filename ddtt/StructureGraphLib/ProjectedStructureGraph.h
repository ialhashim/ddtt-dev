#pragma once
#include "StructureGraph.h"
#include "GraphCorresponder.h"

#include "SoftwareRenderer.h"
#include "MarchingSquares.h"

// Triangulation library
#include "polypartition.hpp"

#include "glhelper.h"

struct NodeProjection{
	Array2D_Vector3 skeletonFrames;
	Array1D_Vector3 boundary;
};

class ProjectedStructureGraph
{
public:
	ProjectedStructureGraph(Structure::Graph *inputGraph, int screenWidth) : graph(inputGraph), screenWidth(screenWidth)
	{
		if(!inputGraph) return;

		// Setup projection and camera
		updateCamera();

		// Encode projected geometries
		viewport = Eigen::Vector4i(0, 0, screenWidth, screenWidth);
		mvp = projectionMatrix * cameraMatrix;

		// Encode
		projectNodes();
	}
	
    Structure::Graph * graph;

	// Projection parameters
	int screenWidth;
    Eigen::Matrix4d projectionMatrix, cameraMatrix;
	Eigen::Matrix4d mvp;
	Eigen::Vector4i viewport;

	// Node projections
	QMap<QString, NodeProjection> projections;
	void projectNodes()
	{
		bool isBuildSkeletonFrames = false;

		projections.clear();

		for(auto n : graph->nodes)
		{
			Array2D_Vector3 skeletonFrames;
			Array1D_Vector3 boundary;
		
			// Build skeleton frames
			if( isBuildSkeletonFrames )
			{
				// Sample node skeleton
				double res = n->bbox().sizes().norm() * 0.1;
				Array1D_Vector3 skeletonPoints = n->discretizedAsCurve( res );

				// Project skeleton to screen space
				for(auto & p : skeletonPoints)
				{
					p = (mvp * p.homogeneous()).colwise().hnormalized();
					p[0] = viewport[0] + viewport[2] * (p[0]+1)/2;
					p[1] = viewport[1] + viewport[3] * (-p[1]+1)/2; // notice '-y' to get from top left
					p[2] = p[2];
				}

				// Re-sample to 'N' points
				int N = 10;
				NURBS::NURBSCurved::createCurveFromPoints(skeletonPoints).SubdivideByLength(N, skeletonPoints);

				for(size_t i = 0; i + 1 < skeletonPoints.size(); i++)
				{
					Vector3 p0 = skeletonPoints[i], p1 = skeletonPoints[i+1];

					Array1D_Vector3 frame(3, Vector3(0,0,0));
					frame[0] = p0; // position
					frame[1] = (p1-p0).normalized(); // tangent
					frame[2] = rotatedVec(frame[1], M_PI_2, Vector3(0,0,1)); // normal

					skeletonFrames.push_back( frame );

					// Last frame
					if(i + 2 == skeletonPoints.size()){
						skeletonFrames.push_back(frame);
						skeletonFrames.back().front() = skeletonPoints.back();
					}
				}
			}

			// Render geometry as one solid blob to a 'buffer'
			if( true )
			{
				QVector< QVector<Vector3> > triangles;
				auto nodeMesh = n->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >();
				for(auto f : nodeMesh->faces()){
					QVector<Vector3> face;
					for(auto v : nodeMesh->vertices(f)){
						Vector3 p = nodeMesh->vertex_coordinates()[v];
						p = (mvp * p.homogeneous()).colwise().hnormalized();
						p[0] = viewport[0] + viewport[2] * (p[0]+1)/2;
						p[1] = viewport[1] + viewport[3] * (-p[1]+1)/2; // notice '-y' to get from top left
						face.push_back( p );
					}
					triangles.push_back( face );
				}
				MatrixXd buffer = SoftwareRenderer::renderTriangles2D(triangles, viewport[2], viewport[3]);
				for(auto p : MarchingSquares::march(buffer, 1.0))
					boundary.push_back( Vector3(p.x(), p.y(), 0) );
			}

			// Experiment: fixed number points + top corner start
			{
				boundary = refineByNumber(boundary, 300);

				size_t closestStart = 0;
				double minDist = DBL_MAX;

				for(size_t i = 0; i < boundary.size(); i++){
					double dist = boundary[i].norm();
					if(dist < minDist){
						minDist = dist;
						closestStart = i;
					}
				}

				std::rotate(boundary.begin(), boundary.begin() + closestStart, boundary.end());

				// Fix orientation if needed
				TPPLPoly poly; poly.Init((long)boundary.size());
				for (int i = 0; i < (int)boundary.size(); i++){
					poly[i].x = boundary[i].x(); poly[i].y = boundary[i].y();
				}
				if (poly.GetOrientation() == TPPL_CW)
					std::reverse(boundary.begin(), boundary.end());
			}

			// Encode boundary
			{

			}

			// Debug skeleton frames
			/*
			{
				QPainter painter(&myimg);
				painter.setRenderHint(QPainter::Antialiasing);

				for(size_t i = 0; i + 1 < skeletonFrames.size(); i++)
				{
					double length = 40;

					Vector3 p0 = skeletonFrames[i][0], p1 = skeletonFrames[i+1][0];
					Vector3 p2 = p0 + skeletonFrames[i][2] * length;

					drawArrow(p0, p1, Qt::red, painter);
					drawArrow(p0, p2, Qt::green, painter);

					// Last frame
					if(i + 2 == skeletonPoints.size()){
						p0 = p1;
						Vector3 p2 = p0 + skeletonFrames[i+1][2] * length;
						drawArrow(p0, p2, Qt::green, painter);
					}
				}
			}*/

			NodeProjection nproj;
			nproj.boundary = boundary;
			nproj.skeletonFrames = skeletonFrames;
			projections[n->id] = nproj;
		}
	}

	void updateCamera(double rotationAngle = M_PI * 1.3)
	{
		Eigen::Vector3d center, eye;
		Eigen::Vector3d dir( cos(rotationAngle), sin(rotationAngle), 0.25);
		Eigen::AlignedBox3d bbox = graph->bbox();
		center = bbox.center();
    
		//double radius = bbox.sizes().norm() * 3;
		double radius = 4.5; // Fixed

		eye = center + (dir.normalized() * radius);
		projectionMatrix = perspective<double>(20, 1.0, 0.01, 1000);
		cameraMatrix = lookAt< Eigen::Vector3d >(eye, center, Eigen::Vector3d(0,0,1));
	}

	void drawBlended(ProjectedStructureGraph * pgOther, GraphCorresponder * gcorr, double alpha)
	{
		if(!graph || graph->nodes.empty() || !gcorr || !pgOther) return;

		QColor fillColor(0,0,0);
		QColor borderColor(255,0,255);
		bool isDrawBorder = false;
		bool isTriangulate = false;

		glDisable(GL_LIGHTING);
		glLineWidth( 3.0 );

		glDisable(GL_DEPTH_TEST);

		// 2D mode
		{
			Eigen::Vector4i viewport;
			glGetIntegerv(GL_VIEWPORT, viewport.data());
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			glOrtho(0.0f, viewport[2], viewport[3], 0.0f, 0.0f, 1.0f);
			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
		}

		for(auto node : graph->nodes)
		{
			std::vector<QString> correspond = gcorr->correspondingNodesTarget(node->id);
			if(correspond.empty()) continue;

			QString tid = correspond.front(); // for now..

			if(!projections.contains(node->id)) continue;
			Array1D_Vector3 src = projections[node->id].boundary;
			Array1D_Vector3 tgt = pgOther->projections[tid].boundary;

			Array1D_Vector3 blended;
			for(size_t i = 0; i < src.size(); i++)
			{
				Vector3 p = AlphaBlend( qRanged(0.0, alpha, 1.0) , src[i], tgt[i]);
				blended.push_back(p);
			}

			if( isTriangulate )
			{
				TPPLPoly poly; poly.Init((long)blended.size());
				for(int i = 0; i < (int)blended.size(); i++){
					poly[i].x = blended[i].x();
					poly[i].y = blended[i].y();
				}
				TPPLPartition pp;
				std::list<TPPLPoly> inpolys, outpolys;
				inpolys.push_back(poly);

				pp.Triangulate_MONO(&inpolys, &outpolys);

				glColorQt(QColor(200,0,200));
				glBegin(GL_TRIANGLES);
				for(auto poly : outpolys)
					for(int i = 0; i < poly.GetNumPoints(); i++)
						glVertex2d(poly[i].x, poly[i].y);
				glEnd();
			}
			else
			{
				glClear(GL_STENCIL_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

				glEnable(GL_STENCIL_TEST);
				glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE);
				glDepthMask(GL_FALSE);
				glEnable(GL_STENCIL_TEST);
				glStencilFunc(GL_ALWAYS,0x1,0x1);
				glStencilOp(GL_KEEP,GL_KEEP,GL_INVERT);
				glBegin(GL_TRIANGLE_FAN);
				for(size_t i = 0; i < blended.size(); i++)
				{
					Vector3 p = blended[i];
					glVector3(p);
				}
				glEnd();


				glColorQt( fillColor );

				glDepthMask(GL_TRUE);
				glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
				glStencilFunc(GL_EQUAL,0x1,0x1);
				glStencilOp(GL_KEEP,GL_KEEP,GL_KEEP);
				glBegin(GL_TRIANGLE_FAN);
				for(size_t i = 0; i < blended.size(); i++){
					Vector3 p = blended[i];
					glVector3(p);
				}
				glEnd();
			
				glDisable(GL_STENCIL_TEST);
			}

			// Border
			if( isDrawBorder )
			{
				glColorQt( borderColor );
				glBegin(GL_LINE_LOOP);
				for(size_t i = 0; i < blended.size(); i++)
				{
					Vector3 p = blended[i];
					glVector3(p);
				}
				glEnd();
			}
		}
		glEnable(GL_DEPTH_TEST);
	}

	QImage drawBlendedImage(ProjectedStructureGraph * pgOther, GraphCorresponder * gcorr, double alpha)
	{
		QImage img(screenWidth, screenWidth, QImage::Format_RGB888);
		img.fill( qRgba(255,255,255,255) ); // white background
		if(!graph || !pgOther || !pgOther->graph || !gcorr) return img;

		// Draw blended boundaries
		{
			QPainter painter( &img );
			painter.setBrush(QBrush(Qt::black));
			painter.setPen(Qt::NoPen);

			for(auto node : graph->nodes)
			{
				std::vector<QString> correspond = gcorr->correspondingNodesTarget(node->id);
				if(correspond.empty()) continue;

				QString tid = correspond.front(); // for now..

				if(!projections.contains(node->id)) continue;
				Array1D_Vector3 src = projections.value(node->id).boundary;
				Array1D_Vector3 tgt = pgOther->projections.value(tid).boundary;

				QPainterPath path;

				for(size_t i = 0; i < src.size(); i++)
				{
					Vector3 p = AlphaBlend( qRanged(0.0, alpha, 1.0) , src[i], tgt[i]);
					if(i == 0) 
						path.moveTo(p.x(), p.y());
					else
						path.lineTo(p.x(), p.y());
				}

				painter.drawPath( path );
			}
		}

		return img;
	}

	QImage drawBoundaryImage()
	{
		QImage img(screenWidth, screenWidth, QImage::Format_RGB888);
		img.fill( qRgba(255,255,255,255) ); // white background
		if(!graph) return img;

		QPainter painter( &img );
		painter.setBrush(QBrush(Qt::black));
		painter.setPen(Qt::NoPen);

		for(auto node : graph->nodes)
		{
			if(!projections.contains(node->id)) continue;
			Array1D_Vector3 src = projections.value(node->id).boundary;
			
			QPainterPath path;

			for(size_t i = 0; i < src.size(); i++)
			{
				Vector3 p = src[i];

				if(i == 0) path.moveTo(p.x(), p.y());
				else path.lineTo(p.x(), p.y());
			}

			painter.drawPath( path );
		}

		return img;
	}
};
