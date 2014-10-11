#pragma once

#include "tetgen.h"
#include "SurfaceMeshModel.h"

struct TetGen{
    TetGen( SurfaceMesh::SurfaceMeshModel * mesh, QString options = "epq1.414" ){
        //options += QString("a%1").arg( mesh->bbox().diagonal().maxCoeff() * 0.1 );

        tetgenio in;
        in.firstnumber = 0;

        // Push the vertices
        in.pointlist = new REAL[mesh->n_vertices() * 3];
        in.numberofpoints = mesh->n_vertices();
        Vector3VertexProperty points = mesh->vertex_coordinates();
        for(auto v : mesh->vertices()){
            int idx = v.idx();
            for(int i = 0; i < 3; i++){
                in.pointlist[(idx * 3) + i] = points[v][i];
            }
        }

        // Push the faces
        in.facetlist = new tetgenio::facet[mesh->n_faces()];
        in.numberoffacets = mesh->n_faces();
        foreach(Face f, mesh->faces()){
            int idx = f.idx();

            in.facetlist[idx].polygonlist = new tetgenio::polygon[1];
            in.facetlist[idx].numberofpolygons = 1;

            in.facetlist[idx].polygonlist->numberofvertices = 3;
            in.facetlist[idx].polygonlist->vertexlist = new int[3];

            std::vector<int> fverts;
            for(auto v : mesh->vertices(f)) fverts.push_back( v.idx() );

            for(int vi = 0; vi < 3; vi++) in.facetlist[idx].polygonlist->vertexlist[vi] = fverts[vi];

            in.facetlist[idx].numberofholes = 0;
            in.facetlist[idx].holelist = NULL;
        }

        // Execute
        tetgenbehavior tetoptions;
        tetoptions.parse_commandline(options.toLatin1().data());
        tetrahedralize(&tetoptions, &in, &out);
    }

    std::vector<SurfaceMesh::Vector3> getPoints()
    {
        std::vector<SurfaceMesh::Vector3> pnts( out.numberofpoints );
        for(int i = 0; i < out.numberofpoints; i++){
            for(int c = 0; c < 3; c++)
                pnts[i][c] = out.pointlist[(i * 3) + c];
        }
        return pnts;
    }

    std::vector< std::vector<int> > getCells()
    {
        std::vector< std::vector<int> > cells( out.numberoftetrahedra );

        for(int i = 0; i < out.numberoftetrahedra; i++)
        {
            std::vector<int> c_vertices(4);
            for(int j = 0; j < 4; j++) c_vertices[j] = out.tetrahedronlist[i * 4 + j];
            cells[i] = c_vertices;
        }

        return cells;
    }

    tetgenio out;
};
