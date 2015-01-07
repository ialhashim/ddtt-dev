#include "MeshModel.h"
//#include "SurfaceMeshHelper.h"
#include "SurfaceMeshNormalsHelper.h"

#include "utility.h"
#include "QuickMeshDraw.h"

#define _unused(x) ((void)x)

bool modded_read_obj(SurfaceMeshModel& mesh, const std::string& filename){
	char   s[200];
	float  x, y, z;
	float  nx, ny, nz;
	std::vector<Surface_mesh::Vertex>  vertices;

	// open file (in ASCII mode)
	FILE* in = fopen(filename.c_str(), "r");
	if (!in) return false;

	// if mesh is not empty we need an offset for vertex indices
	// also take into accout that OBJ indices start at 1 (not 0)
	const int voffset = mesh.n_vertices() - 1;

	// clear line once
	memset(&s, 0, 200);

	bool hasnormals = false;

	// parse line by line (currently only supports vertex positions & faces
	while (in && !feof(in) && fgets(s, 200, in)) {
		// comment
		if (s[0] == '#' || isspace(s[0])) continue;

		// vertex
		else if (strncmp(s, "v ", 2) == 0) {
			if (sscanf(s, "v %f %f %f", &x, &y, &z)) {
				mesh.add_vertex(Surface_mesh::Point(x, y, z));
			}
		}

		// face
		else if (strncmp(s, "f ", 2) == 0) {
			int component(0), nV(0);
			bool endOfVertex(false);
			char *p0, *p1(s + 1);

			vertices.clear();

			// skip white-spaces
			while (*p1 == ' ') ++p1;

			while (p1) {
				p0 = p1;

				// overwrite next separator

				// skip '/', '\n', ' ', '\0', '\r' <-- don't forget Windows
				while (*p1 != '/' && *p1 != '\r' && *p1 != '\n' && *p1 != ' ' && *p1 != '\0') ++p1;

				// detect end of vertex
				if (*p1 != '/') {
					endOfVertex = true;
				}

				// replace separator by '\0'
				if (*p1 != '\0') {
					*p1 = '\0';
					p1++; // point to next token
				}

				// detect end of line and break
				if (*p1 == '\0' || *p1 == '\n') {
					p1 = 0;
				}

				// read next vertex component
				if (*p0 != '\0') {
					switch (component) {
					case 0: // vertex
						vertices.push_back(Surface_mesh::Vertex(atoi(p0) + voffset));
						break;

					case 1: // texture coord
						break;

					case 2: // normal
						break;
					}
				}

				++component;

				if (endOfVertex) {
					component = 0;
					nV++;
					endOfVertex = false;
				}
			}

			mesh.add_face(vertices);
		}
		else if (!hasnormals && (strncmp(s, "vn ", 3) == 0)){
			hasnormals = true;
		}


		// clear line
		memset(&s, 0, 200);
	}

	/// Normals
	if (hasnormals){
		/// Go back to beginning of file
		rewind(in);
		/// And start reading
		unsigned int ncounter = 0;
		Vector3VertexProperty normals = mesh.vertex_normals(true);
		while (in && !feof(in) && fgets(s, 200, in)) {
			if (strncmp(s, "vn ", 3) == 0) {
				int nread = sscanf(s, "vn %f %f %f", &nx, &ny, &nz);
				assert(nread == 3);
				_unused(nread);
				if (ncounter >= mesh.n_vertices()) continue; // skip duplicated normals
				normals[Vertex(ncounter)] = Vector3(nx, ny, nz);
				// qDebug() << normals[Vertex(ncounter)];                
				ncounter++;
			}
		}

		if (ncounter != mesh.n_vertices())
			qWarning("Warning: while reading file #normals=%d while #vertices=%d", (int)ncounter, mesh.n_vertices());
	}

	fclose(in);
	return true;
}

MeshModel::MeshModel() : m_mesh(0)
{
	
}


MeshModel::~MeshModel()
{
	if (m_mesh) delete m_mesh;
}

void MeshModel::loadObj(QString filename) {
	m_mesh = new SurfaceMeshModel(filename);
	modded_read_obj(*m_mesh, qPrintable(filename)); ///< ~copy/paste from Surface_mesh

	/// If they are not loaded from file, compute normals
	NormalsHelper h(m_mesh);
	if (!m_mesh->has_face_normals())
		h.compute_face_normals();
	if (!m_mesh->has_vertex_normals())
		h.compute_vertex_normals();
}

void MeshModel::buildDisplayList()
{
	if (glIsList(m_displayListID))
	{
		glDeleteLists(m_displayListID, 1);
	}

	m_displayListID = glGenLists(1);

	QColor c(128, 128, 128, 255);
	
	glNewList(m_displayListID, GL_COMPILE);
	QuickMeshDraw::drawMeshSolid(m_mesh, c);
	glEndList();
}
void MeshModel::draw()
{
	glCallList(m_displayListID);
}

void MeshModel::normalization(double scale)
{
	//int nVerts = m_mesh->n_vertices();
	m_mesh->updateBoundingBox();
	Eigen::AlignedBox3d bbox= m_mesh->bbox();
	if (bbox.isNull() || bbox.isEmpty())
		return;

	double len = bbox.diagonal().norm();
	double rs = 2.0 * scale / len;// 
	double x = bbox.center().x();
	double y = bbox.center().y();
	double z = bbox.center().z();
	Surface_mesh::Vertex_property<Surface_mesh::Point> points = m_mesh->vertex_property<Point>("v:point");
	Surface_mesh::Vertex_iterator vit, vend = m_mesh->vertices_end();
	// update vertex coordinates 
	int i(0);
	for (vit = m_mesh->vertices_begin(); vit != vend; ++vit)
	{
		points[vit] = Eigen::Vector3d((points[vit].x()-x) * rs, (points[vit].y()-y) * rs, 
			(points[vit].z()-z) * rs);
		//points[vit] = points[vit];
		//++i;
	}
}
void MeshModel::translate(double x, double y, double z)
{
	Surface_mesh::Vertex_property<Surface_mesh::Point> points = m_mesh->vertex_property<Point>("v:point");
	Surface_mesh::Vertex_iterator vit, vend = m_mesh->vertices_end();
	// update vertex coordinates 
	for (vit = m_mesh->vertices_begin(); vit != vend; ++vit)
	{
		points[vit] = Eigen::Vector3d(points[vit].x() + x, points[vit].y() + y, points[vit].z()+z);
	}
}