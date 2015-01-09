#include "MeshModel.h"
//#include "SurfaceMeshHelper.h"
#include "SurfaceMeshNormalsHelper.h"
#include "utility.h"
#include "QuickMeshDraw.h"
#include "QuickPointsDraw.h"
#include <QDir>



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

MeshModel::MeshModel() : m_currentPartNo(-1),m_mesh(0)
{
	
}
void MeshModel::clearMesh()
{
	if (m_mesh) delete m_mesh;
	m_mesh = 0;
}
void MeshModel::clearParts()
{
	for (QVector<PointCloud*>::iterator it = m_parts.begin(); it != m_parts.end(); ++it)
	{
		delete *it;
	}
	m_parts.clear();
}
void MeshModel::clear()
{
	clearMesh();
	clearParts();
}
MeshModel::~MeshModel()
{
	clear();
}

void MeshModel::loadObj(const QString &filename) 
{
	clearMesh();
	m_mesh = new SurfaceMeshModel(filename);
	modded_read_obj(*m_mesh, qPrintable(filename)); ///< ~copy/paste from Surface_mesh

	/// If they are not loaded from file, compute normals
	NormalsHelper h(m_mesh);
	if (!m_mesh->has_face_normals())
		h.compute_face_normals();
	if (!m_mesh->has_vertex_normals())
		h.compute_vertex_normals();
}
void MeshModel::loadParts(const QString &rootPathname, const QString &modelname)
{
	m_name = modelname;

	QDir dir(rootPathname + m_name);
	if (!dir.exists())
	{
		return;
	}
	clearParts();


	QStringList filters;
	filters << "*.pts";
	dir.setNameFilters(filters);
	//dir.setSorting(QDir::Name);
	QFileInfoList list = dir.entryInfoList();
	m_parts.resize(list.size());
	for (int i = 0; i < list.size(); ++i)
	{
		QFileInfo file_info = list.at(i);
		PointCloud* pts = new PointCloud();
		pts->load(file_info.absoluteFilePath());
		int idx = file_info.baseName().toInt();
		m_parts[idx] = pts;
	}

}
void MeshModel::setCurrentPartNo(int partNo)
{
	m_currentPartNo = partNo;
}
void MeshModel::loadCorrespondence(const QString& matchPathname, int startNo)
{
	QDir dir(matchPathname);
	QString str = QString::number(startNo);
	QStringList filters;
	filters << str + "_*.txt";
	dir.setNameFilters(filters);
	dir.setSorting(QDir::Name);

	m_matches.clear();
	m_matchColIDs.clear();

	QFileInfoList list = dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i)
	{
		if (i == startNo) continue;

		m_matchColIDs.push_back(i);

		QFileInfo file_info = list.at(i);
		QVector<int> match;
		parsePairwiseMatches(file_info.absoluteFilePath(), match);
		m_matches.push_back(match);
	}
}
void MeshModel::parsePairwiseMatches(const QString& filename, QVector<int>& match)
{
	QFile file(filename);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		QString str("can't open ");
		str.append(filename);
		QMessageBox::warning(0, "Warnning", str, QMessageBox::Yes);
		return;
	}
	match.clear();

	QTextStream in(&file);
	int i, j; float p;
	while (!in.atEnd())
	{
		in >> i >> j >> p;
		match.push_back(j);
	}

	file.close();
}
void MeshModel::buildDisplayList(double ptSize, bool isDrawBBox)
{
	if (glIsList(m_displayListID))
	{
		glDeleteLists(m_displayListID, 1);
	}

	m_displayListID = glGenLists(1);

	QColor m(128, 128, 128, 50);
	glNewList(m_displayListID, GL_COMPILE);
	QuickMeshDraw::drawMeshSolid(m_mesh, m);
	

	//////////////
	//if (isDrawBBox)
	//{
	//	QColor b(0, 255, 0, 50);
	//	QuickMeshDraw::drawAABBox(m_mesh->bbox(), b);
	//QuickMeshDraw::drawAABBox(m_bboxParts, b);
	//}

	////////////
	if (m_currentPartNo == -1)
	{
		for (int i = 0; i < m_parts.size(); ++i)
		{
			QuickPointsDraw::drawPoints(m_parts[i], ptSize, qRandomColor());
		}
	}
	else
	{
		QColor p;
		for (int i = 0; i < m_parts.size(); ++i)
		{
			if (m_currentPartNo == i)
			{
				if (isDrawBBox)
					p = QColor(255, 0, 0, 255);
				else
					p = QColor(255, 255, 0, 255);
			}
			else
			{
				p = QColor(125, 125, 125, 50);
			}
			QuickPointsDraw::drawPoints(m_parts[i], ptSize, p);
		}
	}

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

	this->updateBBoxParts();
	x = m_bboxParts.center().x();
	y = m_bboxParts.center().y();
	z = m_bboxParts.center().z();
	double len1 = m_bboxParts.diagonal().norm();
	rs = rs * len / len1;
	for (int i = 0; i < m_parts.size(); ++i)
	{
		PointCloud* pc = m_parts[i];
		for (int j = 0; j < pc->size(); ++j)
		{
			Surface_mesh::Point & pt = pc->points(j);
			pt = Eigen::Vector3d((pt.x() - x) * rs, (pt.y() - y) * rs,
				(pt.z() - z) * rs);
		}
	}

	m_mesh->updateBoundingBox();
	this->updateBBoxParts();
}
void MeshModel::updateBBoxParts()
{
	m_bboxParts.setNull();
	for (int i = 0; i < m_parts.size(); ++i)
	{
		PointCloud* pc = m_parts[i];
		pc->updateBoundingBox();
		m_bboxParts.extend(pc->bbox());
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

	for (int i = 0; i < m_parts.size(); ++i)
	{
		PointCloud* pc = m_parts[i];
		for (int j = 0; j < pc->size(); ++j)
		{
			Surface_mesh::Point & pt = pc->points(j);
			pt = Eigen::Vector3d(pt.x() + x, pt.y() + y, pt.z() + z);
		}
	}

	m_mesh->updateBoundingBox();
	this->updateBBoxParts();
}