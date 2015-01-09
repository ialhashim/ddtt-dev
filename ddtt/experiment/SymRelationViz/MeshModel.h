#pragma once
#include "SurfaceMeshModel.h"
using namespace SurfaceMesh;
#include "PointCloud.h"
#include <qgl.h>

class MeshModel
{
public:
	MeshModel();
	~MeshModel();
	void loadObj(const QString &filename);
	void loadParts(const QString &rootPathname, const QString &modelname);
	void setCurrentPartNo(int partNo);
	void loadCorrespondence(const QString& matchPathname, int startNo);
	void parsePairwiseMatches(const QString& filename, QVector<int>& match);
	void normalization(double scale);
	void translate(double x, double y, double z);
	void draw();
	void buildDisplayList(double ptSize=10.0, bool isDrawBBox=false);
	void updateBBoxParts();
private:
	void clear();
	void clearMesh();
	void clearParts();
public:
	GLuint m_displayListID;
	SurfaceMeshModel *m_mesh;
	QVector<PointCloud*> m_parts;
	Eigen::AlignedBox3d m_bboxParts;
	int m_currentPartNo;
	// a m_parts.size()*Scene.m_modelNum matrix, 
	//each column is a set of correspondence between this model and the model represented by the column
	QVector<QVector<int> > m_matches;
	// m_matches[i] is match between this model to the m_matchColIDs[i] model
	// e.g. m_matchColIDs[0] = 1, then m_matches[0] is match between this model with the 1st model
	QVector<int> m_matchColIDs;
	QString m_name;
};

