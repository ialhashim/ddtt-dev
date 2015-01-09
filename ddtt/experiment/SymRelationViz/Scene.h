#pragma once
#include "MeshModel.h"

class Scene
{
public:
	Scene();
	~Scene();
	void loadScene(const QString dirname);
	void loadConf(const QString filename);
	void setCurrentModelAndPart(int modelNo, int partNo);
	void clearScene();
	void draw();
	void buildModelDislayList();
	void layout();
	QVector<MeshModel*> m_modelList;
	double m_modelScale;
	double m_ptSize;
	QString m_name;
	int m_modelNum; // number of models == m_modelList.size(), but read from the conf.txt before m_modelList is constructed
	int m_currentModelNo;
	int m_currentPartNo;
};

