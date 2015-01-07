#pragma once
#include <QDir>

#include "MeshModel.h"

class Scene
{
public:
	Scene();
	~Scene();
	void loadScene(const QString dirname);
	void clearScene();
	void draw();
	void buildModelDislayList();
	void layout();
	QVector<MeshModel*> m_modelList;
	double m_modelScale;
};

