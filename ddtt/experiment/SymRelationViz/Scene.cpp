#include "Scene.h"
#include <cmath>

Scene::Scene() : m_modelScale(0.6)
{
}


Scene::~Scene()
{
	clearScene();
}
void Scene::clearScene()
{
	for (QVector<MeshModel*>::iterator it = m_modelList.begin(); it != m_modelList.end(); ++it)
	{
		delete *it;
	}
	m_modelList.clear();
		
}
void Scene::loadScene(const QString dirname)
{
	QDir dir(dirname);
	if (!dir.exists())
	{
		return;
	}
	dir.setFilter(QDir::Files | QDir::NoSymLinks);
	QFileInfoList list = dir.entryInfoList();

	int file_count = list.count();
	if (file_count <= 0)
	{
		return;
	}

	QStringList string_list;
	for (int i = 0; i < list.size(); ++i)
	{
		QFileInfo file_info = list.at(i);
		QString suffix = file_info.suffix();
		if (QString::compare(suffix, QString("obj"), Qt::CaseInsensitive) == 0)
		{
			MeshModel *mm = new MeshModel();
			mm->loadObj(file_info.absoluteFilePath());
			//mm->normalization(m_modelScale);
			mm->m_mesh->updateBoundingBox();
			mm->m_mesh->update_face_normals();
			mm->m_mesh->update_vertex_normals();
			m_modelList.push_back(mm);
		}
	}
}
void Scene::layout()
{
	int ncol = 3;
	int size = m_modelList.size();
	int nrow = std::ceil(size / 3);
	int i(0);
	foreach(MeshModel *m, m_modelList)
	{	
		m->normalization(m_modelScale);
		int row = i / ncol, col = i % ncol;
		m->translate(col,row,0.0);
		++i;
	}
}
void Scene::draw()
{
	foreach(MeshModel *m, m_modelList)
	{
		m->draw();
	}
}
void Scene::buildModelDislayList()
{
	foreach(MeshModel *m, m_modelList)
	{
		m->buildDisplayList();
	}
}
