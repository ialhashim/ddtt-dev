#include "Scene.h"
#include <cmath>
#include <QDir>

Scene::Scene() : m_currentModelNo(0), m_currentPartNo(-1), m_modelScale(0.6), m_ptSize(10.0)
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
void Scene::loadConf(const QString filename)
{
	QFile file(filename);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		QString str("can't open ");
		str.append(filename);
		QMessageBox::warning(0, "Warnning", str, QMessageBox::Yes);
		return;
	}

	QTextStream in(&file);
	m_name = in.readLine();
	m_name = in.readLine(); // 2nd line
	m_modelNum = in.readLine().toInt(); //3rd line

	file.close();
}
void Scene::loadScene(const QString dirname)
{
	QDir dir(dirname);
	if (!dir.exists())
	{
		return;
	}

	////////////////// parse the conf.txt
	QString conf = dir.absolutePath();
	conf.append("/exports/conf.txt");
	loadConf(conf);
	
	////////////////////// parse 
	QStringList filters;
	filters << "*.obj";
	dir.setNameFilters(filters);
	dir.setSorting(QDir::Name);
	QFileInfoList list = dir.entryInfoList();

	QString pathname = dir.absolutePath();
	pathname.append("/exports/");
	for (int i = 0; i < list.size(); ++i)
	{
		QFileInfo file_info = list.at(i);
		MeshModel *mm = new MeshModel();

		// load model
		mm->loadObj(file_info.absoluteFilePath());
		//mm->normalization(m_modelScale);
		mm->m_mesh->updateBoundingBox();
		mm->m_mesh->update_face_normals();
		mm->m_mesh->update_vertex_normals();

		// load parts
		QString modelname = m_name;
		modelname.append("_");
		modelname.append(QString::number(i));
		mm->loadParts(pathname, modelname);

		// load correspondence
		QString matchPathname = pathname;
		matchPathname.append("pairwise_matches/");
		mm->loadCorrespondence(matchPathname,i);

		m_modelList.push_back(mm);
	}

	if (m_modelNum != m_modelList.size())
		QMessageBox::warning(0, "Warnning", "num of models are not equal with the num specified in the conf.txt", QMessageBox::Yes);
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
	MeshModel *mm = this->m_modelList[m_currentModelNo];
	QVector<QVector<int> >::iterator itor = mm->m_matches.begin();
	int i(0);
	foreach(MeshModel *m, m_modelList)
	{
		if (i == m_currentModelNo)
		{
			m->setCurrentPartNo(m_currentPartNo);
		}
		else
		{
			if (m_currentPartNo == -1)
			{
				m->setCurrentPartNo(m_currentPartNo);
			}
			else
			{
				QVector<int> match = *itor;
				m->setCurrentPartNo(match[m_currentPartNo]);
			}
			++itor;
		}
		m->buildDisplayList(m_ptSize, i==m_currentModelNo);
		++i;
	}
}
void Scene::setCurrentModelAndPart(int modelNo, int partNo)
{
	m_currentModelNo = modelNo;
	m_currentPartNo = partNo;
}
