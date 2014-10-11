#pragma once

#include "SurfaceMeshPlugins.h"
#include "interfaces/ModePluginDockWidget.h"

#include "StructureGraph.h"
namespace Structure{
    struct Landmark : public Eigen::Vector3d{
        Landmark(size_t id = -1, const Eigen::Vector3d vec = Eigen::Vector3d(0,0,0)) :
            id(id), Eigen::Vector3d(vec){u = v = -1; partid = "none";}
        double u,v;
        QString partid;
        size_t id;
    };
    struct ShapeGraph : public Graph{
        ShapeGraph(QString path):Graph(path){}
        QVector<Landmark> landmarks;
    };
	void saveLandmarks(QString filename){}
	void loadLandmarks(QString filename){}
}

class experiment: public SurfaceMeshModePlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "experiment.plugin.starlab")
    Q_INTERFACES(ModePlugin)

    QIcon icon(){ return QIcon(":/images/experiment.png"); }
	
public:
    QString name() { return "Deformation Correspondence"; }
    QString description() { return "Deformation Correspondence"; }

	bool isApplicable(){ return true; }
	
	void create();
    void destroy(){}
    void decorate();
	
	bool keyPressEvent(QKeyEvent*);
	bool mouseMoveEvent(QMouseEvent*);
	bool mousePressEvent(QMouseEvent* event);
	
    experiment() : widget(NULL), dockwidget(NULL) {}
    QWidget * widget;
	ModePluginDockWidget * dockwidget;

    QVector<Structure::ShapeGraph*> graphs;
};
