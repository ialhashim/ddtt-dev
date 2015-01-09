include($$[STARLAB])
include($$[SURFACEMESH])
StarlabTemplate(plugin)

QT += gui opengl

HEADERS +=  SymRelationViz.h \ 
    SymRelationViz-widget.h \
    Scene.h \
	MeshModel.h \
	PointCloud.h \
	QuickMeshDraw.h \
	QuickPointsDraw.h \
	utility.h
SOURCES +=  SymRelationViz.cpp \ 
    SymRelationViz-widget.cpp \
    Scene.cpp \
	MeshModel.cpp \
	PointCloud.cpp

FORMS       += SymRelationViz-widget.ui
RESOURCES   += SymRelationViz.qrc
