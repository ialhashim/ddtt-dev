include($$[STARLAB])
include($$[SURFACEMESH])
StarlabTemplate(plugin)

QT += gui opengl

HEADERS +=  SymRelationViz.h \ 
    SymRelationViz-widget.h \
    Scene.h \
	MeshModel.h \
	QuickMeshDraw.h \
	utility.h
SOURCES +=  SymRelationViz.cpp \ 
    SymRelationViz-widget.cpp \
    Scene.cpp \
	MeshModel.cpp

FORMS       += SymRelationViz-widget.ui
RESOURCES   += SymRelationViz.qrc
