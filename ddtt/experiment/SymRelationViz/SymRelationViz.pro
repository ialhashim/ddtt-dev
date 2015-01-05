include($$[STARLAB])
include($$[SURFACEMESH])
StarlabTemplate(plugin)

QT += gui opengl

HEADERS +=  SymRelationViz.h \ 
    SymRelationViz-widget.h
SOURCES +=  SymRelationViz.cpp \ 
    SymRelationViz-widget.cpp

FORMS       += SymRelationViz-widget.ui
RESOURCES   += SymRelationViz.qrc
