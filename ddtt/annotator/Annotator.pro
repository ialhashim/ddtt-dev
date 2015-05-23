include($$[STARLAB])
include($$[SURFACEMESH])
include($$[OCTREE])
StarlabTemplate(appbundle)

QT += core gui opengl svg network

TARGET = Annotator

HEADERS += Annotator.h mydrawarea.h
SOURCES += Annotator.cpp main.cpp  mydrawarea.cpp
FORMS += Annotator.ui

RC_FILE = Annotator.rc

# Build flag
CONFIG(debug, debug|release) {
    CFG = debug
} else {
    CFG = release
}

## Required libraries
# NURBS library
LIBS += -L$$PWD/../NURBS/$$CFG/lib -lNURBS
INCLUDEPATH += ../NURBS

# Surface Reconstruction library
LIBS += -L$$PWD/../Reconstruction/$$CFG/lib -lReconstruction
INCLUDEPATH += ../Reconstruction

# Splat Rendering library
LIBS += -L$$PWD/../GlSplatRendererLib/$$CFG/lib -lGlSplatRendererLib
INCLUDEPATH += ../GlSplatRendererLib

# StructureGraph library
LIBS += -L$$PWD/../StructureGraphLib/$$CFG/lib -lStructureGraphLib
INCLUDEPATH += ../StructureGraphLib

