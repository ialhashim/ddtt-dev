include($$[STARLAB])
include($$[SURFACEMESH])

include($$[NANOFLANN])
StarlabTemplate(app)

#CONFIG *= console
TARGET = objToSkeleton

# Build flag
CONFIG(debug, debug|release) {
    CFG = debug
} else {
    CFG = release
}

# NURBS library
LIBS += -L$$PWD/../NURBS/$$CFG/lib -lNURBS
INCLUDEPATH += ../NURBS

# Surface Reconstruction library
LIBS += -L$$PWD/../Reconstruction/$$CFG/lib -lReconstruction
INCLUDEPATH += ../Reconstruction

# Splat Rendering library
LIBS += -L$$PWD/../GLSplatRendererLib/$$CFG/lib -lGLSplatRendererLib
INCLUDEPATH += ../GLSplatRendererLib

# TopoBlender library
LIBS += -L$$PWD/../StructureGraphLib/$$CFG/lib -lStructureGraphLib
INCLUDEPATH += ../StructureGraphLib

SOURCES +=  objToSkeleton.cpp
HEADERS += objToSkeleton.h

mac:LIBS += -framework CoreFoundation # We need this for GLee

DESTDIR = $$PWD/$$CFG
