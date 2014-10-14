include($$[STARLAB])
include($$[SURFACEMESH])
StarlabTemplate(plugin)

# Build flag
CONFIG(debug, debug|release) {
    CFG = debug
} else {
    CFG = release
}

QT += gui opengl xml svg

HEADERS +=  experiment.h \
            experiment-widget.h \
    Deformer.h
SOURCES +=  experiment.cpp \
            experiment-widget.cpp \
    Deformer.cpp

FORMS       += experiment-widget.ui
RESOURCES   += experiment.qrc

# NURBS library
LIBS += -L$$PWD/../NURBS/$$CFG/lib -lNURBS
INCLUDEPATH += ../NURBS

# StructureGraph library
LIBS += -L$$PWD/../StructureGraphLib/$$CFG/lib -lStructureGraphLib
INCLUDEPATH += ../StructureGraphLib

# ShapeOp
LIBS += -L$$PWD/libShapeOp/$$CFG/lib -llibShapeOp
INCLUDEPATH += ./libShapeOp
DEFINES += SHAPEOP_EXPORT
