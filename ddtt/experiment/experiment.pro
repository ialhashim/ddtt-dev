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
            experiment-widget.h
SOURCES +=  experiment.cpp \
            experiment-widget.cpp

FORMS       += experiment-widget.ui
RESOURCES   += experiment.qrc

# NURBS library
LIBS += -L$$PWD/../NURBS/$$CFG/lib -lNURBS
INCLUDEPATH += ../NURBS

# StructureGraph library
LIBS += -L$$PWD/../StructureGraphLib/$$CFG/lib -lStructureGraphLib
INCLUDEPATH += ../StructureGraphLib
