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

# NURBS library
LIBS += -L$$PWD/../NURBS/$$CFG/lib -lNURBS
INCLUDEPATH += ../NURBS

# StructureGraph library
LIBS += -L$$PWD/../StructureGraphLib/$$CFG/lib -lStructureGraphLib
INCLUDEPATH += ../StructureGraphLib

HEADERS += ddtt-plugin.h ddtt_widget.h Corresponder.h ShapeCorresponder.h
SOURCES += ddtt-plugin.cpp ddtt_widget.cpp Corresponder.cpp ShapeCorresponder.cpp

RESOURCES += ddtt-plugin.qrc
FORMS += ddtt_widget.ui
