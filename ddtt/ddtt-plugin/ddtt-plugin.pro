include($$[STARLAB])
include($$[SURFACEMESH])
include($$[OCTREE])
StarlabTemplate(plugin)

QT += gui opengl xml svg

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
LIBS += -L$$PWD/../GlSplatRendererLib/$$CFG/lib -lGlSplatRendererLib
INCLUDEPATH += ../GlSplatRendererLib

# StructureGraph library
LIBS += -L$$PWD/../StructureGraphLib/$$CFG/lib -lStructureGraphLib
INCLUDEPATH += ../StructureGraphLib

# AuctionLIB library
LIBS += -L$$PWD/../AuctionLIB/$$CFG/lib -lAuctionLib
INCLUDEPATH += ../AuctionLIB

HEADERS += ddtt-plugin.h ddtt_widget.h Corresponder.h ShapeCorresponder.h DeformPathItem.h DeformPathItemWidget.h DeformScene.h DeformationPath.h \
    ImageCompare.h \
    PathsGenerator.h
SOURCES += ddtt-plugin.cpp ddtt_widget.cpp Corresponder.cpp ShapeCorresponder.cpp DeformPathItem.cpp DeformPathItemWidget.cpp DeformScene.cpp DeformationPath.cpp \
    ImageCompare.cpp \
    PathsGenerator.cpp
RESOURCES += ddtt-plugin.qrc
FORMS += ddtt_widget.ui
