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

# AuctionLIB library
LIBS += -L$$PWD/../AuctionLIB/$$CFG/lib -lAuctionLib
INCLUDEPATH += ../AuctionLIB

HEADERS += ddtt-plugin.h ddtt_widget.h Corresponder.h ShapeCorresponder.h \
    DeformPathItem.h \
    DeformPathItemWidget.h \
    DeformScene.h
SOURCES += ddtt-plugin.cpp ddtt_widget.cpp Corresponder.cpp ShapeCorresponder.cpp \
    DeformPathItem.cpp \
    DeformPathItemWidget.cpp \
    DeformScene.cpp

RESOURCES += ddtt-plugin.qrc
FORMS += ddtt_widget.ui
