include($$[STARLAB])
include($$[SURFACEMESH])
include($$[NANOFLANN])
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
            Deformer.h \
            DeformEnergy.h \
            CorrespondenceGenerator.h \
            CorrespondenceSearch.h \
            ShapeGraph.h \
            DeformEnergy2.h \
            EnergyGuidedDeformation.h \
            DeformToFit.h \
    Propagate.h \
    PropagateProximity.h \
    StructureAnalysis.h

SOURCES +=  experiment.cpp \
            experiment-widget.cpp \
            Deformer.cpp \
            DeformEnergy.cpp \
            CorrespondenceGenerator.cpp \
            CorrespondenceSearch.cpp \
            DeformEnergy2.cpp \
            EnergyGuidedDeformation.cpp \
            DeformToFit.cpp \
    Propagate.cpp \
    PropagateProximity.cpp \
    StructureAnalysis.cpp

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
