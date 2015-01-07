include($$[STARLAB])
include($$[SURFACEMESH])
include($$[NANOFLANN])
include($$[QHULL])
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
            #Deformer.h \
            DeformEnergy.h \
            CorrespondenceGenerator.h \
            CorrespondenceSearch.h \
            ShapeGraph.h \
            DeformEnergy2.h \
            StructureAnalysis.h \
            EnergyGuidedDeformation.h \
            DeformToFit.h \
            PropagateProximity.h \
            PropagateSymmetry.h \
            EvaluateCorrespondence.h \
            BatchProcess.h \
            ComputeApproxMVBB.h \
            convexhull.h

SOURCES +=  experiment.cpp \
            experiment-widget.cpp \
            #Deformer.cpp \
            DeformEnergy.cpp \
            CorrespondenceGenerator.cpp \
            CorrespondenceSearch.cpp \
            DeformEnergy2.cpp \
            StructureAnalysis.cpp \
            EnergyGuidedDeformation.cpp \
            DeformToFit.cpp \
            PropagateProximity.cpp \
            PropagateSymmetry.cpp \
            EvaluateCorrespondence.cpp \
            BatchProcess.cpp \
            ComputeApproxMVBB.cpp

FORMS       += experiment-widget.ui
RESOURCES   += experiment.qrc

# NURBS library
LIBS += -L$$PWD/../NURBS/$$CFG/lib -lNURBS
INCLUDEPATH += ../NURBS

# StructureGraph library
LIBS += -L$$PWD/../StructureGraphLib/$$CFG/lib -lStructureGraphLib
INCLUDEPATH += ../StructureGraphLib

# ShapeOp
#LIBS += -L$$PWD/libShapeOp/$$CFG/lib -llibShapeOp
#INCLUDEPATH += ./libShapeOp
#DEFINES += SHAPEOP_EXPORT

# QDigraph library
QT += webkitwidgets
LIBS += -L$$PWD/libQDigraph/$$CFG/lib -lQDigraphLib
INCLUDEPATH += ./libQDigraph
