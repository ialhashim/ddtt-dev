include($$[STARLAB])
include($$[SURFACEMESH])
include($$[OCTREE])
include($$[NANOFLANN])
include($$[QHULL])
StarlabTemplate(plugin)

QT += gui opengl xml svg

HEADERS += \
    particles.h \
    particles-widget.h \
    ParticleMesh.h \
    Particle.h \
    Raytracing.h \
    BasicTable.h \
    StructureAnalysis.h

SOURCES += \
    particles.cpp \
    particles-widget.cpp \
    ParticleMesh.cpp \
    Raytracing.cpp \
    StructureAnalysis.cpp

# SVG viewer
HEADERS += svgview.h
SOURCES += svgview.cpp
	
FORMS       += particles-widget.ui
RESOURCES   += particles.qrc

# Build options
CONFIG(debug, debug|release) {CFG = debug} else {CFG = release}

# Sphere library
LIBS += -L$$PWD/../spherelib/$$CFG/lib -lspherelib
INCLUDEPATH += ../spherelib

# External library
win32:LIBS += -L"$$_PRO_FILE_PWD_/embree2/"

# Qhull linker warning in debug mode
QMAKE_LFLAGS += /ignore:4099
