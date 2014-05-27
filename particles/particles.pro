include($$[STARLAB])
include($$[SURFACEMESH])
include($$[OCTREE])
include($$[NANOFLANN])
StarlabTemplate(plugin)

QT += gui opengl xml svg

HEADERS += \
    particles.h \
    particles-widget.h \
    ParticleMesh.h \
    Particle.h

SOURCES += \
    particles.cpp \
    particles-widget.cpp \
    ParticleMesh.cpp \
    Particle.cpp

FORMS       += particles-widget.ui
RESOURCES   += particles.qrc
