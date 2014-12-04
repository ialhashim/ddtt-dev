include($$[STARLAB])
include($$[SURFACEMESH])
StarlabTemplate(none)

Qt += opengl

TARGET = libLFD
TEMPLATE = lib
CONFIG += staticlib

SOURCES += LFD.cpp LFDWidget.cpp
HEADERS += LFD.h  LFDWidget.h

# Build options
CONFIG(debug, debug|release) {CFG = debug} else {CFG = release}
DESTDIR = $$PWD/$$CFG/lib

FORMS += LFDWidget.ui
