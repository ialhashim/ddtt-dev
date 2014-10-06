include($$[STARLAB])
StarlabTemplate(none)

QT -= gui

TARGET = libShapeOp
TEMPLATE = lib
CONFIG += staticlib

SOURCES += API.cpp
HEADERS += API.h \
    Common.h \
    Constraint.h \
    Force.h \
    LSSolver.h \
    Solver.h \
    Types.h

# Build options
CONFIG(debug, debug|release) {CFG = debug} else {CFG = release}
DESTDIR = $$PWD/$$CFG/lib

DEFINES += SHAPEOP_EXPORT SHAPEOP_HEADER_ONLY
