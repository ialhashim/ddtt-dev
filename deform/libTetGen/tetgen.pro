include($$[STARLAB])
include($$[SURFACEMESH])
StarlabTemplate(none)

TEMPLATE = lib
CONFIG += staticlib

## Build flag
CONFIG(debug, debug|release) {
    CFG = debug
} else {
    CFG = release
}

## Library name and destination
TARGET = libTetGen
DESTDIR = $$PWD/$$CFG/lib

HEADERS += tetgen.h \
    tetgenLib.h
SOURCES += tetgen.cxx predicates.cxx

