include($$[STARLAB])
include($$[SURFACEMESH])
StarlabTemplate(plugin)

Qt += opengl

SOURCES += LFD_driver.cpp 
HEADERS += LFD_driver.h

# Build options
CONFIG(debug, debug|release) {CFG = debug} else {CFG = release}

LIBS += -L$$PWD/$$CFG/lib -llibLFD
INCLUDEPATH += ./
