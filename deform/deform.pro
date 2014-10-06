include($$[STARLAB])
include($$[SURFACEMESH])
StarlabTemplate(plugin)

QT += gui opengl xml svg

HEADERS += \
    deform.h \
    deform-widget.h

SOURCES += \
    deform.cpp \
    deform-widget.cpp

FORMS       += deform-widget.ui
RESOURCES   += deform.qrc

# Build options
CONFIG(debug, debug|release) {CFG = debug} else {CFG = release}

## Libraries:
# ShapeOp
INCLUDEDIR += libShapeOp
DEFINES += SHAPEOP_HEADER_ONLY

