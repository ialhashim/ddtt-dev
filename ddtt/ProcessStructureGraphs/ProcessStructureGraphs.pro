#-------------------------------------------------
#
# Project created by QtCreator 2015-05-13T11:07:39
#
#-------------------------------------------------
include($$[STARLAB])
include($$[SURFACEMESH])
StarlabTemplate(appbundle)

QT       += core gui opengl xml

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = ProcessStructureGraphs
TEMPLATE = app


SOURCES += main.cpp\
           mainwindow.cpp \
           Evaluator.cpp

HEADERS  += mainwindow.h \
            globals.h \
            Evaluator.h

FORMS    += mainwindow.ui \
            Evaluator.ui

CONFIG(debug, debug|release) {
    CFG = debug
} else {
    CFG = release
}

# NURBS library
LIBS += -L$$PWD/../NURBS/$$CFG/lib -lNURBS
INCLUDEPATH += $$PWD/../NURBS

# StructureGraph library
LIBS += -L$$PWD/../StructureGraphLib/$$CFG/lib -lStructureGraphLib
INCLUDEPATH += $$PWD/../StructureGraphLib
