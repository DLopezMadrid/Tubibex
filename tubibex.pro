#-------------------------------------------------
#
# Project created by QtCreator 2015-07-22T23:00:33
#
# Using QT 4.8.6  /  gdb 7.6.1  /  gcc 3.4.4  /  ibex 2.1.16  /  soplex 1.7.2  /  VIBes 0.2  /  dynibex
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = tubibex
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp \
    sivia.cpp \
    init.cpp


#Using ibex 2.1.16 in windows
#IBEXHOME = C:\MinGW\msys\1.0\home\Dani\Ibex
#INCLUDEPATH += $$IBEXHOME/ibex-2.1.16/include/ibex $$IBEXHOME/ibex-2.1.16/include $$IBEXHOME/soplex-1.7.2/src -frounding-math -msse2 -mfpmath=sse
#LIBS += -L$$IBEXHOME/ibex-2.1.16/lib -L$$IBEXHOME/soplex-1.7.2/lib -libex -lsoplex -lprim


#Linux
CONFIG += link_pkgconfig
PKGCONFIG += ibex

HEADERS += \
    sivia.h \
    init.h


#This program uses a modified version of ibex 2.1.16 where a method has been made public instead of private and a modified version of
# dynibex where the run_simulation fn returns a qlist of all the solutions not just the last one
