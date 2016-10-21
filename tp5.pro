QMAKE_MAC_SDK = macosx10.11

QMAKE_CXXFLAGS += -stdlib=libc++ -std=c++11 -fopenmp

INCLUDEPATH += /opt/local/include
INCLUDEPATH += /usr/local/include/libiomp

LIBS += -stdlib=libc++

SOURCES += $$PWD/openglwindow.cpp \
    main.cpp \
    gamewindow.cpp
HEADERS += $$PWD/openglwindow.h \
    gamewindow.h

target.path = .

INSTALLS += target

RESOURCES += gestionnaire.qrc
