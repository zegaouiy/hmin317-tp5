QMAKE_MAC_SDK = macosx10.11


INCLUDEPATH += /opt/local/include
INCLUDEPATH += /usr/local/include/libiomp

SOURCES += $$PWD/openglwindow.cpp \
    main.cpp \
    gamewindow.cpp
HEADERS += $$PWD/openglwindow.h \
    gamewindow.h

target.path = .

INSTALLS += target

RESOURCES += gestionnaire.qrc
