win32 {
    #QMAKE_CXXFLAGS += -openmp
    #QMAKE_CFLAGS += -openmp
    QMAKE_CXXFLAGS += /openmp
    QMAKE_CFLAGS += /openmp
    DEFINES += _OPENMP
}

unix {
    QMAKE_CXXFLAGS += -fopenmp
    QMAKE_CFLAGS += -fopenmp
    QMAKE_LFLAGS += -fopenmp
    LIBS += -fopenmp
}
