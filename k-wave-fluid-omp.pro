#-------------------------------------------------
#
# Project created by QtCreator
#
# k-wave-fluid-omp application
#
#-------------------------------------------------

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG -= debug_and_release
CONFIG += c++11

TARGET = k-wave-fluid-omp

win32:QMAKE_LFLAGS += /ignore:4099

INCLUDEPATH += $$PWD
DEPENDPATH += $$PWD
INCLUDEPATH += $$(MKLROOT)\include\
INCLUDEPATH += $$(MKLROOT)\include\fftw
DEPENDPATH += $$(MKLROOT)\include\fftw

LIBS += -L$$(MKLROOT)\lib\intel64 \
    -lmkl_intel_ilp64 \
    -lmkl_intel_thread \
    -lmkl_core \
    -llibiomp5md \

win32 {
    QMAKE_CXXFLAGS += /DMKL_ILP64
    QMAKE_CFLAGS += /DMKL_ILP64
#    QMAKE_CXXFLAGS += /Qmkl:parallel
#    QMAKE_CFLAGS += /Qmkl:parallel
    DEFINES -= _UNICODE
}

unix {
    QMAKE_CXXFLAGS += -mkl=parallel
    QMAKE_CFLAGS += -mkl=parallel
}

# Check Qt version
lessThan(QT_VERSION, 4.8): error(Qt version is too old)

# HDF5 library
include($$PWD/qtproject/hdf5.pri)

# OpenMP library
include($$PWD/qtproject/openmp.pri)

SOURCES += \
    main.cpp \
    Compression/CompressHelper.cpp \
    Containers/MatrixContainer.cpp \
    Containers/MatrixRecord.cpp \
    Containers/OutputStreamContainer.cpp \
    GetoptWin64/Getopt.cpp \
    Hdf5/Hdf5File.cpp \
    Hdf5/Hdf5FileHeader.cpp \
    KSpaceSolver/KSpaceFirstOrderSolver.cpp \
    Logger/Logger.cpp \
    MatrixClasses/BaseFloatMatrix.cpp \
    MatrixClasses/BaseIndexMatrix.cpp \
    MatrixClasses/ComplexMatrix.cpp \
    MatrixClasses/FftwComplexMatrix.cpp \
    MatrixClasses/IndexMatrix.cpp \
    MatrixClasses/RealMatrix.cpp \
    OutputStreams/BaseOutputStream.cpp \
    OutputStreams/CuboidOutputStream.cpp \
    OutputStreams/IndexOutputStream.cpp \
    OutputStreams/WholeDomainOutputStream.cpp \
    Parameters/CommandLineParameters.cpp \
    Parameters/Parameters.cpp \

HEADERS += \
    Compression/CompressHelper.h \
    Containers/MatrixContainer.h \
    Containers/MatrixRecord.h \
    Containers/OutputStreamContainer.h \
    GetoptWin64/Getopt.h \
    Hdf5/Hdf5File.h \
    Hdf5/Hdf5FileHeader.h \
    KSpaceSolver/KSpaceFirstOrderSolver.h \
    Logger/ErrorMessages.h \
    Logger/ErrorMessagesLinux.h \
    Logger/ErrorMessagesWindows.h \
    Logger/Logger.h \
    Logger/OutputMessages.h \
    Logger/OutputMessagesLinux.h \
    Logger/OutputMessagesWindows.h \
    MatrixClasses/BaseFloatMatrix.h \
    MatrixClasses/BaseIndexMatrix.h \
    MatrixClasses/BaseMatrix.h \
    MatrixClasses/ComplexMatrix.h \
    MatrixClasses/FftwComplexMatrix.h \
    MatrixClasses/IndexMatrix.h \
    MatrixClasses/RealMatrix.h \
    OutputStreams/BaseOutputStream.h \
    OutputStreams/CuboidOutputStream.h \
    OutputStreams/IndexOutputStream.h \
    OutputStreams/WholeDomainOutputStream.h \
    Parameters/CommandLineParameters.h \
    Parameters/Parameters.h \
    Utils/DimensionSizes.h \
    Utils/ErrorMessages.h \
    Utils/MatrixNames.h \
    Utils/TimeMeasure.h \

#QMAKE_POST_LINK += $${QMAKE_COPY} \"$$(ROOT)\redist\intel64_win\mkl\mkl_avx2.dll\" \"$$OUT_PWD/\" &
#QMAKE_POST_LINK += $${QMAKE_COPY} \"$$(ROOT)\redist\intel64_win\mkl\mkl_core.dll\" \"$$OUT_PWD/\" &
#QMAKE_POST_LINK += $${QMAKE_COPY} \"$$(ROOT)\redist\intel64_win\mkl\mkl_def.dll\" \"$$OUT_PWD/\" &
#QMAKE_POST_LINK += $${QMAKE_COPY} \"$$(ROOT)\redist\intel64_win\mkl\mkl_intel_thread.dll\" \"$$OUT_PWD/\"

