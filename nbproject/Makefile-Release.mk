#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/HDF5/HDF5_File.o \
	${OBJECTDIR}/KSpaceSolver/KSpaceFirstOrder3DSolver.o \
	${OBJECTDIR}/MatrixClasses/BaseFloatMatrix.o \
	${OBJECTDIR}/MatrixClasses/BaseLongMatrix.o \
	${OBJECTDIR}/MatrixClasses/ComplexMatrix.o \
	${OBJECTDIR}/MatrixClasses/FFTWComplexMatrix.o \
	${OBJECTDIR}/MatrixClasses/LongMatrix.o \
	${OBJECTDIR}/MatrixClasses/MatrixContainer.o \
	${OBJECTDIR}/MatrixClasses/OutputHDF5Stream.o \
	${OBJECTDIR}/MatrixClasses/RealMatrix.o \
	${OBJECTDIR}/MatrixClasses/UXYZ_SGXYZMatrix.o \
	${OBJECTDIR}/Parameters/CommandLineParameters.o \
	${OBJECTDIR}/Parameters/Parameters.o \
	${OBJECTDIR}/main.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-m64 -Wall -fopenmp -msse4.2 -ffast-math -fassociative-math -ftree-vectorizer-verbose=0
CXXFLAGS=-m64 -Wall -fopenmp -msse4.2 -ffast-math -fassociative-math -ftree-vectorizer-verbose=0

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L/home/jiri/Software/hdf5-1.8.10-serial/lib /home/jiri/Software/fftw3/lib/libfftw3f_omp.a /home/jiri/Software/fftw3/lib/libfftw3f.a /home/jiri/Software/hdf5-1.8.10-serial/lib/libhdf5_hl.a /home/jiri/Software/hdf5-1.8.10-serial/lib/libhdf5.a -lz

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/code

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/code: /home/jiri/Software/fftw3/lib/libfftw3f_omp.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/code: /home/jiri/Software/fftw3/lib/libfftw3f.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/code: /home/jiri/Software/hdf5-1.8.10-serial/lib/libhdf5_hl.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/code: /home/jiri/Software/hdf5-1.8.10-serial/lib/libhdf5.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/code: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	g++ -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/code ${OBJECTFILES} ${LDLIBSOPTIONS} -fopenmp -fftw3f -lm

${OBJECTDIR}/HDF5/HDF5_File.o: HDF5/HDF5_File.cpp 
	${MKDIR} -p ${OBJECTDIR}/HDF5
	${RM} $@.d
	$(COMPILE.cc) -O3 -Wall -I/home/jiri/Software/hdf5-1.8.10-serial/include -I/home/jiri/Software/fftw3/include -I. -Wall -fopenmp -msse4.2 -ffast-math -fassociative-math -ftree-vectorizer-verbose=0 -MMD -MP -MF $@.d -o ${OBJECTDIR}/HDF5/HDF5_File.o HDF5/HDF5_File.cpp

${OBJECTDIR}/KSpaceSolver/KSpaceFirstOrder3DSolver.o: KSpaceSolver/KSpaceFirstOrder3DSolver.cpp 
	${MKDIR} -p ${OBJECTDIR}/KSpaceSolver
	${RM} $@.d
	$(COMPILE.cc) -O3 -Wall -I/home/jiri/Software/hdf5-1.8.10-serial/include -I/home/jiri/Software/fftw3/include -I. -Wall -fopenmp -msse4.2 -ffast-math -fassociative-math -ftree-vectorizer-verbose=0 -MMD -MP -MF $@.d -o ${OBJECTDIR}/KSpaceSolver/KSpaceFirstOrder3DSolver.o KSpaceSolver/KSpaceFirstOrder3DSolver.cpp

${OBJECTDIR}/MatrixClasses/BaseFloatMatrix.o: MatrixClasses/BaseFloatMatrix.cpp 
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} $@.d
	$(COMPILE.cc) -O3 -Wall -I/home/jiri/Software/hdf5-1.8.10-serial/include -I/home/jiri/Software/fftw3/include -I. -Wall -fopenmp -msse4.2 -ffast-math -fassociative-math -ftree-vectorizer-verbose=0 -MMD -MP -MF $@.d -o ${OBJECTDIR}/MatrixClasses/BaseFloatMatrix.o MatrixClasses/BaseFloatMatrix.cpp

${OBJECTDIR}/MatrixClasses/BaseLongMatrix.o: MatrixClasses/BaseLongMatrix.cpp 
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} $@.d
	$(COMPILE.cc) -O3 -Wall -I/home/jiri/Software/hdf5-1.8.10-serial/include -I/home/jiri/Software/fftw3/include -I. -Wall -fopenmp -msse4.2 -ffast-math -fassociative-math -ftree-vectorizer-verbose=0 -MMD -MP -MF $@.d -o ${OBJECTDIR}/MatrixClasses/BaseLongMatrix.o MatrixClasses/BaseLongMatrix.cpp

${OBJECTDIR}/MatrixClasses/ComplexMatrix.o: MatrixClasses/ComplexMatrix.cpp 
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} $@.d
	$(COMPILE.cc) -O3 -Wall -I/home/jiri/Software/hdf5-1.8.10-serial/include -I/home/jiri/Software/fftw3/include -I. -Wall -fopenmp -msse4.2 -ffast-math -fassociative-math -ftree-vectorizer-verbose=0 -MMD -MP -MF $@.d -o ${OBJECTDIR}/MatrixClasses/ComplexMatrix.o MatrixClasses/ComplexMatrix.cpp

${OBJECTDIR}/MatrixClasses/FFTWComplexMatrix.o: MatrixClasses/FFTWComplexMatrix.cpp 
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} $@.d
	$(COMPILE.cc) -O3 -Wall -I/home/jiri/Software/hdf5-1.8.10-serial/include -I/home/jiri/Software/fftw3/include -I. -Wall -fopenmp -msse4.2 -ffast-math -fassociative-math -ftree-vectorizer-verbose=0 -MMD -MP -MF $@.d -o ${OBJECTDIR}/MatrixClasses/FFTWComplexMatrix.o MatrixClasses/FFTWComplexMatrix.cpp

${OBJECTDIR}/MatrixClasses/LongMatrix.o: MatrixClasses/LongMatrix.cpp 
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} $@.d
	$(COMPILE.cc) -O3 -Wall -I/home/jiri/Software/hdf5-1.8.10-serial/include -I/home/jiri/Software/fftw3/include -I. -Wall -fopenmp -msse4.2 -ffast-math -fassociative-math -ftree-vectorizer-verbose=0 -MMD -MP -MF $@.d -o ${OBJECTDIR}/MatrixClasses/LongMatrix.o MatrixClasses/LongMatrix.cpp

${OBJECTDIR}/MatrixClasses/MatrixContainer.o: MatrixClasses/MatrixContainer.cpp 
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} $@.d
	$(COMPILE.cc) -O3 -Wall -I/home/jiri/Software/hdf5-1.8.10-serial/include -I/home/jiri/Software/fftw3/include -I. -Wall -fopenmp -msse4.2 -ffast-math -fassociative-math -ftree-vectorizer-verbose=0 -MMD -MP -MF $@.d -o ${OBJECTDIR}/MatrixClasses/MatrixContainer.o MatrixClasses/MatrixContainer.cpp

${OBJECTDIR}/MatrixClasses/OutputHDF5Stream.o: MatrixClasses/OutputHDF5Stream.cpp 
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} $@.d
	$(COMPILE.cc) -O3 -Wall -I/home/jiri/Software/hdf5-1.8.10-serial/include -I/home/jiri/Software/fftw3/include -I. -Wall -fopenmp -msse4.2 -ffast-math -fassociative-math -ftree-vectorizer-verbose=0 -MMD -MP -MF $@.d -o ${OBJECTDIR}/MatrixClasses/OutputHDF5Stream.o MatrixClasses/OutputHDF5Stream.cpp

${OBJECTDIR}/MatrixClasses/RealMatrix.o: MatrixClasses/RealMatrix.cpp 
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} $@.d
	$(COMPILE.cc) -O3 -Wall -I/home/jiri/Software/hdf5-1.8.10-serial/include -I/home/jiri/Software/fftw3/include -I. -Wall -fopenmp -msse4.2 -ffast-math -fassociative-math -ftree-vectorizer-verbose=0 -MMD -MP -MF $@.d -o ${OBJECTDIR}/MatrixClasses/RealMatrix.o MatrixClasses/RealMatrix.cpp

${OBJECTDIR}/MatrixClasses/UXYZ_SGXYZMatrix.o: MatrixClasses/UXYZ_SGXYZMatrix.cpp 
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} $@.d
	$(COMPILE.cc) -O3 -Wall -I/home/jiri/Software/hdf5-1.8.10-serial/include -I/home/jiri/Software/fftw3/include -I. -Wall -fopenmp -msse4.2 -ffast-math -fassociative-math -ftree-vectorizer-verbose=0 -MMD -MP -MF $@.d -o ${OBJECTDIR}/MatrixClasses/UXYZ_SGXYZMatrix.o MatrixClasses/UXYZ_SGXYZMatrix.cpp

${OBJECTDIR}/Parameters/CommandLineParameters.o: Parameters/CommandLineParameters.cpp 
	${MKDIR} -p ${OBJECTDIR}/Parameters
	${RM} $@.d
	$(COMPILE.cc) -O3 -Wall -I/home/jiri/Software/hdf5-1.8.10-serial/include -I/home/jiri/Software/fftw3/include -I. -Wall -fopenmp -msse4.2 -ffast-math -fassociative-math -ftree-vectorizer-verbose=0 -MMD -MP -MF $@.d -o ${OBJECTDIR}/Parameters/CommandLineParameters.o Parameters/CommandLineParameters.cpp

${OBJECTDIR}/Parameters/Parameters.o: Parameters/Parameters.cpp 
	${MKDIR} -p ${OBJECTDIR}/Parameters
	${RM} $@.d
	$(COMPILE.cc) -O3 -Wall -I/home/jiri/Software/hdf5-1.8.10-serial/include -I/home/jiri/Software/fftw3/include -I. -Wall -fopenmp -msse4.2 -ffast-math -fassociative-math -ftree-vectorizer-verbose=0 -MMD -MP -MF $@.d -o ${OBJECTDIR}/Parameters/Parameters.o Parameters/Parameters.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O3 -Wall -I/home/jiri/Software/hdf5-1.8.10-serial/include -I/home/jiri/Software/fftw3/include -I. -Wall -fopenmp -msse4.2 -ffast-math -fassociative-math -ftree-vectorizer-verbose=0 -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/code

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
