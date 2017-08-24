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
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Debug-pcjaros-static
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/Hdf5/Hdf5File.o \
	${OBJECTDIR}/KSpaceSolver/KSpaceFirstOrder3DSolver.o \
	${OBJECTDIR}/MatrixClasses/BaseFloatMatrix.o \
	${OBJECTDIR}/MatrixClasses/BaseIndexMatrix.o \
	${OBJECTDIR}/MatrixClasses/ComplexMatrix.o \
	${OBJECTDIR}/MatrixClasses/FFTWComplexMatrix.o \
	${OBJECTDIR}/MatrixClasses/IndexMatrix.o \
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
CCFLAGS=-m64 -std=c++11 -O3 -g3 -gdwarf-2 -Wall -fopenmp -mavx -mtune=native -march=native -ffast-math -fassociative-math
CXXFLAGS=-m64 -std=c++11 -O3 -g3 -gdwarf-2 -Wall -fopenmp -mavx -mtune=native -march=native -ffast-math -fassociative-math

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L./ -L${EBROOTHDF5}/lib -L${EBROOTFFTW}/lib -Wl,-rpath,'./' -Wl,-rpath,'${EBROOTHDF5}/lib' -Wl,-rpath,'${EBROOTFFTW}/lib' ${EBROOTHDF5}/lib/libhdf5_hl.a ${EBROOTHDF5}/lib/libhdf5.a ${EBROOTSZIP}/lib/libsz.a ${EBROOTZLIB}/lib/libz.a ${EBROOTFFTW}/lib/libfftw3f.a ${EBROOTFFTW}/lib/libfftw3f_omp.a -ldl

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/k-wave-fluid-omp

${CND_DISTDIR}/${CND_CONF}/k-wave-fluid-omp: ${EBROOTHDF5}/lib/libhdf5_hl.a

${CND_DISTDIR}/${CND_CONF}/k-wave-fluid-omp: ${EBROOTHDF5}/lib/libhdf5.a

${CND_DISTDIR}/${CND_CONF}/k-wave-fluid-omp: ${EBROOTSZIP}/lib/libsz.a

${CND_DISTDIR}/${CND_CONF}/k-wave-fluid-omp: ${EBROOTZLIB}/lib/libz.a

${CND_DISTDIR}/${CND_CONF}/k-wave-fluid-omp: ${EBROOTFFTW}/lib/libfftw3f.a

${CND_DISTDIR}/${CND_CONF}/k-wave-fluid-omp: ${EBROOTFFTW}/lib/libfftw3f_omp.a

${CND_DISTDIR}/${CND_CONF}/k-wave-fluid-omp: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}
	g++ -o ${CND_DISTDIR}/${CND_CONF}/k-wave-fluid-omp ${OBJECTFILES} ${LDLIBSOPTIONS} -O3 -g3 -gdwarf-2 -Wall -fopenmp -mavx -mtune=native -march=native

${OBJECTDIR}/Hdf5/Hdf5File.o: Hdf5/Hdf5File.cpp
	${MKDIR} -p ${OBJECTDIR}/Hdf5
	${RM} "$@.d"
	$(COMPILE.cc) -w -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Hdf5/Hdf5File.o Hdf5/Hdf5File.cpp

${OBJECTDIR}/KSpaceSolver/KSpaceFirstOrder3DSolver.o: KSpaceSolver/KSpaceFirstOrder3DSolver.cpp
	${MKDIR} -p ${OBJECTDIR}/KSpaceSolver
	${RM} "$@.d"
	$(COMPILE.cc) -w -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/KSpaceSolver/KSpaceFirstOrder3DSolver.o KSpaceSolver/KSpaceFirstOrder3DSolver.cpp

${OBJECTDIR}/MatrixClasses/BaseFloatMatrix.o: MatrixClasses/BaseFloatMatrix.cpp
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} "$@.d"
	$(COMPILE.cc) -w -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MatrixClasses/BaseFloatMatrix.o MatrixClasses/BaseFloatMatrix.cpp

${OBJECTDIR}/MatrixClasses/BaseIndexMatrix.o: MatrixClasses/BaseIndexMatrix.cpp
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} "$@.d"
	$(COMPILE.cc) -w -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MatrixClasses/BaseIndexMatrix.o MatrixClasses/BaseIndexMatrix.cpp

${OBJECTDIR}/MatrixClasses/ComplexMatrix.o: MatrixClasses/ComplexMatrix.cpp
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} "$@.d"
	$(COMPILE.cc) -w -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MatrixClasses/ComplexMatrix.o MatrixClasses/ComplexMatrix.cpp

${OBJECTDIR}/MatrixClasses/FFTWComplexMatrix.o: MatrixClasses/FFTWComplexMatrix.cpp
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} "$@.d"
	$(COMPILE.cc) -w -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MatrixClasses/FFTWComplexMatrix.o MatrixClasses/FFTWComplexMatrix.cpp

${OBJECTDIR}/MatrixClasses/IndexMatrix.o: MatrixClasses/IndexMatrix.cpp
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} "$@.d"
	$(COMPILE.cc) -w -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MatrixClasses/IndexMatrix.o MatrixClasses/IndexMatrix.cpp

${OBJECTDIR}/MatrixClasses/MatrixContainer.o: MatrixClasses/MatrixContainer.cpp
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} "$@.d"
	$(COMPILE.cc) -w -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MatrixClasses/MatrixContainer.o MatrixClasses/MatrixContainer.cpp

${OBJECTDIR}/MatrixClasses/OutputHDF5Stream.o: MatrixClasses/OutputHDF5Stream.cpp
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} "$@.d"
	$(COMPILE.cc) -w -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MatrixClasses/OutputHDF5Stream.o MatrixClasses/OutputHDF5Stream.cpp

${OBJECTDIR}/MatrixClasses/RealMatrix.o: MatrixClasses/RealMatrix.cpp
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} "$@.d"
	$(COMPILE.cc) -w -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MatrixClasses/RealMatrix.o MatrixClasses/RealMatrix.cpp

${OBJECTDIR}/MatrixClasses/UXYZ_SGXYZMatrix.o: MatrixClasses/UXYZ_SGXYZMatrix.cpp
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} "$@.d"
	$(COMPILE.cc) -w -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MatrixClasses/UXYZ_SGXYZMatrix.o MatrixClasses/UXYZ_SGXYZMatrix.cpp

${OBJECTDIR}/Parameters/CommandLineParameters.o: Parameters/CommandLineParameters.cpp
	${MKDIR} -p ${OBJECTDIR}/Parameters
	${RM} "$@.d"
	$(COMPILE.cc) -w -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Parameters/CommandLineParameters.o Parameters/CommandLineParameters.cpp

${OBJECTDIR}/Parameters/Parameters.o: Parameters/Parameters.cpp
	${MKDIR} -p ${OBJECTDIR}/Parameters
	${RM} "$@.d"
	$(COMPILE.cc) -w -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Parameters/Parameters.o Parameters/Parameters.cpp

${OBJECTDIR}/main.o: main.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -w -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
