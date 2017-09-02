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
CND_CONF=Debug-generic-dynamic
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/Containers/MatrixContainer.o \
	${OBJECTDIR}/Containers/MatrixRecord.o \
	${OBJECTDIR}/Containers/OutputStreamContainer.o \
	${OBJECTDIR}/Hdf5/Hdf5File.o \
	${OBJECTDIR}/Hdf5/Hdf5FileHeader.o \
	${OBJECTDIR}/KSpaceSolver/KSpaceFirstOrder3DSolver.o \
	${OBJECTDIR}/Logger/Logger.o \
	${OBJECTDIR}/MatrixClasses/BaseFloatMatrix.o \
	${OBJECTDIR}/MatrixClasses/BaseIndexMatrix.o \
	${OBJECTDIR}/MatrixClasses/ComplexMatrix.o \
	${OBJECTDIR}/MatrixClasses/FftwComplexMatrix.o \
	${OBJECTDIR}/MatrixClasses/IndexMatrix.o \
	${OBJECTDIR}/MatrixClasses/RealMatrix.o \
	${OBJECTDIR}/OutputStreams/BaseOutputStream.o \
	${OBJECTDIR}/OutputStreams/CuboidOutputStream.o \
	${OBJECTDIR}/OutputStreams/IndexOutputStream.o \
	${OBJECTDIR}/OutputStreams/WholeDomainOutputStream.o \
	${OBJECTDIR}/Parameters/CommandLineParameters.o \
	${OBJECTDIR}/Parameters/Parameters.o \
	${OBJECTDIR}/main.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-m64 -std=c++11 -O3 -g3 -gdwarf-2 -Wall -fopenmp
CXXFLAGS=-m64 -std=c++11 -O3 -g3 -gdwarf-2 -Wall -fopenmp

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L./ -L${EBROOTHDF5}/lib -L${EBROOTFFTW}/lib -Wl,-rpath,'./' -Wl,-rpath,'${EBROOTHDF5}/lib' -Wl,-rpath,'${EBROOTFFTW}/lib' -lfftw3f -lfftw3f_omp -lhdf5 -lhdf5_hl -lsz -lz

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/k-wave-fluid-omp

${CND_DISTDIR}/${CND_CONF}/k-wave-fluid-omp: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}
	g++ -o ${CND_DISTDIR}/${CND_CONF}/k-wave-fluid-omp ${OBJECTFILES} ${LDLIBSOPTIONS} -O3 -g3 -gdwarf-2 -Wall -fopenmp

${OBJECTDIR}/Containers/MatrixContainer.o: Containers/MatrixContainer.cpp
	${MKDIR} -p ${OBJECTDIR}/Containers
	${RM} "$@.d"
	$(COMPILE.cc) -Wall -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Containers/MatrixContainer.o Containers/MatrixContainer.cpp

${OBJECTDIR}/Containers/MatrixRecord.o: Containers/MatrixRecord.cpp
	${MKDIR} -p ${OBJECTDIR}/Containers
	${RM} "$@.d"
	$(COMPILE.cc) -Wall -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Containers/MatrixRecord.o Containers/MatrixRecord.cpp

${OBJECTDIR}/Containers/OutputStreamContainer.o: Containers/OutputStreamContainer.cpp
	${MKDIR} -p ${OBJECTDIR}/Containers
	${RM} "$@.d"
	$(COMPILE.cc) -Wall -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Containers/OutputStreamContainer.o Containers/OutputStreamContainer.cpp

${OBJECTDIR}/Hdf5/Hdf5File.o: Hdf5/Hdf5File.cpp
	${MKDIR} -p ${OBJECTDIR}/Hdf5
	${RM} "$@.d"
	$(COMPILE.cc) -Wall -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Hdf5/Hdf5File.o Hdf5/Hdf5File.cpp

${OBJECTDIR}/Hdf5/Hdf5FileHeader.o: Hdf5/Hdf5FileHeader.cpp
	${MKDIR} -p ${OBJECTDIR}/Hdf5
	${RM} "$@.d"
	$(COMPILE.cc) -Wall -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Hdf5/Hdf5FileHeader.o Hdf5/Hdf5FileHeader.cpp

${OBJECTDIR}/KSpaceSolver/KSpaceFirstOrder3DSolver.o: KSpaceSolver/KSpaceFirstOrder3DSolver.cpp
	${MKDIR} -p ${OBJECTDIR}/KSpaceSolver
	${RM} "$@.d"
	$(COMPILE.cc) -Wall -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/KSpaceSolver/KSpaceFirstOrder3DSolver.o KSpaceSolver/KSpaceFirstOrder3DSolver.cpp

${OBJECTDIR}/Logger/Logger.o: Logger/Logger.cpp
	${MKDIR} -p ${OBJECTDIR}/Logger
	${RM} "$@.d"
	$(COMPILE.cc) -Wall -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Logger/Logger.o Logger/Logger.cpp

${OBJECTDIR}/MatrixClasses/BaseFloatMatrix.o: MatrixClasses/BaseFloatMatrix.cpp
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} "$@.d"
	$(COMPILE.cc) -Wall -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MatrixClasses/BaseFloatMatrix.o MatrixClasses/BaseFloatMatrix.cpp

${OBJECTDIR}/MatrixClasses/BaseIndexMatrix.o: MatrixClasses/BaseIndexMatrix.cpp
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} "$@.d"
	$(COMPILE.cc) -Wall -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MatrixClasses/BaseIndexMatrix.o MatrixClasses/BaseIndexMatrix.cpp

${OBJECTDIR}/MatrixClasses/ComplexMatrix.o: MatrixClasses/ComplexMatrix.cpp
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} "$@.d"
	$(COMPILE.cc) -Wall -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MatrixClasses/ComplexMatrix.o MatrixClasses/ComplexMatrix.cpp

${OBJECTDIR}/MatrixClasses/FftwComplexMatrix.o: MatrixClasses/FftwComplexMatrix.cpp
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} "$@.d"
	$(COMPILE.cc) -Wall -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MatrixClasses/FftwComplexMatrix.o MatrixClasses/FftwComplexMatrix.cpp

${OBJECTDIR}/MatrixClasses/IndexMatrix.o: MatrixClasses/IndexMatrix.cpp
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} "$@.d"
	$(COMPILE.cc) -Wall -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MatrixClasses/IndexMatrix.o MatrixClasses/IndexMatrix.cpp

${OBJECTDIR}/MatrixClasses/RealMatrix.o: MatrixClasses/RealMatrix.cpp
	${MKDIR} -p ${OBJECTDIR}/MatrixClasses
	${RM} "$@.d"
	$(COMPILE.cc) -Wall -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MatrixClasses/RealMatrix.o MatrixClasses/RealMatrix.cpp

${OBJECTDIR}/OutputStreams/BaseOutputStream.o: OutputStreams/BaseOutputStream.cpp
	${MKDIR} -p ${OBJECTDIR}/OutputStreams
	${RM} "$@.d"
	$(COMPILE.cc) -Wall -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/OutputStreams/BaseOutputStream.o OutputStreams/BaseOutputStream.cpp

${OBJECTDIR}/OutputStreams/CuboidOutputStream.o: OutputStreams/CuboidOutputStream.cpp
	${MKDIR} -p ${OBJECTDIR}/OutputStreams
	${RM} "$@.d"
	$(COMPILE.cc) -Wall -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/OutputStreams/CuboidOutputStream.o OutputStreams/CuboidOutputStream.cpp

${OBJECTDIR}/OutputStreams/IndexOutputStream.o: OutputStreams/IndexOutputStream.cpp
	${MKDIR} -p ${OBJECTDIR}/OutputStreams
	${RM} "$@.d"
	$(COMPILE.cc) -Wall -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/OutputStreams/IndexOutputStream.o OutputStreams/IndexOutputStream.cpp

${OBJECTDIR}/OutputStreams/WholeDomainOutputStream.o: OutputStreams/WholeDomainOutputStream.cpp
	${MKDIR} -p ${OBJECTDIR}/OutputStreams
	${RM} "$@.d"
	$(COMPILE.cc) -Wall -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/OutputStreams/WholeDomainOutputStream.o OutputStreams/WholeDomainOutputStream.cpp

${OBJECTDIR}/Parameters/CommandLineParameters.o: Parameters/CommandLineParameters.cpp
	${MKDIR} -p ${OBJECTDIR}/Parameters
	${RM} "$@.d"
	$(COMPILE.cc) -Wall -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Parameters/CommandLineParameters.o Parameters/CommandLineParameters.cpp

${OBJECTDIR}/Parameters/Parameters.o: Parameters/Parameters.cpp
	${MKDIR} -p ${OBJECTDIR}/Parameters
	${RM} "$@.d"
	$(COMPILE.cc) -Wall -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Parameters/Parameters.o Parameters/Parameters.cpp

${OBJECTDIR}/main.o: main.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -Wall -I./ -I${EBROOTHDF5}/include -I${EBROOTFFTW}/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

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
