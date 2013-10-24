# /**
# * @File        Makefile
# * @Author      Jiri Jaros
# * @Affiliation FIT, Brno University of Technology
# * @Email       jarosjir@fit.vutbr.cz
# * @Comments    Linux makefile for Raijin
# * 
# * @Tool        kspaceFirstOrder3D 2.14
# * @Created     24 September 2012, 10:57 AM
# * @LastModif   24 October 2013, 22:50 PM
#
# * @License: 
# * This file is part of the C++ extension of the k-Wave Toolbox (http://www.k-wave.org)
# * Copyright (C) 2013 Jiri Jaros and Bradley Treeby
# * 
# * This file is part of the k-Wave. k-Wave is free software: you can redistribute it 
# * and/or modify it under the terms of the GNU Lesser General Public License as 
# * published by the Free Software Foundation, either version 3 of the License, 
# * or (at your option) any later version.
# * 
# * k-Wave is distributed in the hope that it will be useful, but 
# * WITHOUT ANY WARRANTY; without even the implied warranty of 
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
# * See the GNU Lesser General Public License for more details. 
# * 
# * You should have received a copy of the GNU Lesser General Public License 
# * along with k-Wave. If not, see <http://www.gnu.org/licenses/>.
# */

#################################################################################
#	The source codes can be compiled ONLY under Linux x64 			#
#	by GNU g++ 4.3 and newer OR Intel Compiler icpc 11			#
#	The code uses STATIC Linking by default					#
#										#
#	Necessary libraries:							#
#		- FFTW 3.0 and newer OR Intel MKL 11 and newer			#
#		- HDF5 version 1.8 and newer					#
#										#
#	How to compile libraries						#
#		- FFTW : download from "http://www.fftw.org/" 			#
#			run configure script with following parameters:		#
#			--enable-float --enable-sse --enable-openmp		#
#		- MKL  : Only if not using FFTW. Dowload from 			#
#			 "http://software.intel.com/en-us/intel-mkl" 	        #
#		- HDF5 : download from http://www.hdfgroup.org/HDF5/ 		#
#			run configure script with these parameters:	        #		
#			--enable-hl 						#
#									      	#
#	This makefile enables static and dynamic compilation using 		#
#		1) GNU gcc + FFTW  						#
#		2) Intel compiler + Intel MKL					#
#									      	#
#									      	#
#################################################################################


#################################################################################
#	  Set following flags based on your compiler and library paths 		#
#################################################################################

# Select compiler
COMPILER = GNU
#COMPILER = Intel

# static lining is deafult
LINKING = STATIC
#LINKING = DYNAMIC

#instruction set (plain,sse,avx)
FFTW_INSTR_SET = avx




#set up paths: FFT_DIR for FFTW or MKL for MKL
FFT_DIR=$(HOME)/Software/fftw3
#MKL_DIR=/opt/intel/composer_xe_2011_sp1/mkl
#HDF5_DIR=$(HOME)/Software/hdf5_serial



#################################################################################



############################## GNU g++ + FFTW ###################################
ifeq ($(COMPILER),GNU)
  #GNU version 4.4.7
  CXX	   = g++


#Generic AMD and Intel CPUS
ifeq ($(FFTW_INSTR_SET),plain)
  CPU_FLAGS=-msse2 -m64
  FFT_DIR:= $(FFT_DIR)/plain
endif

#Intel Nehalem/Westmere and similar
ifeq ($(FFTW_INSTR_SET),sse)
  CPU_FLAGS=-msse4.2 -m64
  FFT_DIR:= $(FFT_DIR)/sse2
endif

#Intel Sandy Bridge 
ifeq ($(FFTW_INSTR_SET),avx)
  CPU_FLAGS=-mavx -m64
  FFT_DIR:= $(FFT_DIR)/avx
endif

  CXXFLAGS = -O3 -fopenmp $(CPU_FLAGS) -ffast-math -fassociative-math -Wall \
		     -I$(HDF5_DIR)/include -I$(FFT_DIR)/include -I.
  
  ifeq ($(LINKING),STATIC)

	LDFLAGS  = -fopenmp $(CPU_FLAGS) -L$(HDF5_DIR)/lib -L$(FFT_DIR)/lib 

	LIBS     = $(FFT_DIR)/lib/libfftw3f.a  		\
		   $(FFT_DIR)/lib/libfftw3f_omp.a 	\
	  	    -lhdf5 -lhdf5_hl
#	  	   $(HDF5_DIR)/lib/libhdf5_hl.a 	\
#		   $(HDF5_DIR)/lib/libhdf5.a 		\
#	   	   $(ZLIB_DIR)/lib/libz.a		\

    else	
	LDFLAGS  = -fopenmp $(CPU_FLAGS) -L$(HDF5_DIR)/lib -L$(FFT_DIR)/lib
	LIBS     = -lfftw3f  -lfftw3f_omp -lhdf5 -lhdf5_hl
  endif
endif


########################### Intel Compiler + MKL ################################
ifeq ($(COMPILER),Intel)

  #intel 2013 compilers
  CXX	   = icpc

  #Generic AMD and Intel CPUS
ifeq ($(FFTW_INSTR_SET),plain)
  CPU_FLAGS=-msse2 -m64
endif

#Intel Nehalem/Westmere and similar
ifeq ($(FFTW_INSTR_SET),sse)
  CPU_FLAGS=-msse4.2 -m64
endif

#Intel Sandy Bridge 
ifeq ($(FFTW_INSTR_SET),avx)
  CPU_FLAGS=-xhost
  
  # -fno-alias -xhost
endif



		
  CXXFLAGS = -fast -openmp  $(CPU_FLAGS)  -ipo -Wall \
  	     -I$(HDF5_DIR)/include -I$(INTEL_MKL_ROOT)/include -I$(INTEL_MKL_ROOT)/include/fftw -I./

  ifeq ($(LINKING),STATIC)
	LDFLAGS  = -fast -openmp $(CPU_FLAGS) -L$(HDF5_DIR)/lib -mkl=parallel -ipo

	LIBS = 	$(HDF5_DIR)/lib/libhdf5_hl.a 			\
		$(HDF5_DIR)/lib/libhdf5.a 			\
		-lz
  else	
	LDFLAGS  = -openmp $(CPU_FLAGS) -L$(HDF5_DIR)/lib  -mkl=parallel -ipo
  	LIBS     = -lhdf5 -lhdf5_hl -lz
  endif
endif


################################# Compile #####################################


TARGET		= kspaceFirstOrder3D-OMP

all:		$(TARGET)	


$(TARGET):	main.o 					\
		HDF5/HDF5_File.o			\
		KSpaceSolver/KSpaceFirstOrder3DSolver.o	\
		MatrixClasses/BaseFloatMatrix.o		\
		MatrixClasses/BaseLongMatrix.o		\
		MatrixClasses/ComplexMatrix.o		\
		MatrixClasses/FFTWComplexMatrix.o	\
		MatrixClasses/LongMatrix.o		\
		MatrixClasses/MatrixContainer.o		\
		MatrixClasses/OutputHDF5Stream.o	\
		MatrixClasses/RealMatrix.o		\
		MatrixClasses/UXYZ_SGXYZMatrix.o	\
		Parameters/CommandLineParameters.o	\
		Parameters/Parameters.o			\


	$(CXX) $(LDFLAGS) main.o 			\
		HDF5/HDF5_File.o			\
		KSpaceSolver/KSpaceFirstOrder3DSolver.o	\
		MatrixClasses/BaseFloatMatrix.o		\
		MatrixClasses/BaseLongMatrix.o		\
		MatrixClasses/ComplexMatrix.o		\
		MatrixClasses/FFTWComplexMatrix.o	\
		MatrixClasses/LongMatrix.o		\
		MatrixClasses/MatrixContainer.o		\
		MatrixClasses/OutputHDF5Stream.o	\
		MatrixClasses/RealMatrix.o		\
		MatrixClasses/UXYZ_SGXYZMatrix.o	\
		Parameters/CommandLineParameters.o	\
		Parameters/Parameters.o			\
		$(LIBS)					\
		-o $@


$(TARGET).o : $(TARGET).cpp
	$(CXX) $(CXXFLAGS) -c $(TARGET).cpp


clean:
	rm -f *.o HDF5/*.o KSpaceSolver/*.o MatrixClasses/*.o Parameters/*.o $(TARGET)


