# /**
# * @File        Makefile
# * @Author      Jiri Jaros
# *              Faculty of Information Technology
# *              Brno University of Technology 
# * @Email       jarosjir@fit.vutbr.cz
# * @Comments    Linux makefile for Ubuntu 16:04
# * 
# * @Tool        kspaceFirstOrder3D 2.15
# * @Created     18 September 2014, 12:32 
# * @LastModif   02 September 2017, 22:32
#
# * @License: 
# * This file is part of the C++ extension of the k-Wave Toolbox (http://www.k-wave.org)
# * Copyright (C) 2014 Jiri Jaros and Bradley Treeby
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
#	by GNU g++ 5.0 and newer OR Intel Compiler icpc 11 and newer 		#
#	The newer compiler, the more advanced instruction set can be used	#
#	We recomend compilation with g++ 5.3/4.9 or icpc 16/17                  #
#                                                                               #
#	The code uses GNU compiler and STATIC linking by default		#
#										#
#	Necessary libraries:							#
#		- FFTW 3.0 and newer OR Intel MKL 11 and newer			#
#		- HDF5 version 1.8 and newer					#
#										#
#	How to compile libraries						#
#		- FFTW : download from "http://www.fftw.org/" 			#
#			run configure script with following parameters:		#
#			--enable-float --enable-sse --enable-openmp		#
#                       or                                                      # 
#			--enable-float --enable-avx --enable-openmp		#
#		- MKL  : Only if not using FFTW. Dowload from 			#
#			 "http://software.intel.com/en-us/intel-mkl" 	        #
#		- HDF5 : download from http://www.hdfgroup.org/HDF5/ 		#
#			run configure script with these parameters:	        #		
#			--enable-hl --enable-static --enable-shared     	#
#									      	#
#	This makefile enables static and dynamic compilation using 		#
#		1) GNU g++ + FFTW  						#
#		2) Intel icpc + Intel MKL					#
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
#LINKING = STATIC
LINKING = DYNAMIC

#Set up paths: FFT_DIR for FFTW or MKL for MKL
FFT_DIR=$(EBROOTFFTW)
#MKL_DIR=/opt/intel/composer_xe_2013/mkl
#HDF5_DIR=/usr/local/hdf5-1.8.10-serial-icc
HDF5_DIR=$(EBROOTHDF5)

# Select CPU architecture (what instruction set to be used). 
# The libraries such as FFTW, HDF5 and MKL are to be compiled under the same architecture
# e.g. if you want to use AVX in K-Wave, compile FFTW with --enable-avx
CPU_ARCH = native
#CPU_ARCH = SSE3
#CPU_ARCH = SSE4
#CPU_ARCH = AVX
#CPU_ARCH = AVX2

#Get GIT hash (only if you biuld form Gitlab repository)
KWAVE_GIT_HASH=$(shell git rev-parse HEAD)
#Otherwise, use the last official build hash
#KWAVE_GIT_HASH=43e7bcfd05ffa7497ee6a895e34df9a424f56a9e
#################################################################################


############################## GNU g++ + FFTW ###################################
ifeq ($(COMPILER),GNU)
  CXX = g++

  #Set CPU architecture flags
  ifeq ($(CPU_ARCH),SSE3)      #Intel Core 2
    CPU_FLAGS=-m64 -msse3

  else ifeq ($(CPU_ARCH),SSE4) #Intel Nehalem, Westmere
    CPU_FLAGS=-m64 -msse4

  else ifeq ($(CPU_ARCH),AVX)  #Sandy Bridge, Ivy Bridge
    CPU_FLAGS=-m64 -mavx

  else ifeq ($(CPU_ARCH),AVX2) #Haswell, Broadwell
    CPU_FLAGS=-m64 -mavx2

  else #defautl is native - the max perfromance for this CPU
     CPU_FLAGS=-m64 -march=native -mtune=native
  endif


  # Set compiler flags and header files directories
  CXXFLAGS = $(CPU_FLAGS) -std=c++11 -O3 -fopenmp -ffast-math -fassociative-math -Wall \
	     -I$(HDF5_DIR)/include -I$(FFT_DIR)/include -I. \
             -D__KWAVE_GIT_HASH__=\"$(KWAVE_GIT_HASH)\"
  
  ifeq ($(LINKING),STATIC)
	# Static link
	LDFLAGS  = $(CPU_FLAGS)  -fopenmp               \
                   -L$(HDF5_DIR)/lib -L$(FFT_DIR)/lib  --static

	LIBS     = $(FFT_DIR)/lib/libfftw3f.a  		\
		   $(FFT_DIR)/lib/libfftw3f_omp.a 	\
   	  	   $(HDF5_DIR)/lib/libhdf5_hl.a 	\
   		   $(HDF5_DIR)/lib/libhdf5.a 		\
		  -lz  					\
		  -ldl  # we need this for HDF5-1.8.11 and newer

    else	
        # Dynamic link
	LDFLAGS  = $(CPU_FLAGS) -O3 -fopenmp                 \
                   -L$(HDF5_DIR)/lib -L$(FFT_DIR)/lib	     \
		   -Wl,-rpath,$(HDF5_DIR)/lib:$(FFT_DIR)/lib

	LIBS     = -lfftw3f  -lfftw3f_omp -lhdf5 -lhdf5_hl -lz
  endif
endif


########################### Intel Compiler + MKL ################################
ifeq ($(COMPILER),Intel)	
  CXX	   = icpc

  #Set CPU architecture flags
  ifeq ($(CPU_ARCH),SSE3)	#Intel, Core 2
    CPU_FLAGS=-m64 -march=core2

  else ifeq ($(CPU_ARCH),SSE4) #Intel Nehalem, Westmere
    CPU_FLAGS=-m64 -march=corei7

  else ifeq ($(CPU_ARCH),AVX)  #Sandy Bridge, Ivy Bridge
    CPU_FLAGS=-m64 -march=corei7-avx

  else ifeq ($(CPU_ARCH),AVX2) #Haswell, Broadwell
    CPU_FLAGS=-m64 -march=core-avx2

  else #defautl is native - the max perfromance for this CPU
     CPU_FLAGS=-m64 -march=native
  endif


  # Set compiler flags and header files directories		
  CXXFLAGS = -O3 -ipo -openmp  $(CPU_FLAGS) -Wall                                     \
  	     -I$(HDF5_DIR)/include -I$(MKL_DIR)/include -I$(MKL_DIR)/include/fftw -I. \
             -D__KWAVE_GIT_HASH__=\"$(KWAVE_GIT_HASH)\"


  ifeq ($(LINKING),STATIC)
	# Static link
	LDFLAGS  = $(CPU_FLAGS) -O3 -ipo -openmp  -static    \
                  -L$(MKL_DIR)/lib/intel64 -L$(HDF5_DIR)/lib

	LIBS = 	$(HDF5_DIR)/lib/libhdf5_hl.a 			              \
		$(HDF5_DIR)/lib/libhdf5.a 			              \
 	        -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a \
                      $(MKL_DIR)/lib/intel64/libmkl_core.a 		      \
                      $(MKL_DIR)/lib/intel64/libmkl_intel_thread.a            \
                -Wl,--end-group                                               \
                -lpthread -lm -lz					      
#		-ldl  # we need this for HDF5-1.8.11                                           

  else	
        # Dynamic link
	LDFLAGS  = $(CPU_FLAGS) -O3 -ipo -openmp             \
                   -L$(HDF5_DIR)/lib                         \
		   -Wl,-rpath,$(HDF5_DIR)/lib:$(FFT_DIR)/lib \
		   -mkl=parallel
  	LIBS     = -lhdf5 -lhdf5_hl -lz
  endif
endif


################################# Compile #####################################

TARGET		= kspaceFirstOrder3D-OMP

all:		$(TARGET)	


$(TARGET):	main.o 					 \
	        Containers/MatrixContainer.o             \
		Containers/MatrixRecord.o                \
		Containers/OutputStreamContainer.o       \
		Hdf5/Hdf5File.o				 \
		Hdf5/Hdf5FileHeader.o			 \
		KSpaceSolver/KSpaceFirstOrder3DSolver.o	 \
		Logger/Logger.o				 \
		MatrixClasses/BaseFloatMatrix.o		 \
		MatrixClasses/BaseIndexMatrix.o		 \
		MatrixClasses/ComplexMatrix.o		 \
		MatrixClasses/FftwComplexMatrix.o	 \
		MatrixClasses/IndexMatrix.o		 \
		MatrixClasses/RealMatrix.o		 \
		OutputStreams/BaseOutputStream.o         \
		OutputStreams/IndexOutputStream.o        \
		OutputStreams/CuboidOutputStream.o       \
		OutputStreams/WholeDomainOutputStream.o  \
		Parameters/CommandLineParameters.o	 \
		Parameters/Parameters.o			 \


	$(CXX) $(LDFLAGS) main.o 			 \
		Containers/MatrixContainer.o             \
		Containers/MatrixRecord.o                \
		Containers/OutputStreamContainer.o       \
		Hdf5/Hdf5File.o				 \
		Hdf5/Hdf5FileHeader.o			 \
		KSpaceSolver/KSpaceFirstOrder3DSolver.o	 \
		Logger/Logger.o				 \
		MatrixClasses/BaseFloatMatrix.o		 \
		MatrixClasses/BaseIndexMatrix.o		 \
		MatrixClasses/ComplexMatrix.o		 \
		MatrixClasses/FftwComplexMatrix.o	 \
		MatrixClasses/IndexMatrix.o		 \
		MatrixClasses/RealMatrix.o		 \
		OutputStreams/BaseOutputStream.o         \
		OutputStreams/IndexOutputStream.o        \
		OutputStreams/CuboidOutputStream.o       \
		OutputStreams/WholeDomainOutputStream.o  \
		Parameters/CommandLineParameters.o	 \
		Parameters/Parameters.o			 \
		$(LIBS)	\
		-o $@


$(TARGET).o : $(TARGET).cpp
	$(CXX) $(CXXFLAGS) -c $(TARGET).cpp


clean:
	rm -f *.o Hdf5/*.o KSpaceSolver/*.o MatrixClasses/*.o Parameters/*.o Containers/*.o  OutputStreams/*.o Logger/*.o $(TARGET)


