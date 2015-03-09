# /**
# * @File        Makefile
# * @Author      Jiri Jaros
# * @Affiliation FIT, Brno University of Technology
# * @Email       jarosjir@fit.vutbr.cz
# * @Comments    Linux makefile for Anselm
# * 
# * @Tool        kspaceFirstOrder3D 2.15
# * @Created     24 September 2012, 10:57 AM
# * @LastModif   11 February 2014, 17:55 PM
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


#load following modues for GNU
# module load gcc/4.8.1
# module load hdf5/1.8.13
# module load fftw3/3.3.3-gcc


#load following modules for Intel
#module load intel/13.5.192(default)
#module load hdf5/1.8.13

#HDF5=$(HOME)/Software/hdf5-1.8.10-serial-icc


#Get GIT hash (only if you biuld form Gitlab repository)
KWAVE_GIT_HASH=$(shell git rev-parse HEAD)
#Otherwise, use the last official build hash

#################################################################################



############################## GNU g++ + FFTW ###################################
ifeq ($(COMPILER),GNU)
  #GNU version
  CXX	   = g++


  CXXFLAGS = -mtune=native -march=native -mavx -O3 -fopenmp -ffast-math -fassociative-math -Wall \
		     -I$(HDF5_INC_DIR) -I$(FFTW3_INC_DIR) -I. \
		     -D__KWAVE_GIT_HASH__=\"$(KWAVE_GIT_HASH)\"
  
  ifeq ($(LINKING),STATIC)

	LDFLAGS  = -fopenmp $(CPU_FLAGS) -L$(HDF5_LIB_DIR) -L$(FFTW3_LIB_DIR) -Wl,-rpath,$(HDF5_LIB_DIR)

	LIBS     = $(FFTW3_LIB_DIR)/libfftw3f.a  	\
		   $(FFTW3_LIB_DIR)/libfftw3f_omp.a 	\
		     -lhdf5 -lhdf5_hl
    else	
	LDFLAGS  = -fopenmp $(CPU_FLAGS) -L$(HDF5_LIB_DIR) -L$(FFT_LIB_DIR)
	LIBS     = -lfftw3f  -lfftw3f_omp -lhdf5 -lhdf5_hl
  endif
endif


########################### Intel Compiler + MKL ################################
ifeq ($(COMPILER),Intel)

  #intel 2013 compilers
  CXX	   = icpc


		
  CXXFLAGS = -fast -xhost -xavx -openmp  $(CPU_FLAGS) -Wall -ipo \
    	     -I$(HDF5_INC_DIR) -I$(MKL_INC_DIR) -I$(MKL_INC_DIR)/fftw -I./ 
#    	     -I$(HDF5)/include -I$(MKL_INC_DIR) -I$(MKL_INC_DIR)/fftw -I./ 


  ifeq ($(LINKING),STATIC)
	LDFLAGS  = -fast -openmp $(CPU_FLAGS) -L$(HDF5)/lib \
		   -mkl=parallel	
		   
#		   -L/apps/intel/inspector_xe_2013/lib64/pinruntime/glibc \
#		   -L/apps/libs/zlib/intel/1.2.8/lib		\
#		   -L/apps/libs/szip/intel/2.1/lib		\
#		   -mkl=parallel

	LIBS = 	$(HDF5)/lib/libhdf5_hl.a 			\
		$(HDF5)/lib/libhdf5.a 				\
		/apps/libs/zlib/intel/1.2.8/lib/libz.a 		\
		-lm


#		/apps/libs/zlib/intel/1.2.8/lib/libz.a 		\
#		/apps/libs/szip/intel/2.1/lib/libsz.a		\
#		-lm

  else	
	LDFLAGS  = -openmp $(CPU_FLAGS) -L$(HDF5_LIB_DIR) -L. \
                   -Wl,--start-group $(MKL_LIB_DIR)/libmkl_intel_ilp64.a $(MKL_LIB_DIR)/libmkl_core.a

#	-L/apps/libs/zlib/intel/1.2.8/lib/	\
#	-L/apps/libs/szip/intel/2.1/lib/	\

  	LIBS     = -lhdf5_hl -lhdf5
  endif
endif


################################# Compile #####################################


TARGET		= kspaceFirstOrder3D-OMP

all:		$(TARGET)	


$(TARGET):	main.o 					\
		HDF5/HDF5_File.o			\
		KSpaceSolver/KSpaceFirstOrder3DSolver.o	\
		MatrixClasses/BaseFloatMatrix.o		\
		MatrixClasses/BaseIndexMatrix.o		\
		MatrixClasses/ComplexMatrix.o		\
		MatrixClasses/FFTWComplexMatrix.o	\
		MatrixClasses/IndexMatrix.o		\
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
		MatrixClasses/BaseIndexMatrix.o		\
		MatrixClasses/ComplexMatrix.o		\
		MatrixClasses/FFTWComplexMatrix.o	\
		MatrixClasses/IndexMatrix.o		\
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


