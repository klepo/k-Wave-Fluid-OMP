# /**
# * @File        Makefile
# * @Author      Jiri Jaros
# *              Faculty of Information Technology
# *              Brno University of Technology 
# * @Email       jarosjir@fit.vutbr.cz
# * @Comments    Linux makefile for Salomon supercomputer
# * 
# * @Tool        kspaceFirstOrder 2.17
# * @Created     04 September 2017, 12:50 
# * @LastModif   20 February  2009, 14:33
#
# * @copyright Copyright (C) 2019 Jiri Jaros and Bradley Treeby.
# *
# * This file is part of the C++ extension of the [k-Wave Toolbox]
# *  (http://www.k-wave.org).
# *
# * This file is part of the k-Wave. k-Wave is free software: you can 
# * redistribute it and/or modify it under the terms
# * of the GNU Lesser General Public License as published by the Free Software 
# * Foundation, either version 3 of the License, or (at your option) any later 
# * version.
# *
# * k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY 
# * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
# * more details.
# *
# * You should have received a copy of the GNU Lesser General Public License 
# * along with k-Wave.
# * If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).
# */

################################################################################
#                                                                              #
#   The source codes can be compiled ONLY under Linux x64                      #
#   by GNU g++ 5.0 and newer OR Intel Compiler icpc 2016 and newer             #
#   The newer compiler, the more advanced instruction set can be used          #
#   We recomend compilation with g++ 7.3 or icpc 2018                          #
#                                                                              #
#   The code uses GNU compiler and STATIC linking by default                   #
#                                                                              #
#       Necessary libraries:                                                   #
#         - FFTW 3.3.x and newer OR Intel MKL 11 and newer                     #
#         - HDF5 version 1.8.x and newer                                       #
#                                                                              #
#    How to compile libraries                                                  #
#         - FFTW : download from "http://www.fftw.org/"                        #
#                  run configure script with following parameters:             #
#                  --enable-float --enable-sse --enable-openmp                 #
#                       or                                                     # 
#                  --enable-float --enable-avx --enable-openmp                 #
#         - MKL  : Only if not using FFTW. Dowload from                        #
#                  "http://software.intel.com/en-us/intel-mkl"                 #
#         - HDF5 : download from                                               #
#                   https://support.hdfgroup.org/HDF5/release/obtain518.html   #
#                  run configure script with these parameters:                 #
#                  --enable-hl --enable-static --enable-shared                 #
#                                                                              #
#   This makefile enables static and dynamic compilation using                 #
#        1) GNU g++ + FFTW                                                     #
#        2) Intel icpc + Intel MKL                                             #
#                                                                              #
################################################################################


################################################################################
#         Set following flags based on your compiler and library paths         #
################################################################################

# Select compiler
#COMPILER = GNU
COMPILER = Intel

# static lining is deafult
LINKING = STATIC
#LINKING = DYNAMIC


#load following modues for GNU
#there are no modules for GCC version, You have to first compile 
# Zlib, SZIP, HDF5, FFTW 

#load following modules for Intel. Intel 2017 not supported!!!
#ml HDF5/1.8.18-intel-2017a-serial
#ml intel/2017.00


#Get GIT hash (only if you biuld form Gitlab repository)
KWAVE_GIT_HASH=$(shell git rev-parse HEAD)
#Otherwise, use the last official build hash

FFTW3_DIR=
HDF5_DIR=$(EBROOTHDF5)
ZLIB_DIR=$(EBROOTZLIB)
SZIP_DIR=$(EBROOTSZIP)


#################################################################################



############################## GNU g++ + FFTW ###################################
ifeq ($(COMPILER),GNU)
  #GNU version
  CXX	   = g++


  CXXFLAGS = -mtune=native -march=native -O3 -fopenmp -ffast-math \
             -fassociative-math -Wall \
             -I$(HDF5_DIR)/include -I$(FFTW3_DIR)/include -I. \
	     -D__KWAVE_GIT_HASH__=\"$(KWAVE_GIT_HASH)\"
  
  ifeq ($(LINKING),STATIC)

	LDFLAGS  = -fopenmp -O3 -L$(HDF5_DIR)/lib \
	           -L$(FFTW3_DIR)/lib -Wl,-rpath,$(HDF5_LIB_DIR)

	LIBS     = $(FFTW3_DIR)/lib/libfftw3f.a  	\
		   $(FFTW3_DIR)/lib/libfftw3f_omp.a 	\
   		   $(ZLIB_DIR)/lib/libz.a 		\
   		   $(SZIP_DIR)/lib/libsz.a		\
 		   -lm -ldl
    else	
	LDFLAGS  = -fopenmp -O3  -L$(HDF5_DIR)/lib -L$(ZLIB_DIR)/lib \
	           -L$(SZIP_DIR)/lib L$(FFT_DIR)/lib  \
	           -Wl,-rpath,$(HDF5_DIR)/lib:$(ZLIB_DIR)/lib:$(SZIP_DIR)/lib
	LIBS     = -lfftw3f -lfftw3f_omp -lhdf5 -lhdf5_hl -lsz -lz
  endif
endif


########################### Intel Compiler + MKL ################################
ifeq ($(COMPILER),Intel)

  #intel 2016 compilers
  CXX	   = icpc


		
  CXXFLAGS = -O3 -xhost -ipo -qopenmp -Wall -std=c++11 \
    	     -I$(HDF5_DIR)/include -I$(MKLROOT)/inlcude \
    	     -I$(MKLROOT)/include/fftw -I./ 


  ifeq ($(LINKING),STATIC)
	LDFLAGS  = -O3 -qopenmp -xhost -ipo -L$(HDF5_DIR)/lib -L$(MKLROOT)/lib/intel64 \
		   -mkl=parallel -static-intel

	LIBS = 	$(HDF5_DIR)/lib/libhdf5_hl.a 	\
		$(HDF5_DIR)/lib/libhdf5.a 	\
		$(ZLIB_DIR)/lib/libz.a 		\
		$(SZIP_DIR)/lib/libsz.a		\
		-lm -ldl

  else	
        # Dynamic link
	LDFLAGS  = $(CPU_FLAGS) -O3 -qopenmp             \
                   -L$(HDF5_DIR)/lib -L$(ZLIB_DIR)/lib -L$(SZIP_DIR)/lib \
            	   -Wl,-rpath,$(HDF5_DIR)/lib:$(MKLROOT)/lib \
    		   -mkl=parallel
        LIBS     = -lhdf5 -lhdf5_hl -lz -lsz

  endif
endif


################################# Compile #####################################


TARGET		= kspaceFirstOrder-OMP

all:		$(TARGET)	


$(TARGET):	main.o 					 \
		Compression/CompressHelper.o             \
		Containers/MatrixContainer.o             \
		Containers/MatrixRecord.o                \
		Containers/OutputStreamContainer.o       \
		Hdf5/Hdf5File.o				 \
		Hdf5/Hdf5FileHeader.o			 \
		KSpaceSolver/KSpaceFirstOrderSolver.o	 \
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
		Compression/CompressHelper.o             \
		Containers/MatrixContainer.o             \
		Containers/MatrixRecord.o                \
		Containers/OutputStreamContainer.o       \
		Hdf5/Hdf5File.o				 \
		Hdf5/Hdf5FileHeader.o			 \
		KSpaceSolver/KSpaceFirstOrderSolver.o	 \
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
		$(LIBS)					 \
		-o $@


$(TARGET).o : $(TARGET).cpp
	$(CXX) $(CXXFLAGS) -c $(TARGET).cpp


clean:
	rm -f *.o Compression/*.o Hdf5/*.o KSpaceSolver/*.o MatrixClasses/*.o Parameters/*.o \
	       Containers/*.o  OutputStreams/*.o Logger/*.o $(TARGET)


