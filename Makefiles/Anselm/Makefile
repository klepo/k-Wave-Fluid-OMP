# /**
# * @File        Makefile
# * @Author      Jiri Jaros
# *              Faculty of Information Technology
# *              Brno University of Technology 
# * @Email       jarosjir@fit.vutbr.cz
# * @Comments    Linux makefile for Anselm supercomputer
# * 
# * @Tool        kspaceFirstOrder3D 2.17
# * @Created     04 September 2017, 12:50 
# * @LastModif   09 January   2019, 11:53
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
#   by GNU g++ 5.0 and newer OR Intel Compiler icpc 2015 and newer             #
#   The newer compiler, the more advanced instruction set can be used          #
#   We recomend compilation with g++ 5.3/5.4 or icpc 16/17                     #
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

#load following modules for Intel
#ml HDF5/1.8.18-intel-2017a-serial


#Get GIT hash (only if you biuld form Gitlab repository)
KWAVE_GIT_HASH=$(shell git rev-parse HEAD)
#Otherwise, use the last official build hash

################################################################################



############################# GNU g++ + FFTW ###################################
ifeq ($(COMPILER),GNU)
  #GNU version
  CXX	   = g++


  # Set compiler flags and header files directories
  CXXFLAGS = -m64 -march=native -mtune=native -std=c++11 -\
             -O3 -fopenmp -ffast-math -fassociative-math -Wall \
	     -I$(HDF5_DIR)/include -I$(FFT_DIR)/include -I. \
             -D__KWAVE_GIT_HASH__=\"$(KWAVE_GIT_HASH)\"
  
  ifeq ($(LINKING),STATIC)
	# Static link
	LDFLAGS  = $(CPU_FLAGS)  -fopenmp               \
                   -L$(HDF5_DIR)/lib -L$(FFT_DIR)/lib -static


	LIBS     = $(FFT_DIR)/lib/libfftw3f.a  		\
		   $(FFT_DIR)/lib/libfftw3f_omp.a 	\
   	  	   $(HDF5_DIR)/lib/libhdf5_hl.a 	\
   		   $(HDF5_DIR)/lib/libhdf5.a 		\
		   $(ZLIB_DIR)/lib/libz.a		\
		   $(SZIP_DIR)/lib/libsz.a		\
		   -lm	-lmvec				\
		   -ldl  # we need this for HDF5-1.8.11 and newer

    else	
        # Dynamic link
	LDFLAGS  = $(CPU_FLAGS) -O3 -fopenmp                 \
                   -L$(HDF5_DIR)/lib -L$(FFT_DIR)/lib	     \
		   -Wl,-rpath,$(HDF5_DIR)/lib:$(FFT_DIR)/lib

	LIBS     = -lfftw3f -lfftw3f_omp -lhdf5 -lhdf5_hl -lm -lz
  endif
endif


########################### Intel Compiler + MKL ###############################
ifeq ($(COMPILER),Intel)

  #intel 2015 compilers
  CXX	   = icpc


		
  CXXFLAGS = -O3 -xhost -qopenmp -ipo -Wall \
    	     -I$(EBROOTHDF5)/include -I./ 
#    	     -I$(HDF5)/include -I$(MKL_INC_DIR) -I$(MKL_INC_DIR)/fftw -I./ 


  ifeq ($(LINKING),STATIC)
	LDFLAGS  = -O3 -qopenmp -xhost -ipo -L$(EBROOTHDF5)/lib \
		   -mkl=parallel -static-intel

	LIBS = 	$(EBROOTHDF5)/lib/libhdf5_hl.a 	\
		$(EBROOTHDF5)/lib/libhdf5.a 	\
		$(EBROOTZLIB)/lib/libz.a 	\
		$(EBROOTSZIP)/lib/libsz.a 	\
		-lm -ldl

  else	
        # Dynamic link
	LDFLAGS  = $(CPU_FLAGS) -O3 -xhost -qopenmp \
                   -L$(EBROOTHDF5)/lib		\
            	   -Wl,-rpath,$(EBROOTHDF5)/lib \
    		   -mkl=parallel
    		   
        LIBS     = -lhdf5 -lhdf5_hl -lz -lsz

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
		$(LIBS)					 \
		-o $@


$(TARGET).o : $(TARGET).cpp
	$(CXX) $(CXXFLAGS) -c $(TARGET).cpp


clean:
	rm -f Hdf5/*.o KSpaceSolver/*.o MatrixClasses/*.o Parameters/*.o Containers/*.o  OutputStreams/*.o Logger/*.o  $(TARGET)


