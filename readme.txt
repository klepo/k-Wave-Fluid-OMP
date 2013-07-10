This file is part of the C++ extension of the k-Wave Toolbox (http://www.k-wave.org).\n
Copyright (C) 2013 Jiri Jaros and Bradley Treeby
 
This file is part of k-Wave. k-Wave is free software: you can redistribute it 
and/or modify it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of the License, 
or (at your option) any later version.
 
k-Wave is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU Lesser General Public License for more details. 

You should have received a copy of the GNU Lesser General Public License 
along with k-Wave. If not, see http://www.gnu.org/licenses/.
 

--------------------------------------------------------------------------------
Content of the repo
--------------------------------------------------------------------------------
Code - C++ code
	Doxygen - Full doxygen documentation
	HDF5 - C++ sources of the HDF5 module
	KSpaceSolver - the core of the C++ code
	MatrixClasses - C++ sources of different matrices classes
	Parameters - C++ sources of commandline and simulation parameters
	Utils - C++ sources of useful functions and constants
	build - NetBeans build directory
	dist - NetBeans distribution directory
	nbproject - NetBeans project

Data    - Simple input data
GNUMakefile - GNU Make file. If you want to use this makefile, copy it into the Code directory, and set the correct paths to HDF5, FFTW/MKL	
Manual - PDF version of the k-Wave manual
MatlabScripts 
	- GenerateInputData - use saveHDF5Script.m to generate input data.
	- ValidationTests
		- compareH5Files.m - compares two different HDF5 files and prints differences
		- V1p1_TEST_ALL.m - exhaustive validation tests (ensure the binary is under k-Wave Toolbox/binary) between MATLAB and C++
k-Wave Toolbox - MATLAB version of the toolbox. If you want to use the C++ code from within MATLAB or run the validation test, copy the binary into k-Wave Toolbox/binaries

