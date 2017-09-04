## Overview

k-Wave is an open source MATLAB toolbox designed for the time-domain simulation
of propagating acoustic waves in 1D, 2D, or 3D. The toolbox has a wide range of 
functionality, but at its heart is an advanced numerical model that can account 
for both linear or nonlinear wave propagation, an arbitrary distribution of 
weakly heterogeneous material parameters, and power law acoustic  absorption, 
see the http://www.k-wave.org.

This project is a part of the k-Wave toolbox accelerating 3D simulations using 
an optimized C++ implementation to run moderate to big grid sizes (64x64x64 to 
512x512x512).

## Repository structure

    .
    +--Containers        - Matrix and output stream containers
    +--Data              - Small test data
    +--GetoptWin64       - Windows version of the getopt routine
    +--Hdf5              - HDF5 classes (file access)
    +--KSpaceSolver      - Solver classes with all the kernels
    +--Logger            - Logger class to report progress and errors
    +--MatrixClasses     - Matrix classes to hold data
    +--Makefiles         - GNU makefiles for different systems
    +--OutputHDF5Streams - Output streams to sample data
    +--Parameters        - Parameters of the simulation
    +--Utils             - Utility routines
    +--nbproject         - NetBeans IDE 8.2 project
    Changelog.md         - Changelog
    License.md           - License file
    Makefile             - NetBeans 8.2 makefile
    Readme.md            - Readme
    Doxyfile             - Doxygen generator file
    header_bg.png        - Doxygen logo
    main.cpp             - Main file of the project


## Compilation

The source codes of `kpsaceFirstOrder3D-OMP` are written using the C++-11 
standard, the OpenMP 4.0 library, FFTW 3.3.6 or MKL 11 and HDF5 1.8.x library

There are variety of different C++ compilers that can be used to compile the 
source codes. We recommend using either the GNU C++ compiler (gcc/g++) version 
5.0 and higher, or the Intel C++ compiler version 2015 and higher. Please note
that Visual Studio compilers do not support OpenMP 4.0 standard and cannot be 
used thus. Be also aware of Intel compilers 2017 and their MKL bug producing 
incorrect results when AVX2 is enabled. The codes can be compiled on 64-bit 
Linux and Windows.  32-bit systems are not supported!

This section describes the compilation procedure using GNU and Intel compilers 
on Linux.  (The Windows users are encouraged to download the Visual Studio 2015 
project and compile it using Intel Compiler from within Visual Studio.)

Before compiling the code, it is necessary to install a C++ compiler and several
libraries.  The GNU compiler is usually part of Linux distributions and 
distributed as open source.  It can be downloaded from http://gcc.gnu.org/) 
if necessary.

The Intel compiler can be downloaded from 
http://software.intel.com/en-us/intel-composer-xe/. This package also includes 
the Intel MKL (Math Kernel Library)library that contains FFT. The Intel compiler
is only free for non-commercial use.

The code also relies on several libraries that are to be installed before 
compiling:

 1. HDF5 library - Mandatory I/O library, version 1.8.x, 
         https://support.hdfgroup.org/HDF5/release/obtain518.html.
 1. FFTW library - Optional library for FFT, version 3.3.x,
         http://www.fftw.org/.
 1. MKL library  - Optional library for FFT, version 2015 or higher
         http://software.intel.com/en-us/intel-composer-xe/.

Although it is possible to use any combination of the FFT library and the 
compiler, the best performance is observed when using GNU compiler and FFTW, 
or Intel Compiler and Intel MKL.

### The HDF5 library installation procedure

 1. Download the 64-bit HDF5 library 
 https://support.hdfgroup.org/HDF5/release/obtain518.html. Please use version 
 1.8.x, the version 1.10.x is not compatible with MATLAB yet.
  
 2. Configure the HDF5 distribution. Enable the high-level library and specify 
 an installation folder by typing:
    ```bash
    ./configure --enable-hl --prefix=folder_to_install
    ```
 3. Make the HDF library by typing:
    ```bash
    make -j
    ```
 4. Install the HDF5 library by typing:
    ```bash
    make install


### The FFTW library installation procedure

 1. Download the FFTW library package for your platform, 
    http://www.fftw.org/download.html).

 2. Configure the FFTW distribution. Enable OpenMP support, and desired SIMD 
    instruction set, single precision floating point arithmetic, and specify an
    installation folder:
    ```bash
    ./configure --enable-single --enable-sse --enable-openmp --enable-shared \
                --prefix=folder_to_install
    ```

    if you intend to use the FFTW library (and the C++ code) only on a selected 
    machine and want to get the best possible performance, you may also add 
    processor specific optimisations and AVX or AVX2 instructions set. Note,
    the compiled binary code is not likely to be portable on different CPUs. 
    SSE2 version will work on any processor, AVX on Intel Sandy Bridge and 
    newer, AVX2 on Intel Haswell and newer. 
    ```bash
    ./configure --enable-single --enable-avx --enable-openmp  --enable-shared \
                --with-gcc-arch=native --prefix=folder_to_install
    ```
   
    More information about the installation and customization can be found at 
    http://www.fftw.org/fftw3_doc/Installation-and-Customization.htm. 
    For recent CPUs based on Sandy Bridge, Ivy Bridge, Haswell and Broadwell with
    strongly recommend to use the AVX and AVX2 support.

 3. Make the FFTW library by typing:
    ```bash
    make -j
    ```
 4. Install the FFTW library by typing:
    ```bash
    make install
    ````


### The Intel Compiler and MKL installation procedure

 1. Download  the Intel Composer XE package for your platform 
    http://software.intel.com/en-us/intel-compilers.
 
 2. Run the installation script and follow the procedure by typing:
    ```bash
    ./install.sh
    ```

### Compiling the C++ code on Linux

After the libraries and the compiler have been installed, you are ready to 
compile the `kspaceFirstOrder3D-OMP` code.

 1. Select the most appropriate makefile. 
    We recommend `Makefiles/Release/Makefile`
     
 2. The Makefile supports code compilation under GNU compiler and FFTW, or Intel
    compiler with MKL. Uncomment the desired compiler by removing character `#`.
    ```bash
    #COMPILER = GNU
    #COMPILER = Intel
    ```
 
 3. Select how to link the libraries. Static linking is preferred as it may be 
    a bit faster, however, on some systems (HPC clusters) it may be better to 
    use dynamic linking and use the system specific libraries at runtime.
    ```bash
    #LINKING = STATIC
    #LINKING = DYNAMIC
    ```	

 4. Select the instruction set and the CPU architecture.
    For users who will only use the binary on the same machine as compiled, 
    the best choice is `CPU_ARCH=native`. If you are about to run the same 
    binary on different machines or you want to cross-compile the code, you are 
    free to use any of the possible choices, where SSE 3 is the most general but
    slowest and AVX2 is the most recent instruction set while believed to be the
    fastest one.
    ```bash
    #CPU_ARCH = native
    #CPU_ARCH = SSE3
    #CPU_ARCH = SSE4
    #CPU_ARCH = AVX
    #CPU_ARCH = AVX2
    ```	
 5. Set installation paths of the libraries (an example is shown bellow). Zlib 
    and SZIP may be required if the compression is switched on.
    ```bash
    MKL_DIR=
    FFT_DIR=
    HDF5_DIR=
    ZLIB_DIR=
    SZIP_DIR=
    ```	
 6. Compile the source code by typing:
    ```bash
    make -j
    ```
    If you want to clean the distribution, type:
    ```bash
    make clean
    ```
 
Alternatively, the code can be compiled using Netbeans IDE and attached project.

## Usage

The C++ codes offers a lot of parameters and output flags to be used. For more 
information, please type:

```bash
./kspaceFirstOrder3D-OMP --help
```
