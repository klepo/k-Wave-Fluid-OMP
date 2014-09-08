/**
 * @file        main.cpp
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The main file
 *
 * @version     kspaceFirstOrder3D 2.15
 *
 * @date        11 July      2012, 10:57 (created) \n
 *              02 September 2014, 16:49 (revised)
 *
 *
 *
 * @section License
 * This file is part of the C++ extension of the k-Wave Toolbox (http://www.k-wave.org).\n
 * Copyright (C) 2014 Jiri Jaros and Bradley Treeby
 *
 * This file is part of k-Wave. k-Wave is free software: you can redistribute it
 * and/or modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * k-Wave is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with k-Wave. If not, see http://www.gnu.org/licenses/.
 *
 *
 *
 * @mainpage kspaceFirstOrder3D-OMP
 *
 * @section Overview 1 Overview
 *
 * k-Wave is an open source MATLAB toolbox designed for the time-domain simulation of propagating acoustic waves in
 * 1D, 2D, or 3D. The toolbox has a wide range of functionality, but at its heart is an advanced numerical model that
 * can account for both linear or nonlinear wave propagation, an arbitrary distribution of weakly heterogeneous material
 * parameters, and power law acoustic absorption (http://www.k-wave.org).
 *
 * This project is a part of the k-Wave toolbox accelerating 3D simulations using an optimized C++ implementation to run moderate to big grid sizes.
 * Compiled binaries of the C++ code for x86 architectures are available from (http://www.k-wave.org/download). Both 64-bit Linux (Ubuntu / Debian) and 64-bit Windows
 * versions are provided. This Doxygen documentation was created based on the Linux version and provides details on the implementation of the C++ code.
 *
 *
 *
 * @section Compilation 2 Compilation
 *
 *
 * The source codes of <tt>kpsaceFirstOrder3D-OMP</tt> are written using the C++03 standard and the OpenMP 2.0 library.
 * There are variety of different C++ compilers that can be used to compile the source codes. We recommend using either the GNU C++ compiler (gcc/g++)
 * version 4.3 and higher, or the Intel C++ compiler version 11.0 and higher. The codes can be compiled on 64-bit Linux and Windows.
 * 32-bit systems are not supported. This section describes the compilation procedure using GNU and Intel compilers on Linux.
 * (The Windows users are encouraged to download the Visual Studio 2010 project and compile it using Intel Compiler from within Visual Studio.)
 *
 * Before compiling the code, it is necessary to install a C++ compiler and several libraries.  The GNU compiler is usually part of
 * Linux distributions and distributed as open source.  It can be downloaded from http://gcc.gnu.org/ if necessary. The Intel compiler can be downloaded from
 * http://software.intel.com/en-us/intel-composer-xe/. This package also includes the Intel MKL (Math Kernel Library) library that contains FFT. The Intel compiler
 *  is only free for non-commercial use.
 *
 * The code also relies on several libraries that is to be installed before compiling:
 *
 * \li HDF5 library - Mandatory I/O library, version 1.8 or higher (http://www.hdfgroup.org/HDF5/).
 * \li FFTW library - Optional library for FFT, version 3.0 or higher (http://www.fftw.org/).
 * \li MKLlibrary - Optional library for FFT, version 11.0 or higher (http://software.intel.com/en-us/intel-composer-xe/).

 *
 * Although it is possible to use any combination of the FFT library and the compiler,
 * the best performance is observed when using GNU compiler and FFTW, or Intel Compiler and Intel MKL.

 * <b> 2.1 The HDF5 library installation procedure </b>

 * 1. Download the 64-bit HDF5 library  package for your platform (http://www.hdfgroup.org/HDF5/release/obtain5.html).
 *
 * 2. Configure the HDF5 distribution. Enable the high-level library and specify an installation folder by typing:
\verbatim
        ./configure --enable-hl --prefix=folder_to_install
\endverbatim
 * 3. Make the HDF library by typing:
\verbatim
        make
\endverbatim
 * 4. Install the HDF library by typing:
\verbatim
	make install
\endverbatim
 *
 *
 *
 * <b> 2.2 The FFTW library installation procedure </b>
 *
 * 1. Download the FFTW library  package for your platform (http://www.fftw.org/download.html).
 *
 * 2. Configure the FFTW distribution. Enable OpenMP support, SSE instruction set, single precision data type, and specify an installation folder:
\verbatim
        ./configure --enable-single --enable-sse --enable-openmp --enable-shared --prefix=folder_to_install
\endverbatim
     if you intend to use the FFTW library (and the C++ code) only on a selected machine and want to get the best possible performance, you may also add processor specific optimisations and AVX instructions set. Note, the compiled binary code is not likely to be portable on different CPUs (e.g. even from Intel Sandy Bridge to Intel Nehalem).
\verbatim
        ./configure --enable-single --enable-avx --enable-openmp  --enable-shared --with-gcc-arch=<arch> --prefix=folder_to_install
\endverbatim
 * More information about the installation and customization can be found at http://www.fftw.org/fftw3_doc/Installation-and-Customization.htm.
 *
 * 3. Make the FFTW library by typing:
\verbatim
	make
\endverbatim
 * 4. Install the FFTW library by typing:
\verbatim
	make install
\endverbatim

 * <b> 2.3 The Intel Compiler and MKL installation procedure </b>
 *
 *
 * 1. Download  the Intel Composer XE package for your platform (http://software.intel.com/en-us/intel-compilers).
 *
 * 2. Run the installation script and follow the procedure by typing:
\verbatim
	./install.sh
\endverbatim
 *
 *
 *
 * <b> 2.4 Compiling the C++ code </b>
 *
 * When the libraries and the compiler have been installed, you are ready to compile the <tt>kspaceFirstOrder3D-OMP </tt> code.
 *

 * 1. Download the <tt>kspaceFirstOrder3D-OMP </tt> souce codes.
 *
 * 2. Open the \c Makefile file.  The Makefile supports code compilation under GNU compiler and FFTW, or Intel compiler with MKL. Uncomment the desired compiler by removing the character of `<tt>#</tt>'.
\verbatim
	#COMPILER = GNU
	#COMPILER = Intel
\endverbatim
 *3. Select how to link the libraries. Static linking is preferred, however, on some systems (HPC clusters) it may be better to use dynamic linking and use the system specific library at runtime.
\verbatim
	#LINKING = STATIC
	#LINKING = DYNAMIC
\endverbatim
 *
 * 4. Set installation paths of the libraries (an example is shown bellow).
\verbatim
	FFT_DIR=/usr/local
	MKL_DIR=/opt/intel/composer_xe_2011_sp1/mkl
	HDF5_DIR=/usr/local/hdf5-1.8.9
\endverbatim
 * 5. The code will be built to work on all CPUs implementing the Streaming SSE2 instruction set (Intel Pentium IV, AMD Athlon XP 64 and newer). Higher performance may be obtained by providing more information about the target CPU (SSE 4.2 instructions, AVX instruction, architecture optimization). For more detail, check on the compiler documentation.
 *
 * 6. Compile the source code by typing:
\verbatim
        make
\endverbatim
   If you want to clean the distribution, type:
\verbatim
        make clean
\endverbatim
 *
 *
 *
 *
 * @section Parameters 3 Command Line Parameters
 The  C++ code requires two mandatory parameters and accepts a few optional parameters and flags.
 * The mandatory parameters \c -i and \c -o specify the input and output file. The file names respect the path conventions for particular operating system.
 * If any of the files is not specified, cannot be found or created, an error message is shown.
 *
 * The \c -t parameter sets the number of threads used, which defaults the system maximum. On CPUs with Intel Hyper-Threading (HT),
 * performance will generally be better if HT is disabled in the BIOS settings. If HT is switched on, the default will be to create
 * twice as many threads as there are physical processor cores, which might slow down the code execution. Therefore, if the HT is on,
 *  try specifying the number of threads manually for best performance (e.g., 4 for Intel Quad-Core). We recommend experimenting
 * with this parameter to find the best configuration. Note, if there are other tasks being executed on the system,
 * it might be useful to further limit the number of threads to prevent system overload.
 *
 * The \c -r parameter specifies how often information about the simulation progress is printed to the command line.
 * By default, the C++ code prints out the progress of the simulation, the elapsed time, and the estimated time of
 * completion in intervals corresponding to 5% of the total number of times steps.
 *
 * The \c -c parameter specifies the compression level used by the ZIP library to reduce the
 * size of the output file. The actual compression rate is highly dependent on the shape of
 * the sensor mask and the range of stored quantities. In general, the output data is very
 * hard to compress, and using higher compression levels can greatly increase the time to
 * save data while not having a large impact on the final file size. The default compression
 * level of 3 represents a balance between compression ratio and performance that is suitable in most cases.
 *
 * The \c --benchmark parameter enables the total length of simulation (i.e., the number of time steps)
 * to be overwritten by setting a new number of time steps to simulate. This is particularly useful
 * for performance evaluation and benchmarking. As the code performance is relatively stable, 50-100 time steps is
 * usually enough to predict the simulation duration. This parameter can also be used to quickly find the ideal number of CPU threads to use.
 *
 * The \c -h and \c --help parameters print all the parameters of the C++ code, while the \c --version parameter reports the code version and internal build number.
 *
 *
 * The remaining flags specify the output quantities to be recorded during the simulation and stored on disk analogous to
 * the sensor.record input. If the \c -p or \c --p raw flags are set (these are equivalent), a time series of the acoustic pressure at the
 * grid points specified by the sensor mask is recorded. If the \c --p rms and/or \c --p max flags are set, the root mean square and/or maximum values
 * of the pressure at the grid points specified by the sensor mask are recorded. Finally, if the \c --p final flag is set, the values for the entire acoustic pressure
 * field in the final time step of the simulation is stored (this will always include the PML, regardless of the setting for <tt> `PMLInside' </tt>). Flags to record the acoustic
 * particular velocity are defined in an analogous fashion.
 *
 *
 * Any combination of \c p, \c u and \c I  fags is admissible. If no output flag is set, a time-series for
 * the acoustic pressure is recorded. If it is not necessary to collect the output quantities over
 * the entire simulation, the starting time step when the collection begins can be specified
 * using the -s parameter. Note, the index for the first time step is 1 (this follows the MATLAB indexing convention).
 *
 *
\verbatim
---------------------------------- Usage ---------------------------------
Mandatory parameters:
  -i <input_file_name>            : HDF5 input file
  -o <output_file_name>           : HDF5 output file

Optional parameters:
  -t <num_threads>                : Number of CPU threads (default = MAX)
  -r <interval_in_%>              : Progress print interval (default = 5%)
  -c <comp_level>                 : Output file compression level <0,9>
                                      (default = 0)
  --benchmark <steps>             : Run a specified number of time steps

  --checkpoint_file <file_name>   : HDF5 checkpoint file
  --checkpoint_interval <seconds> : Stop after a given number of seconds and
                                      store the actual state

  -h                              : Print help
  --help                          : Print help
  --version                       : Print version

Output flags:
  -p                              : Store acoustic pressure
                                      (default if nothing else is on)
                                      (the same as --p_raw)
  --p_raw                         : Store raw time series of p (default)
  --p_rms                         : Store rms of p
  --p_max                         : Store max of p
  --p_min                         : Store min of p
  --p_max_all                     : Store max of p (whole domain)
  --p_min_all                     : Store min of p (whole domain)
  --p_final                       : Store final pressure field

  -u                              : Store ux, uy, uz
                                      (the same as --u_raw)
  --u_raw                         : Store raw time series of ux, uy, uz
  --u_non_staggered_raw           : Store non-staggered raw time series of
                                      ux, uy, uz
  --u_rms                         : Store rms of ux, uy, uz
  --u_max                         : Store max of ux, uy, uz
  --u_min                         : Store min of ux, uy, uz
  --u_max_all                     : Store max of ux ,uy, yz (whole domain)
  --u_min_all                     : Store max of ux ,uy, yz (whole domain)
  --u_final                       : Store final acoustic velocity

  --copy_sensor_mask              : Copy sensor mask to the output file

  -s <timestep>                   : Time step when data collection begins
                                      (default = 1)
--------------------------------------------------------------------------
\endverbatim
 *
 *
 *
 *
 *
 *
 * @section HDF5Files 4 HDF5 File Structure
 *
   * The C++ code has been designed as a standalone application which is not dependent on MATLAB libraries or a MEX interface.
 * This is of particular importance when using servers and supercomputers without MATLAB support.
 * For this reason, simulation data must be transferred between the C++ code and MATLAB using external input and output files.
 * These files are stored using the Hierarchical Data Format HDF5 (http://www.hdfgroup.org/HDF5/). This is a data model,
 * library, and file format for storing and managing data. It supports a variety of datatypes, and is designed for
 *  flexible and efficient I/O and for high volume and complex data. The HDF5 technology suite includes tools and applications for
 * managing, manipulating, viewing, and analysing data in the HDF5 format.
 *
 *
 *
 * Each HDF5 file is a container for storing a variety of scientific data and is composed
 * of two primary types of objects: groups and datasets. A HDF5 group is a structure
 * containing zero or more HDF5 objects, together with supporting metadata. A HDF5
 * group can be seen as a disk folder. A HDF5 dataset is a multidimensional array of
 * data elements, together with supporting metadata. A HDF5 dataset can be seen as a
 * disk file. Any HDF5 group or dataset may also have an associated attribute list. A
 * HDF5 attribute is a user-defined HDF5 structure that provides extra information about a
 * HDF5 object. More information can be obtained from the HDF5 documentation (http://www.hdfgroup.org/HDF5/doc/index.html).

 *
 * The HDF5 input file for the C++ simulation code contains a file header with brief description of the simulation
 * stored in string attributes, and the root group `/' which stores all the simulation properties in the form of 3D datasets
 * (a complete list of input datasets is given in Appendix B). The HDF5 output file contains a file header with the simulation
 * description as well as performance statistics, such as the simulation time and memory consumption, stored in string attributes.
 * The results of the simulation are stored in the root group `/' in the form of 3D datasets.
 *
 *
 * \verbatim
==============================================================================================================
                                        Input File Header
=============================================================================================================
created_by                              Short description of the tool that created this file
creation_date                           Date when the file was created
file_description                        Short description of the content of the file (e.g. simulation name)
file_type                               Type of the file (input)
major_version                           Major version of the file definition (1)
minor_version                           Minor version of the file definition (0)
==============================================================================================================
 \endverbatim
 *
 * \verbatim
==============================================================================================================
                                        Output File Header
==============================================================================================================
created_by                              Short description of the tool that created this file
creation_date                           Date when the file was created
file_description                        Short description of the content of the file (e.g. simulation name)
file_type                               Type of the file (output)
major_version                           Major version of the file definition (1)
minor_version                           Minor version of the file definition (0)
-------------------------------------------------------------------------------------------------------------
host_names                              List of hosts (computer names) the simulation was executed on
number_of_cpu_cores                     Number of CPU cores used for the simulation
data_loading_phase_execution_time       Time taken to load data from the file
pre-processing_phase_execution_time     Time taken to pre-process data
simulation_phase_execution_time         Time taken to run the simulation
post-processing_phase_execution_time    Time taken to complete the post-processing phase
total_execution_time                    Total execution time
peak_core_memory_in_use                 Peak memory required per core during the simulation
total_memory_in_use Total               Peak memory in use
==============================================================================================================
 \endverbatim
 *
 *
 *
 * All input and output parameters are stored as three dimensional datasets with dimension sizes designed by <tt>(Nx, Ny, Nz)</tt>. In order to support scalars and 1D and 2D arrays,
 * the unused dimensions are set to 1 . For example, scalar variables are stored with a dimension size of <tt>(1,1,1)</tt>, 1D vectors oriented in y-direction are stored with a dimension
 * size of <tt>(1, Ny, 1)</tt>, and so on. If the dataset stores a complex variable, the real and
 * imaginary parts are stored in an interleaved layout and the lowest used dimension size is
 * doubled (i.e., Nx for a 3D matrix, Ny for a 1D vector oriented in the y-direction). The
 * datasets are physically stored in row-major order (in contrast to column-major order used
 * by MATLAB) using either the <tt> `H5T_IEEE_F32LE' </tt> data type for floating point datasets or
 * <tt> `H5T_STD_U64LE' </tt> for integer based datasets.
 *
 * In order to enable compression of big datasets (3D variables, output time series), selected
 * datasets are not stored as monolithic blocks but broken into chunks that are compressed
 * by the ZIP library and stored separately. The chunk size is defined as follows:
 *
 * \li <tt> (Nx, Ny, 1) </tt> in the case of 3D variables (one 2D slab).
 * \li <tt> (N_sensor_points, 1, 1) </tt> in the case of the output time series (one time step of the simulation).
 *
 * All datasets have two attributes that specify the content of the dataset. The \c `data_type'
 * attribute specifies the data type of the dataset. The admissible values are either \c `float' or
 * \c `long'. The \c `domain_type' attribute specifies the domain of the dataset. The admissible
 * values are either \c `real' for the real domain or \c `complex' for the complex domain. The
 * C++ code reads these attributes and checks their values.
 *
 *
 *
 * \verbatim
==============================================================================================================
                                        Input File Datasets
==============================================================================================================
Name                            Size           Data type       Domain Type      Condition of Presence
==============================================================================================================
  1. Simulation Flags
--------------------------------------------------------------------------------------------------------------
  ux_source_flag                (1, 1, 1)       long            real
  uy_source_flag                (1, 1, 1)       long            real
  uz_source_flag                (1, 1, 1)       long            real
  p_source_flag                 (1, 1, 1)       long            real
  p0_source_flag                (1, 1, 1)       long            real
  transducer_source_flag        (1, 1, 1)       long            real
  nonuniform_grid_flag          (1, 1, 1)       long            real            must be set to 0
  nonlinear_flag                (1, 1, 1)       long            real
  absorbing_flag                (1, 1, 1)       long            real
--------------------------------------------------------------------------------------------------------------
  2. Grid Properties
--------------------------------------------------------------------------------------------------------------
  Nx                            (1, 1, 1)       long            real
  Ny                            (1, 1, 1)       long            real
  Nz                            (1, 1, 1)       long            real
  Nt                            (1, 1, 1)       long            real
  t_index                       (1, 1, 1)       long            real            File version 1.1
  dt                            (1, 1, 1)       float           real
  dx                            (1, 1, 1)       float           real
  dy                            (1, 1, 1)       float           real
  dz                            (1, 1, 1)       float           real
--------------------------------------------------------------------------------------------------------------
  3 Medium Properties
--------------------------------------------------------------------------------------------------------------
  3.1 Regular Medium Properties

  rho0                          (Nx, Ny, Nz)    float         real              heterogenous
  		                         	(1, 1, 1)       float         real              homogenous
  rho0_sgx                      (Nx, Ny, Nz)    float         real              heterogenous
                                (1, 1, 1)       float         real              homogenous
  rho0_sgy                      (Nx, Ny, Nz)    float         real              heterogenous
                                (1, 1, 1)       float         real              homogenous
  rho0_sgz                      (Nx, Ny, Nz)    float         real              heterogenous
				                        (1, 1, 1)       float         real              homogenous
  c0                            (Nx, Ny, Nz)    float         real              heterogenous
  		       		                (1, 1, 1)       float         real              homogenous
  c_ref                         (1, 1, 1)       float         real

  3.2 Nonlinear Medium Properties (defined if (nonlinear_flag == 1))

  BonA                          (Nx, Ny, Nz)    float         real              heterogenous
                                (1, 1, 1)       float         real              homogenous

  3.3 Absorbing Medium Properties (defined if (absorbing_flag == 1))

  alpha_coef                    (Nx, Ny, Nz)    float         real              heterogenous
                                (1, 1, 1)       float         real              homogenous
  alpha_power                   (1, 1, 1)       float         real
--------------------------------------------------------------------------------------------------------------
  4. Sensor Variables
--------------------------------------------------------------------------------------------------------------
  sensor_mask_type              (1, 1, 1)       long          real              File version 1.1 (0 = index, 1 = corners)
  sensor_mask_index             (Nsens, 1, 1)   long          real              File version 1.0 always, File version 1.1 if sensor_mask_type == 0
  sensor_mask_corners           (Ncubes, 6, 1)  long          real              File version 1.1, if sensor_mask_type == 1
--------------------------------------------------------------------------------------------------------------
  5 Source Properties
--------------------------------------------------------------------------------------------------------------
  5.1 Velocity Source Terms (defined if (ux_source_flag == 1 || uy_source_flag == 1 || uz_source_flag == 1))

  u_source_mode                 (1, 1, 1)          long       real
  u_source_many                 (1, 1, 1)          long       real
  u_source_index                (Nsrc, 1, 1)       long       real
  ux_source_input               (1, Nt_src, 1)     float      real              u_source_many == 0
                                (Nsrc, Nt_src, 1)  float      real              u_source_many == 1
  uy_source_input               (1, Nt_src,  1)    float      real              u_source_many == 0
                                (Nsrc, Nt_src, 1)  float      real              u_source_many == 1
  uz_source_input               (1, Nt_src, 1)     float      real              u_source_many == 0
                                (Nt_src, Nsrc, 1)  float      real              u_source_many == 1

  5.2 Pressure Source Terms (defined if p_source_flag == 1))

  p_source_mode                 (1, 1, 1)          long       real
  p_source_many                 (1, 1, 1)          long       real
  p_source_index                (Nsrc, 1, 1)       long       real
  p_source_input                (Nsrc, Nt_src, 1)  float      real              p_source_many == 1
                                (1, Nt_src, 1)     float      real              p_source_many == 0

  5.3 Transducer Source Terms (defined if (transducer_source_flag == 1))

  u_source_index        	      (Nsrc, 1, 1)       long       real
  transducer_source_input       (Nt_src, 1, 1)     float      real
  delay_mask			              (Nsrc, 1, 1)       float      real

  5.4 IVP Source Terms (defined if ( p0_source_flag ==1)

  p0_source_input               (Nx, Ny, Nz)        float     real

--------------------------------------------------------------------------------------------------------------
  6. K-space and Shift Variables
--------------------------------------------------------------------------------------------------------------
  ddx_k_shift_pos_r     	      (Nx/2 + 1, 1, 1)  float       complex
  ddx_k_shift_neg_r             (Nx/2 + 1, 1, 1)  float       complex
  ddy_k_shift_pos               (1, Ny, 1)        float       complex
  ddy_k_shift_neg               (1, Ny, 1)        float       complex
  ddz_k_shift_pos               (1, 1, Nz)        float       complex
  ddz_k_shift_neg               (1, 1, Nz)        float       complex
--------------------------------------------------------------------------------------------------------------
  7. PML Variables
--------------------------------------------------------------------------------------------------------------
  pml_x_size                    (1, 1, 1)       long          real
  pml_y_size                    (1, 1, 1)       long          real
  pml_z_size                    (1, 1, 1)       long          real
  pml_x_alpha                   (1, 1, 1)       float         real
  pml_y_alpha                   (1, 1, 1)       float         real
  pml_z_alpha                   (1, 1, 1)       float         real

  pml_x                         (Nx, 1, 1)      float         real
  pml_x_sgx                     (Nx, 1, 1)      float         real
  pml_y                         (1, Ny, 1)      float         real
  pml_y_sgy                     (1, Ny, 1)      float         real
  pml_z                         (1, 1, Nz)      float         real
  pml_z_sgz                     (1, 1, Nz)      float         real
==============================================================================================================
 \endverbatim



 \verbatim
==============================================================================================================
                                        Output File Datasets
==============================================================================================================
Name                            Size           Data type        Domain Type     Condition of Presence
==============================================================================================================
  1. Simulation Flags
--------------------------------------------------------------------------------------------------------------
  ux_source_flag                (1, 1, 1)       long          real
  uy_source_flag                (1, 1, 1)       long          real
  uz_source_flag                (1, 1, 1)       long          real
  p_source_flag                 (1, 1, 1)       long          real
  p0_source_flag                (1, 1, 1)       long          real
  transducer_source_flag        (1, 1, 1)       long          real
  nonuniform_grid_flag          (1, 1, 1)       long          real
  nonlinear_flag		            (1, 1, 1)       long          real
  absorbing_flag		            (1, 1, 1)       long          real
  u_source_mode                 (1, 1, 1)       long          real              if u_source
  u_source_many                 (1, 1, 1)       long          real              if u_source
  p_source_mode                 (1, 1, 1)       long          real              if p_source
  p_source_many                 (1, 1, 1)       long          real              if p_source

--------------------------------------------------------------------------------------------------------------
  2. Grid Properties
--------------------------------------------------------------------------------------------------------------
  Nx                            (1, 1, 1)       long          real
  Ny                            (1, 1, 1)       long          real
  Nz                            (1, 1, 1)       long          real
  Nt                            (1, 1, 1)       long          real
  dt                            (1, 1, 1)       float         real
  dx                            (1, 1, 1)       float         real
  dy                            (1, 1, 1)       float         real
  dz                            (1, 1, 1)       float         real
-------------------------------------------------------------------------------------------------------------
  3. PML Variables
--------------------------------------------------------------------------------------------------------------
  pml_x_size                    (1, 1, 1)       long          real
  pml_y_size                    (1, 1, 1)       long          real
  pml_z_size                    (1, 1, 1)       long          real
  pml_x_alpha                   (1, 1, 1)       float         real
  pml_y_alpha                   (1, 1, 1)       float         real
  pml_z_alpha                   (1, 1, 1)       float         real

  pml_x                         (Nx, 1, 1)      float         real
  pml_x_sgx                     (Nx, 1, 1)      float         real
  pml_y                         (1, Ny, 1)      float         real
  pml_y_sgy                     (1, Ny, 1)      float         real
  pml_z                         (1, 1, Nz)      float         real
  pml_z_sgz                     (1, 1, Nz)      float         real
--------------------------------------------------------------------------------------------------------------
  4. Sensor Variables (present if --copy_sensor_mask)
--------------------------------------------------------------------------------------------------------------
  sensor_mask_type              (1, 1, 1)       long          real              File version 1.1 and --copy_sensor_mask
  sensor_mask_index             (Nsens, 1, 1)   long          real              File version 1.1 and if sensor_mask_type == 0
  sensor_mask_corners           (Ncubes, 6, 1)  long          real              File version 1.1 and if sensor_mask_type == 1
--------------------------------------------------------------------------------------------------------------
  5a. Simulation Results: if sensor_mask_type == 0 (index), or File version == 1.0
--------------------------------------------------------------------------------------------------------------

  p                             (Nsens, Nt - s, 1) float      real              -p or --p_raw
  p_rms                         (Nsens, 1, 1)      float      real              --p_rms
  p_max                         (Nsens, 1, 1)      float      real              --p_max
  p_min                         (Nsens, 1, 1)      float      real              --p_min
  p_max_all                     (Nx, Ny, Nz)       float      real              --p_max_all
  p_min_all                     (Nx, Ny, Nz)       float      real              --p_min_all
  p_final                       (Nx, Ny, Nz)       float      real              --p_final


  ux                            (Nsens, Nt - s, 1) float      real              -u or --u_raw
  uy                            (Nsens, Nt - s, 1) float      real              -u or --u_raw
  uz                            (Nsens, Nt - s, 1) float      real              -u or --u_raw

  ux_non_staggered              (Nsens, Nt - s, 1) float      real              --u_non_staggered_raw
  uy_non_staggered              (Nsens, Nt - s, 1) float      real              --u_non_staggered_raw
  uz_non_staggered              (Nsens, Nt - s, 1) float      real              --u_non_staggered_raw

  ux_rms                        (Nsens, 1, 1)      float      real              --u_rms
  uy_rms                        (Nsens, 1, 1)      float      real              --u_rms
  uz_rms                        (Nsens, 1, 1)      float      real              --u_rms

  ux_max                        (Nsens, 1, 1)      float      real              --u_max
  uy_max                        (Nsens, 1, 1)      float      real              --u_max
  uz_max                        (Nsens, 1, 1)      float      real              --u_max

  ux_min                        (Nsens, 1, 1)      float      real              --u_min
  uy_min                        (Nsens, 1, 1)      float      real              --u_min
  uz_min                        (Nsens, 1, 1)      float      real              --u_min

  ux_max_all                    (Nx, Ny, Nz)       float      real              --u_max_all
  uy_max_all                    (Nx, Ny, Nz)       float      real              --u_max_all
  uz_max_all                    (Nx, Ny, Nz)       float      real              --u_max_all

  ux_min_all                    (Nx, Ny, Nz)       float      real              --u_min_all
  uy_min_all                    (Nx, Ny, Nz)       float      real              --u_min_all
  uz_min_all                    (Nx, Ny, Nz)       float      real              --u_min_all

  ux_final                      (Nx, Ny, Nz)       float      real              --u_final
  uy_final                      (Nx, Ny, Nz)       float      real              --u_final
  uz_final                      (Nx, Ny, Nz)       float      real              --u_final


 --------------------------------------------------------------------------------------------------------------
  5b. Simulation Results: if sensor_mask_type == 1 (corners) and or File version == 1.1
--------------------------------------------------------------------------------------------------------------

  /p                            group of datasets, one per cuboid                -p or --p_raw
  /p/1                          (Cx, Cy, Cz, Nt-s) float      real              1st sampled cuboid
  /p/2                          (Cx, Cy, Cz, Nt-s) float      real              2nd sampled cuboid, etc.

  /p_rms                        group of datasets, one per cuboid               --p_rms
  /p_rms/1                      (Cx, Cy, Cz, Nt-s) float      real              1st sampled cuboid

  /p_max                        group of datasets, one per cuboid               --p_max
  /p_max/1                      (Cx, Cy, Cz, Nt-s) float      real              1st sampled cuboid

  /p_min                        group of datasets, one per cuboid               --p_min
  /p_min/1                      (Cx, Cy, Cz, Nt-s) float      real              1st sampled cuboid

  p_max_all                     (Nx, Ny, Nz)       float      real              --p_max_all
  p_min_all                     (Nx, Ny, Nz)       float      real              --p_min_all
  p_final                       (Nx, Ny, Nz)       float      real              --p_final


  /ux                           group of datasets, one per cuboid               -u or --u_raw
  /ux/1                         (Cx, Cy, Cz, Nt-s) float     real               1st sampled cuboid
  /uy                           group of datasets, one per cuboid               -u or --u_raw
  /uy/1                         (Cx, Cy, Cz, Nt-s) float     real               1st sampled cuboid
  /uz                           group of datasets, one per cuboid               -u or --u_raw
  /uz/1                         (Cx, Cy, Cz, Nt-s) float     real               1st sampled cuboid

  /ux_non_staggered              group of datasets, one per cuboid               --u_non_staggered_raw
  /ux_non_staggered/1            (Cx, Cy, Cz, Nt-s) float     real               1st sampled cuboid
  /uy_non_staggered              group of datasets, one per cuboid               --u_non_staggered_raw
  /uy_non_staggered/1            (Cx, Cy, Cz, Nt-s) float     real               1st sampled cuboid
  /uz_non_staggered              group of datasets, one per cuboid               --u_non_staggered_raw
  /uz_non_staggered/1            (Cx, Cy, Cz, Nt-s) float     real               1st sampled cuboid

  /ux_rms                       group of datasets, one per cuboid               --u_rms
  /ux_rms/1                     (Cx, Cy, Cz, Nt-s) float     real               1st sampled cuboid
  /uy_rms                       group of datasets, one per cuboid               --u_rms
  /uy_rms/1                     (Cx, Cy, Cz, Nt-s) float     real               1st sampled cuboid
  /uz_rms                       group of datasets, one per cuboid               --u_rms
  /uy_rms/1                     (Cx, Cy, Cz, Nt-s) float     real               1st sampled cuboid

  /ux_max                       group of datasets, one per cuboid               --u_max
  /ux_max/1                     (Cx, Cy, Cz, Nt-s) float     real               1st sampled cuboid
  /uy_max                       group of datasets, one per cuboid               --u_max
  /ux_max/1                     (Cx, Cy, Cz, Nt-s) float     real               1st sampled cuboid
  /uz_max                       group of datasets, one per cuboid               --u_max
  /ux_max/1                     (Cx, Cy, Cz, Nt-s) float     real               1st sampled cuboid

  /ux_min                       group of datasets, one per cuboid               --u_min
  /ux_min/1                     (Cx, Cy, Cz, Nt-s) float     real               1st sampled cuboid
  /uy_min                       group of datasets, one per cuboid               --u_min
  /ux_min/1                     (Cx, Cy, Cz, Nt-s) float     real               1st sampled cuboid
  /uz_min                       group of datasets, one per cuboid               --u_min
  /ux_min/1                     (Cx, Cy, Cz, Nt-s) float     real               1st sampled cuboid

  ux_max_all                    (Nx, Ny, Nz)       float      real              --u_max_all
  uy_max_all                    (Nx, Ny, Nz)       float      real              --u_max_all
  uz_max_all                    (Nx, Ny, Nz)       float      real              --u_max_all

  ux_min_all                    (Nx, Ny, Nz)       float      real              --u_min_all
  uy_min_all                    (Nx, Ny, Nz)       float      real              --u_min_all
  uz_min_all                    (Nx, Ny, Nz)       float      real              --u_min_all

  ux_final                      (Nx, Ny, Nz)       float      real              --u_final
  uy_final                      (Nx, Ny, Nz)       float      real              --u_final
  uz_final                      (Nx, Ny, Nz)       float      real              --u_final
==============================================================================================================
\endverbatim
 *
 *
 *
 */



#include <cstdlib>
#include <iostream>
#include <exception>

#ifdef _OPENMP
  #include <omp.h>
#endif

#include <KSpaceSolver/KSpaceFirstOrder3DSolver.h>


using namespace std;


/// separator
static const char * FMT_SmallSeparator = "--------------------------------\n";

/**
 *  The main function
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char** argv)
{

  // Create K-Space solver
  TKSpaceFirstOrder3DSolver KSpaceSolver;

 // print header
  fprintf(stdout,"%s",FMT_SmallSeparator);
  fprintf(stdout,"  %s\n",KSpaceSolver.GetCodeName().c_str());
  fprintf(stdout,"%s",FMT_SmallSeparator);

  // Create parameters and parse command line
  TParameters* Parameters = TParameters::GetInstance();

  Parameters->ParseCommandLine(argc,argv);
  if (Parameters->IsVersion())
  {
    KSpaceSolver.PrintFullNameCodeAndLicense(stdout);
    return 0;
  }


  // set number of threads and bind them to cores
  #ifdef _OPENMP
    KSpaceSolver.SetProcessorAffinity();
    omp_set_num_threads(Parameters->GetNumberOfThreads());
  #endif


  fprintf(stdout, "Number of CPU threads:    %6ld\n", Parameters->GetNumberOfThreads());
  KSpaceSolver.PrintParametersOfSimulation(stdout);


  fprintf(stdout,"%s",FMT_SmallSeparator);
  fprintf(stdout,"........ Initialization ........\n");
  fprintf(stdout,"Memory allocation ..........");
  fflush(stdout);


  // allocate memory
  try
  {
    KSpaceSolver.AllocateMemory();
  }
  catch (exception e)
  {
    fprintf(stdout, "Failed!\nK-Wave panic: Not enough memory to run this simulation!\n%s\n",e.what());
    fprintf(stderr, "K-Wave panic: Not enough memory to run this simulation! \n%s\n",e.what());
    return EXIT_FAILURE;
  }
  fprintf(stdout, "Done\n");


  // Load data from disk
  fprintf(stdout, "Data loading................");
  fflush(stdout);
  try
  {
    KSpaceSolver.LoadInputData();
  }catch (ios::failure e)
  {
    fprintf(stdout, "Failed!\nK-Wave panic: Data loading was not successful!\n%s\n",e.what());
    fprintf(stderr, "K-Wave panic: Data loading was not successful! \n%s\n",e.what());
    return EXIT_FAILURE;
  }
  fprintf(stdout, "Done\n");

  fprintf(stdout,"Elapsed time:          %8.2fs\n",KSpaceSolver.GetDataLoadTime());

  if (Parameters->Get_t_index() > 0)
  {
    fprintf(stdout, "Recovered from t_index: %8ld\n", Parameters->Get_t_index());
  }

  // start computation
  fprintf(stdout,"%s",FMT_SmallSeparator);
  fprintf(stdout, ".......... Computation .........\n");

  KSpaceSolver.Compute();

  fprintf(stdout,"%s",FMT_SmallSeparator);
  fprintf(stdout, "............ Summary ...........\n");
  fprintf(stdout, "Peak memory in use:   %8ldMB\n",KSpaceSolver.ShowMemoryUsageInMB());
  if (KSpaceSolver.GetCumulatedTotalTime() != KSpaceSolver.GetTotalTime())
  {
    fprintf(stdout,"This leg execution time:%7.2fs\n",KSpaceSolver.GetTotalTime());
  }
  fprintf(stdout, "Total execution time:  %8.2fs\n",KSpaceSolver.GetCumulatedTotalTime());


  fprintf(stdout,"%s",FMT_SmallSeparator);
  fprintf(stdout,"       End of computation \n");
  fprintf(stdout,"%s",FMT_SmallSeparator);

  return  EXIT_SUCCESS;
}// end of main
//------------------------------------------------------------------------------

