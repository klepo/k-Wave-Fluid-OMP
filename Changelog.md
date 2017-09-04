## Version 2.16  (04 September 2017)  - release 1.2
  - Introduction of SC@FIT coding style.
  - Explicit calculation of nonlinear term.
  - New visual style and verbose levels
  - Better support for C++-11.
  - OpenMP 4.0 SIMD vectorization of kernels.
  - Parallel source injection.
  - Bugfix reading scalar flags.
  - Bugfix with transducer delay mask after checkpoint.
  - Bugfix with percentage report.

## Version 2.15 (02 October 2013) - release 1.1
  - removed -I output
  - added -u_non_staggered_raw
  - merged source bases for Windows, Linux, Intel and GNU compilers
  - updated doxyen documentation
  - used standard math library
  - added git hash and other useful information to the --version report
  - added wisdom import for checkpoint-restart
  - changed index matrices data types to size_t (used to be long)
  - default compression level is now 0
  - added a header to the checkpoint file
  - fixed aggregated values initialization

## Version 2.14 (21 January 2013) - release 1.0
  - converted from Mercurial to Git
  - fixed bug with -r verbose interval
  - added flag --p_min, --u_min, --p_max_all, --p_min_all, --u_max_all, --u_min_all
  - fixed missing <unistd.h> for gcc/4.7
  - use _OPENMP compiler flag for OpenMPI instead of __NO_OMP__
  - added cuboid cornes sensor mask
  - Output quantities are now managed by OutputStream.
  - Check HDF5 file version
  - added proper alignment for SSE and AVX
  - implemented checkpoint-restart functionality
    
## Version 2.13 (19 September 2012)
  - Start index begins from 1
  - p_avg and u_avg changed to p_rms, u_rms
  - added --u_final (p_final, and u_final stores entire matrices)
  - Intensity calculated base on the temporal and spatial staggered grid
  - log modified (FFTW -> FFT)
  - constants and typedefs moved into classes
  - doxygen documentation

## Version 2.12 (09 August 2012)
  - HDF5 input files introduced (HDF Sequential 1.8.9)
  - HDF output files introduced (all outputs in a single file)
  - Added parameter -c (compression level)
  - Matlab is not needed any more
  - terminal output format was changed
  - removed transducer source duplicity bug
  - Variables renaming
  - size scalars removed
  - Kappa, absorb_nabla1, absorb_nabla2,absorb_eta,absorb_tau generated at the beginning.
  - added scalars 
  - added linear and lossless case
  - new commandline parameters
  - error messages in a single file
  - all matrix names in a single file
  - created mercurial repository

## Version 2.11 (11 July 2012)
  - The framework has been rewritten. 
  - Added common ancestor for all matrices and the Matrix container
  - Dimension modification (stored logically according to Matlab)
  - All obsolete parts and unnecessary parts of the code were removed

## Version 2.10 (19 June 2012)
  - Added non-uniform grid (in case of ux_sgx added into dt./rho0_sgx)

## Version 2.9 (16 February 2012)     
  - Added new u source conditions
  - Added new p source conditions
  - Some matrices renamed
  - New parameters
  - Saving of u matrices
  - Saving of max pressure
  - Adding of parameter -f to store only last iteration of p
  - BonA is now a 3D matrix
  - Changed timing
  - Fixed bug in constructor

## Version 2.7 (8 December 2011)     
  - Direct import from Matlab data files
  - Export to a binary file
  - Modified TParameters and Command line parameters

## Version 2.4 (28 September 2011)     
  - the code was tidied up
  - Add c square
  - Affinity, better way

## Version 2.3 (22 September 2011)   
  - Affinity, first touch strategy
  - Pragma for fusion
  - Added three temp real matrices (no memory allocation during execution any more)
  - Removed PAPI interface
  - Removed TAU interface
  

## Version 2.2 (5 September 2011)  - debug version
  - Included PAPI interface
  - Included TAU interface
  - disabled C2C fftw plans
  - fixed progress calculation bug when running with fewer then 100 time steps
  - added flush to progress

## Version 2.1 (23 August 2011)  
  - Complex2Complex fft_3D (ifft_3D) replaced with Real2Complex and Complex2Real ones
  - All complex computation reduced to number of elements of Nx_r
  - Added Nx_r, dimension sizes
  - Altered progress output format 
  - added parameter -t (thead count) -p (print interval) -c (configuration_file)

## Version 2.0 (16 August 2011)
  - Added new code to handle transducer_input_signal, us_index, delay_mask,
        transducer_source
  - Added TOutputStreamMatrix maintaining results and store them on HDD in a stream fashion
  - Modified configuration.xml and the parameters organization
  - Default number of threads set to be the number of physical cores

## Version 1.1 (12 August 2011)  
  - Operation fusion - more tightly closed operations placed 
        under 1 loop cycle
  - Added TRealMatrixData::CopyFromComplexSameSizeAndNormalize (7% improvement)
  - fftwf_plan_with_nthreads() revised
  - fftwf_init_threads() bug removed

## Version 1.0 (10 August 2011)
  - basic implementation of KSpaceFirstOrder3D - non linear
  - based on KSpace3DBareBone 2.1
  
