/**
 * @file      ErrorMessages.h
 *
 * @author    Jiri Jaros \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            jarosjir@fit.vutbr.cz
 *
 * @brief     The header file containing routines for error messages and error messages common for both linux and
 *            windows version. The specific error messages are in separate files ErrorMessagesLinux.h
 *            and ErrorMessagesWindows.h
 *
 * @version   kspaceFirstOrder 2.17
 *
 * @date      30 August    2017, 11:39 (created) \n
 *            08 February  2023, 12:00 (revised)
 *
 * @copyright Copyright (C) 2019 Jiri Jaros and Bradley Treeby.
 *
 * This file is part of the C++ extension of the [k-Wave Toolbox](http://www.k-wave.org).
 *
 * This file is part of the k-Wave. k-Wave is free software: you can redistribute it and/or modify it under the terms
 * of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with k-Wave.
 * If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).
 */

#ifndef ERROR_MESSAGES_H
#define ERROR_MESSAGES_H

#ifdef __linux__
  #include <Logger/ErrorMessagesLinux.h>
#endif

// Windows build
#ifdef _WIN64
  #include <Logger/ErrorMessagesWindows.h>
#endif

//--------------------------------------------------------------------------------------------------------------------//
//-------------------------------- Common error messages for both Linux and Windows ----------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/// delimiters for linux paths
ErrorMessage kErrFmtPathDelimiters = "/\\_,.:-| ()[]{}";
/// error message - out of memory
ErrorMessage kErrFmtOutOfMemory    = "Error: Not enough CPU memory to run this simulation.";
/// Unknown error - unknown error
ErrorMessage kErrFmtUnknownError   = "Error: An unknown error happened. ";

//----------------------------------------------- HDF5 error messages ------------------------------------------------//
/// HDF5 error message
ErrorMessage kErrFmtCannotCreateFile      = "Error: File \"%s\" could not be created.";
/// HDF5 error message
ErrorMessage kErrFmtCannotRecreateFile    = "Error: Cannot recreate an opened file \"%s\".";
/// HDF5 error message
ErrorMessage kErrFmtCannotReopenFile      = "Error: Cannot reopen an opened file \"%s\".";
/// HDF5 error message
ErrorMessage kErrFmtCannotCloseFile       = "Error: File \"%s\" could not be closed.";
/// HDF5 error message
ErrorMessage kErrFmtCannotWriteDataset    = "Error: Could not write into \"%s\" dataset.";
/// HDF5 error message
ErrorMessage kErrFmtCannotReadDataset     = "Error: Could not read from \"%s\" dataset.";
/// HDF5 error message
ErrorMessage kErrFmtCannotGetFileSize     = "Error: Cannot get the size of file \"%s\".";
/// HDF5 error message
ErrorMessage kErrFmtBadDimensionSizes     = "Error: Dataset \"%s\"  has wrong dimension sizes.";
/// HDF5 error message
ErrorMessage kErrFmtFileNotOpen           = "Error: File \"%s\" was not found or could not be opened.";
/// HDF5 error message
ErrorMessage kErrFmtNotHdf5File           = "Error: File \"%s\" is not a valid HDF5 file.";
/// HDF5 error message
ErrorMessage kErrFmtCannotOpenDataset     = "Error: File \"%s\" could not open dataset \"%s\".";
/// HDF5 error message
ErrorMessage kErrFmtCannotCreateDataset   = "Error: File \"%s\" could not create dataset \"%s\".";
/// HDF5 error message
ErrorMessage kErrFmtCannotSetCompression  = "Error: File \"%s\", dataset \"%s\" could set compression level [%ld].";
/// HDF5 error message
ErrorMessage kErrFmtBadAttributeValue     = "Error: Bad attribute value: [%s,%s] = %s.";
/// HDF5 error message
ErrorMessage kErrFmtCannotWriteAttribute  = "Error: Could not write into \"%s\" attribute of \"%s\" dataset.";
/// HDF5 error message
ErrorMessage kErrFmtCannotReadAttribute   = "Error: Could not read from \"%s\" attribute of \"%s\" dataset.";
/// HDF5 error message
ErrorMessage kErrFmtCannotCreateGroup     = "Error: Could not create group \"%s\" in file \"%s\".";
/// HDF5 error message
ErrorMessage kErrFmtCannotOpenGroup       = "Error: Could not open group \"%s\" in file \"%s\".";
/// HDF5 error message
ErrorMessage kErrFmtBadInputFileType      = "Error: The input file has not a valid format.";
/// HDF5 error message
ErrorMessage kErrFmtBadOutputFileType     = "Error: The output file has not a valid format.";
/// HDF5 error message
ErrorMessage kErrFmtBadCheckpointFileType = "Error: The checkpoint file has not a valid format.";

//------------------------------------------------- Matrix Classes ---------------------------------------------------//
/// Matrix class error message
ErrorMessage kErrFmtMatrixNotFloat   = "Error: Matrix [%s] data type is not of single precision floating point.";
/// Matrix class error message
ErrorMessage kErrFmtMatrixNotReal    = "Error: Matrix [%s] domain is not real.";
/// Matrix class error message
ErrorMessage kErrFmtMatrixNotComplex = "Error: Matrix [%s] domain is not complex.";
/// Matrix class error message
ErrorMessage kErrFmtMatrixNotIndex   = "Error: Matrix [%s] data type is not unsigned long.";

//------------------------------------------------ Matrix Container --------------------------------------------------//
/// Matrix container error message
ErrorMessage kErrFmtBadMatrixType = "Error: Matrix [%s] has unknown type in the C++ code. [File, Line] : [%s,%d].";

/// Matrix container error message
ErrorMessage kErrFmtRelocationError = "Error: Matrix [%s] is being reallocated in matrix container.";

//-------------------------------------------- Command line Parameters -----------------------------------------------//
/// Command line parameters error message
ErrorMessage kErrFmtNoProgressPrintInterval = "Error: No or invalid progress print interval.";
/// Command line parameters error message
ErrorMessage kErrFmtInvalidNumberOfThreads  = "Error: No or invalid number of CPU threads.";
/// Command line parameters error message
ErrorMessage kErrFmtNoCompressionLevel      = "Error: No or invalid compression level.";
/// Command line parameters error message
ErrorMessage kErrFmtNoSamplingStartTimeStep = "Error: No or invalid collection start time step.";
/// Command line parameters error message
ErrorMessage kErrFmtNoBenchmarkTimeStep     = "Error: No or invalid number of time step to benchmark.";
/// Command line parameters error message
ErrorMessage kErrFmtNoVerboseLevel          = "Error: No or invalid verbose level.";

/// Error message - input file was not specified
ErrorMessage kErrFmtNoInputFile           = "Error: The input file was not specified.";
/// Command line parameters error message
ErrorMessage kErrFmtNoOutputFile          = "Error: The output file was not specified.";
/// Command line parameters error message
ErrorMessage kErrFmtNoCheckpointFile      = "Error: The checkpoint file was not specified.";
/// Command line parameters error message
ErrorMessage kErrFmtNoCheckpointInterval  = "Error: The checkpoint interval was not specified.";
/// Command line parameters error message
ErrorMessage kErrFmtNoCheckpointTimeSteps = "Error: The checkpoint time steps were not specified.";
/// Command line parameters error message
ErrorMessage kErrFmtNoCheckpointIntervalOrTimeSteps =
  "Error: Neither checkpoint interval in seconds nor in time steps was specified.";
/// Command line parameter error message
ErrorMessage kErrFmtUnknownParameter           = "Error: Unknown command line parameter.";
/// Command line parameter error message
ErrorMessage kErrFmtUnknownParameterOrArgument = "Error: Unknown command line parameter or missing argument.";

/// Command line parameters error message
ErrorMessage kErrFmtInvalidFrequency      = "Error: Invalid frequency value.";
/// Command line parameters error message
ErrorMessage kErrFmtInvalidPeriod         = "Error: Invalid period value.";
/// Command line parameters error message
ErrorMessage kErrFmtInvalidMOS            = "Error: Invalid multiple of overlap size value.";
/// Command line parameters error message
ErrorMessage kErrFmtInvalidHarmonics      = "Error: Invalid harmonics value.";
/// Command line parameters error message
ErrorMessage kErrFmtMissingFrequencyValue = "Error: Missing frequency value.";
/// Command line parameters error message
ErrorMessage kErrFmtMissingPeriodValue    = "Error: Missing period value.";
/// Command line parameters error message
ErrorMessage kErrFmtBadPeriodAndFrequencyValue =
  "Error: Period and frequency values are set concurrently. Set only one of them.";
/// Command line parameters error message
ErrorMessage kErrFmtCannotComputePeriodValue = "Error: Cannot compute period from frequency, missing dt dataset.";
/// Command line parameters error message
ErrorMessage kErrFmtCannotReadCompressionAttributes = "Error: Cannot read compression parameters from attributes.";

/// Command line parameters error message
ErrorMessage kErrFmtInvalidBlockSize = "Error: Invalid reading block size value.";

/// Command line parameters error message
ErrorMessage kErrFmtInvalidPostProcessingFlag =
  "Error: Invalid store flag with combination of the only post-processing flag.";

/// Command line parameters error message
ErrorMessage kErrFmtIllegalAlphaPowerValue = "Error: Illegal value of alpha_power (must not equal to 1.0).";
/// Command line parameters error message
ErrorMessage kErrFmtIllegalSamplingStartTimeStep =
  "Error: The beginning of data sampling is out of the simulation time span <%zu, %zu>.";

/// Command line parameters error message
ErrorMessage kErrFmtBadInputFileFormat  = "Error: Incorrect input file\"%s\" format.";
/// Command line parameters error message
ErrorMessage kErrFmtBadMajorFileVersion = "Error: Incorrect major version of the HDF5 file %s (expected is %s).";
/// Command line parameters error message
ErrorMessage kErrFmtBadMinorFileVersion = "Error: Incorrect minor version of the HDF5 file %s (expected is %s).";
/// Command line parameters error message
ErrorMessage kErrFmtBadSensorMaskType   = "Error: The sensor mask type specified in the input file is not supported.";
/// Command line parameters error message
ErrorMessage kErrFmtNonStaggeredVelocityNotSupportedFileVersion =
  "Error: --u_non_staggered_raw is not supported along with the input file of the version 1.0.";
/// Command line parameters error message
ErrorMessage kErrFmtBadVelocitySourceMode =
  "Error: The velocity source mode type specified in the input file is not supported.";
/// Command line parameters error message
ErrorMessage kErrFmtBadPressureSourceMode =
  "Error: The pressure source mode specified in the input file is not supported.";

//-------------------------------------------- KSpaceFirstOrder3DSolver ----------------------------------------------//
/// KSpaceFirstOrderSolver error message
ErrorMessage kErrFmtBadCheckpointFileFormat = "Error: Incorrect checkpoint file \"%s\" format.";

/// KSpaceFirstOrderSolver error message
ErrorMessage kErrFmtBadOutputFileFormat = "Error: Incorrect output file \"%s\" format.";

/// KSpaceFirstOrderSolver error message
ErrorMessage kErrFmtCheckpointDimensionsMismatch =
  "Error: The dimensions [%ld, %ld, %ld] of the checkpoint file don't match the simulation "
  "dimensions [%ld, %ld, %ld].";

/// KSpaceFirstOrderSolver error message
ErrorMessage kErrFmtOutputDimensionsMismatch =
  "Error: The dimensions [%ld, %ld, %ld] of the output file don't match the simulation "
  "dimensions [%ld, %ld, %ld].";

//-------------------------------------------------- FFTW errors -----------------------------------------------------//
/// FFTW error message
ErrorMessage kErrFmtCreateR2CFftPlanND  = "creating plan for ND real-to-complex fft.";
/// FFTW error message
ErrorMessage kErrFmtCreateC2RFftPlanND  = "creating plan for ND complex-to-real fft.";
/// FFTW error message
ErrorMessage kErrFmtCreateR2CFftPlan1DX = "creating for 1D real-to-complex fft plan in X direction.";
/// FFTW error message
ErrorMessage kErrFmtCreateR2CFftPlan1DY = "creating for 1D real-to-complex fft plan in Y direction.";
/// FFTW error message
ErrorMessage kErrFmtCreateR2CFftPlan1DZ = "creating for 1D real-to-complex fft plan in Z direction.";
/// FFTW error message
ErrorMessage kErrFmtCreateC2RFftPlan1DX = "creating for 1D complex-to-real fft plan in X direction.";
/// FFTW error message
ErrorMessage kErrFmtCreateC2RFftPlan1DY = "creating for 1D complex-to-real fft plan in Y direction.";
/// FFTW error message
ErrorMessage kErrFmtCreateC2RFftPlan1DZ = "creating 1D complex-to-real fft plan in Z direction.";

/// FFTW error message
ErrorMessage kErrFmtCannotCallR2CFftPlan1DZfor2D = "createR2CFftPlan1DZ cannot be called in 2D simulations";
/// FFTW error message
ErrorMessage kErrFmtCannotCallC2RFftPlan1DZfor2D = "createC2RFftPlan1DZ cannot be called in 2D simulations";

/// FFTW error message
ErrorMessage kErrFmtDestroyR2CFftPlan3D  = "destroying plan for 3D real-to-complex fft.";
/// FFTW error message
ErrorMessage kErrFmtDestroyC2RFftPlan3D  = "destroying plan for 3D complex-to-real fft.";
/// FFTW error message
ErrorMessage kErrFmtDestroyR2CFftPlan1DX = "destroying for 1D real-to-complex fft plan in X direction.";
/// FFTW error message
ErrorMessage kErrFmtDestroyR2CFftPlan1DY = "destroying for 1D real-to-complex fft plan in Y direction.";
/// FFTW error message
ErrorMessage kErrFmtDestroyR2CFftPlan1DZ = "destroying for 1D real-to-complex fft plan in Z direction.";
/// FFTW error message
ErrorMessage kErrFmtDestroyC2RFftPlan1DX = "destroying for 1D complex-to-real fft plan in X direction.";
/// FFTW error message
ErrorMessage kErrFmtDestroyC2RFftPlan1DY = "destroying for 1D complex-to-real fft plan in Y direction.";
/// FFTW error message
ErrorMessage kErrFmtDestroyC2RFftPlan1DZ = "destroying 1D complex-to-real fft plan in Z direction.";

/// FFTW error message
ErrorMessage kErrFmtExecuteR2CFftPlan3D  = "executing plan for 3D real-to-complex fft.";
/// FFTW error message
ErrorMessage kErrFmtExecuteC2RFftPlan3D  = "executing plan for 3D complex-to-real fft.";
/// FFTW error message
ErrorMessage kErrFmtExecuteR2CFftPlan1DX = "executing for 1D real-to-complex fft plan in X direction.";
/// FFTW error message
ErrorMessage kErrFmtExecuteR2CFftPlan1DY = "executing for 1D real-to-complex fft plan in Y direction.";
/// FFTW error message
ErrorMessage kErrFmtExecuteR2CFftPlan1DZ = "executing for 1D real-to-complex fft plan in Z direction.";
/// FFTW error message
ErrorMessage kErrFmtExecuteC2RFftPlan1DX = "executing for 1D complex-to-real fft plan in X direction.";
/// FFTW error message
ErrorMessage kErrFmtExecuteC2RFftPlan1DY = "executing for 1D complex-to-real fft plan in Y direction.";
/// FFTW error message
ErrorMessage kErrFmtExecuteC2RFftPlan1DZ = "executing 1D complex-to-real fft plan in Z direction.";

/// FFTW error message
ErrorMessage kErrFmtFftWisdomNotExported = "Warning: Wisdom could not be exported.";
/// FFTW error message
ErrorMessage ErrFmtFftWisdomNotImported  = "Warning: Wisdom could not be imported.";

#endif /* ERROR_MESSAGES_H */
