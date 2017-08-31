/**
 * @file        ErrorMessages.h
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The header file containing all error messages of the project.
 *
 * @version     kspaceFirstOrder3D 2.16
 *
 * @date        09 August    2011, 12:34      (created) \n
 *              24 August    2017, 09:23      (revised)
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
 * along with k-Wave. If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef ERROR_MESSAGES_H
#define ERROR_MESSAGES_H

/**
 * @brief   Datatype for error messages.
 * @details Datatype for error messages.
 */
using ErrorMessage = const char*;

//----------------------------------------------- HDF5 error messages ------------------------------------------------//

/// HDF5 error message
ErrorMessage const kErrFmtCannotCreateFile
 = "Error: File \"%s\" could not be created!\n";
/// HDF5 error message
ErrorMessage const kErrFmtCannotRecreateFile
  = "Error: Cannot recreate an opened file \"%s\"!\n";
/// HDF5 error message
ErrorMessage const kErrFmtCannotReopenFile
  = "Error: Cannot reopen an opened file \"%s\"!\n";
/// HDF5 error message
ErrorMessage const kErrFmtCannotCloseFile
  = "Error: File \"%s\" could not be closed!\n";
/// HDF5 error message
ErrorMessage const kErrFmtCannotWriteDataset
  = "Error: Could not write into \"%s\" dataset!\n";
/// HDF5 error message
ErrorMessage const kErrFmtCannotReadDataset
  = "Error: Could not read from \"%s\" dataset!\n";
/// HDF5 error message
ErrorMessage const kErrFmtBadDimensionSizes
  = "Error: Dataset \"%s\"  has wrong dimension sizes!\n";
/// HDF5 error message
ErrorMessage const kErrFmtFileNotOpen
  = "Error: File \"%s\" could not be opened!\n";
/// HDF5 error message
ErrorMessage const kErrFmtNotHdf5File
  = "Error: File \"%s\" is not a valid HDF5 file!\n";
/// HDF5 error message
ErrorMessage const kErrFmtCannotOpenDataset
  = "Error: File \"%s\" could not open dataset \"%s\"!\n";
/// HDF5 error message
ErrorMessage const kErrFmtCannotSetCompression
  = "Error: File \"%s\", dataset \"%s\" could set compression level [%ld]!\n";
/// HDF5 error message
ErrorMessage const kErrFmtBadAttributeValue
  = "Error: Bad attribute value: [%s,%s] = %s";
/// HDF5 error message
ErrorMessage const kErrFmtCannotWriteAttribute
  = "Error: Could not write into \"%s\" attribute of \"%s\" dataset!\n";
/// HDF5 error message
ErrorMessage const kErrFmtCannotReadAttribute
  = "Error: Could not read from \"%s\" attribute of \"%s\" dataset!\n";
/// HDF5 error message
ErrorMessage const kErrFmtCannotCreateGroup
  = "Error: Could not create group \"%s\" in file \"%s\"!\n";
/// HDF5 error message
ErrorMessage const kErrFmtCannotOpenGroup
  = "Error: Could not open group \"%s\" in file \"%s\"!\n";
/// HDF5 error message
ErrorMessage const kErrFmtBadInputFileType
  = "Error: The input file has not a valid format!\n";
/// HDF5 error message
ErrorMessage const kErrFmtBadOutputFIleType
  = "Error: The output file has not a valid format!\n";
/// HDF5 error message
ErrorMessage const kErrFmtBadCheckpointFileType
  = "Error: The checkpoint file has not a valid format!\n";

//------------------------------------------------- Matrix Classes ---------------------------------------------------//

/// Matrix class error message
ErrorMessage const  kErrFmtNotEnoughMemory
  = "Error: Class %s: Memory allocation failed: Not Enough Memory\n";
/// Matrix class error message
ErrorMessage const  kErrFmtMatrixNotFloat
  = "Error: Matrix [%s] data type is not of single precision floating point!\n";
/// Matrix class error message
ErrorMessage const  kErrFmtMatrixNotReal
  = "Error: Matrix [%s] domain is not real!\n";
/// Matrix class error message
ErrorMessage const  kErrFmtMatrixNotComplex
  = "Error: Matrix [%s] domain is not complex!\n";
/// Matrix class error message
ErrorMessage const  kErrFmtMatrixNotIndex
  = "Error: Matrix [%s] data type is not unsigned long (uint64_t)!\n";

//------------------------------------------------ Matrix Container --------------------------------------------------//
/// Matrix container error message
ErrorMessage const  kErrFmtBadMatrixType
  = "K-Space panic: Matrix [%s] has unknown distribution type in the C++ code!\n\
    [File, line] : [%s,%d]!\n";

/// Matrix container error message
ErrorMessage const  kErrFmtRelocationError
  = "K-Space panic: Matrix [%s] is being reallocated!\n\
    [File, line] : [%s,%d]!\n";


//-------------------------------------------- Command line Parameters -----------------------------------------------//
/// Command line parameters error message
ErrorMessage const kErrFmtNoProgressPrintInterval
  = "Command line parsing error: No or invalid verbose interval provided!\n";
/// Command line parameters error message
ErrorMessage const kErrFmtInvalidNumberOfThreads
  = "Command line parsing error: No or invalid number of CPU threads!\n";
/// Command line parameters error message
ErrorMessage const kErrFmtNoCompressionLevel
  = "Command line parsing error: No or invalid compression level!\n";
/// Command line parameters error message
ErrorMessage const kErrFmtNoSamplingStartTimeStep
  = "Command line parsing error: No or invalid collection start time step!\n";
/// Command line parameters error message
ErrorMessage const kErrFmtNoBenchmarkTimeStep
  = "Command line parsing error: No or invalid benchmark time step count!\n";
/// Command line parameters error message
ErrorMessage const kErrFmtNoInputFile
  = "Error: The input file was not specified!\n";
/// Command line parameters error message
ErrorMessage const kErrFmtNoOutputFile
  = "Error: The output file was not specified!\n";
/// Command line parameters error message
ErrorMessage const kErrFmtNoCheckpointFile
  = "Error: The checkpoint file was not specified!\n";
/// Command line parameters error message
ErrorMessage const kErrFmtNoCheckpointInterval
  = "Error: The checkpoint interval was not specified!\n";

/// Command line parameters error message
ErrorMessage const kErrFmtIllegalAlphaPowerValue
  = "Error: Illegal value of alpha_power!";
/// Command line parameters error message
ErrorMessage const kErrFmtIllegalSamplingStartTimeStep
  = "Error: The start index is out of the simulation span <%ld, %ld>!\n";
/// Command line parameters error message
ErrorMessage const kErrFmtBadInputFileFormat
  = "Error: Incorrect input file\"%s\" format!\n";
/// Command line parameters error message
ErrorMessage const kErrFmtBadMajorFileVersion
  = "Error: Incorrect major version of the HDF5 file %s (expected is %s)!\n";
/// Command line parameters error message
ErrorMessage const kErrFmtBadMinorFileVersion
  = "Error: Incorrect minor version of the HDF5 file %s (expected is %s)!\n";
/// Command line parameters error message
ErrorMessage const kErrFmtBadSensorMaskType
  = "Error: The sensor mask type specified in the input file is not supported! \n";
/// Command line parameters error message
ErrorMessage const kErrFmtNonStaggeredVelocityNotSupportedFileVersion
  = "Error: --u_non_staggered_raw is not supported along with the input file of the version 1.0! \n";

//-------------------------------------------------- FFTW Classes  ---------------------------------------------------//
/// FFTW error message
ErrorMessage const kErrFmtFftPlanNotCreated
  = "Error: The FFTW plan creation for %s failed! \n";
/// FFTW error message
ErrorMessage const kErrFmtFftInvalidPlan
  = "Error: Invalid plan for %s! \n";


//-------------------------------------------- KSpaceFirstOrder3DSolver ----------------------------------------------//
/// KSpaceFirstOrder3DSolver error message
ErrorMessage const kErrFmtBadCheckpointFileFormat
  = "Error: Incorrect checkpoint file \"%s\" format!\n";

/// KSpaceFirstOrder3DSolver error message
ErrorMessage const kErrFmtBadOutputFileFormat
  = "Error: Incorrect output file \"%s\" format!\n";

/// KSpaceFirstOrder3DSolver error message
ErrorMessage const kErrFmtCheckpointDimensionsMismatch
  = "Error: The dimensions [%ld, %ld, %ld] of the checkpoint file don't match the simulation dimensions"
    " [%ld, %ld, %ld] \n";

/// KSpaceFirstOrder3DSolver error message
ErrorMessage const kErrFmtOutputDimensionsMismatch
  = "Error: The dimensions [%ld, %ld, %ld] of the output file don't match the simulation dimensions [%ld, %ld, %ld] \n";

#endif	/* ERROR_MESSAGES_H */

