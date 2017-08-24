/**
 * @file        ComplexMatrix.cpp
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The implementation file with the class for complex matrices.
 *
 * @version     kspaceFirstOrder3D 2.16
 *
 * @date        11 July      2011, 14:02 (created) \n
 *              24 August    2017, 14:42 (revised)
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

#include <iostream>

#include <MatrixClasses/ComplexMatrix.h>

#include <Utils/ErrorMessages.h>

//----------------------------------------------------------------------------//
//                              Constants                                     //
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//                              Definitions                                   //
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
//                              Implementation                                //
//                              public methods                                //
//----------------------------------------------------------------------------//

/**
 * Constructor.
 * @param [in] DimensionSizes - Dimension sizes
 */

TComplexMatrix::TComplexMatrix(const DimensionSizes & DimensionSizes)
                      : TBaseFloatMatrix()
{
  InitDimensions(DimensionSizes);

  AllocateMemory();
} // end of TComplexMatrixData
//-----------------------------------------------------------------------------


/**
 * Read data from HDF5 file (do some basic checks). Only from the root group.
 * \throw ios::failure when there is a problem
 *
 * @param [in] HDF5_File   - HDF5 file
 * @param [in] MatrixName  - HDF5 dataset name
 */
void TComplexMatrix::ReadDataFromHDF5File(Hdf5File & HDF5_File,
                                          const char * MatrixName)
{
  // check data type
  if (HDF5_File.readMatrixDataType(HDF5_File.getRootGroup(), MatrixName) != Hdf5File::MatrixDataType::kFloat)
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage, kErrFmtMatrixNotFloat, MatrixName);
    throw ios::failure(ErrorMessage);
  }

  // check domain type
  if (HDF5_File.readMatrixDomainType(HDF5_File.getRootGroup(), MatrixName) != Hdf5File::MatrixDomainType::kComplex)
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage, kErrFmtMatrixNotComplex, MatrixName);
    throw ios::failure(ErrorMessage);
  }

  // Initialise dimensions
  DimensionSizes ComplexDims = pDimensionSizes;
  ComplexDims.nx = 2 * ComplexDims.nx;

  // Read data from the file
  HDF5_File.readCompleteDataset(HDF5_File.getRootGroup(),
                                MatrixName,
                                ComplexDims,
                                pMatrixData);
}// end of LoadDataFromMatlabFile
//------------------------------------------------------------------------------

/**
 * Write data to HDF5 file (only from the root group).
 * \throw ios::failure an exception what the operation fails
 *
 * @param [in] HDF5_File             - HDF5 file handle
 * @param [in] MatrixName            - HDF5 dataset name
 * @param [in] CompressionLevel      - Compression level for the dataset
 */
void TComplexMatrix::WriteDataToHDF5File(Hdf5File & HDF5_File,
                                         const char * MatrixName,
                                         const size_t CompressionLevel)
{
  // set dimensions and chunks
  DimensionSizes ComplexDims = pDimensionSizes;
  ComplexDims.nx = 2 * ComplexDims.nx;

  DimensionSizes Chunks = ComplexDims;
  ComplexDims.nz = 1;

  // create a dataset
  hid_t HDF5_Dataset_id = HDF5_File.createDataset(HDF5_File.getRootGroup(),
                                                  MatrixName,
                                                  ComplexDims,
                                                  Chunks,
                                                  Hdf5File::MatrixDataType::kFloat,
                                                  CompressionLevel);
 // Write write the matrix at once.
  HDF5_File.writeHyperSlab(HDF5_Dataset_id,
                           DimensionSizes(0, 0, 0),
                           pDimensionSizes,
                           pMatrixData);
  HDF5_File.closeDataset(HDF5_Dataset_id);

 // Write data and domain type
  HDF5_File.writeMatrixDataType(HDF5_File.getRootGroup(),
                                MatrixName,
                                Hdf5File::MatrixDataType::kFloat);

  HDF5_File.writeMatrixDomainType(HDF5_File.getRootGroup(),
                                  MatrixName,
                                  Hdf5File::MatrixDomainType::kComplex);
}// end of WriteDataToHDF5File
//---------------------------------------------------------------------------




//----------------------------------------------------------------------------//
//                              Implementation                                //
//                             protected methods                              //
//----------------------------------------------------------------------------//
/**
 * Initialize matrix dimension sizes.
 * @param [in] DimensionSizes
 */
void TComplexMatrix::InitDimensions(const DimensionSizes & DimensionSizes)
{

  pDimensionSizes = DimensionSizes;

  pTotalElementCount = pDimensionSizes.nx *
                       pDimensionSizes.ny *
                       pDimensionSizes.nz;

  pDataRowSize = (pDimensionSizes.nx << 1);

  p2DDataSliceSize = (pDimensionSizes.nx *
                      pDimensionSizes.ny) << 1;

  // compute actual necessary memory sizes
  pTotalAllocatedElementCount = pTotalElementCount << 1;

}// end of InitDimensions
//------------------------------------------------------------------------------

//----------------------------------------------------------------------------//
//                              Implementation                                //
//                              private methods                               //
//----------------------------------------------------------------------------//
