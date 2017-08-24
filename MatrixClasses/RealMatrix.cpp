/**
 * @file        RealMatrix.cpp
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The implementation file containing the class for real matrices.
 *
 * @version     kspaceFirstOrder3D 2.16
 *
 * @date        11 July      2011, 10:30 (created) \n
 *              24 August    2017, 12:21 (revised)
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
#include <string.h>

#include <MatrixClasses/RealMatrix.h>
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
TRealMatrix::TRealMatrix(const DimensionSizes & DimensionSizes)
                    : TBaseFloatMatrix()
{
  InitDimensions(DimensionSizes);

  AllocateMemory();
}// end of TRealMatrixData
//-----------------------------------------------------------------------------

/**
 * Read data data from HDF5 file (only from the root group).
 *
 * @param [in] HDF5_File  - HDF5 file
 * @param [in] MatrixName - HDF5 dataset name
 *
 * @throw ios::failure if error occurred.
 */
void TRealMatrix::ReadDataFromHDF5File(Hdf5File & HDF5_File,
                                       const char * MatrixName)
{
  // test matrix datatype
  if (HDF5_File.ReadMatrixDataType(HDF5_File.GetRootGroup(), MatrixName) != Hdf5File::hdf5_mdt_float)
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage, kErrFmtMatrixNotFloat, MatrixName);
    throw ios::failure(ErrorMessage);
  }


  if (HDF5_File.ReadMatrixDomainType(HDF5_File.GetRootGroup(), MatrixName) != Hdf5File::hdf5_mdt_real)
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage, kErrFmtMatrixNotReal, MatrixName);
    throw ios::failure(ErrorMessage);
  }

  // Read matrix
  HDF5_File.ReadCompleteDataset(HDF5_File.GetRootGroup(),
                                MatrixName,
                                pDimensionSizes,
                                pMatrixData
                                );
}// end of LoadDataFromMatlabFile
//------------------------------------------------------------------------------

/**
 * Write data to HDF5 file (only from the root group)
 *
 * @param [in] HDF5_File        - HDF5 file
 * @param [in] MatrixName       - HDF5 Matrix name
 * @param [in] CompressionLevel - Compression level
 *
 * @throw ios::failure if an error occurred
 */
void TRealMatrix::WriteDataToHDF5File(Hdf5File & HDF5_File,
                                      const char * MatrixName,
                                      const size_t CompressionLevel)
{
  DimensionSizes Chunks = pDimensionSizes;
  Chunks.nz = 1;

  //1D matrices
  if ((pDimensionSizes.ny == 1) && (pDimensionSizes.nz == 1))
  {
    // Chunk = 4MB
    if (pDimensionSizes.nx > 4 * ChunkSize_1D_4MB)
    {
      Chunks.nx = ChunkSize_1D_4MB;
    }
    else if (pDimensionSizes.nx > 4 * ChunkSize_1D_1MB)
    {
      Chunks.nx = ChunkSize_1D_1MB;
    }
    else if (pDimensionSizes.nx > 4 * ChunkSize_1D_256KB)
    {
      Chunks.nx = ChunkSize_1D_256KB;
    }
  }

  hid_t HDF5_Dataset_id = HDF5_File.CreateFloatDataset(HDF5_File.GetRootGroup(),
                                                       MatrixName,
                                                       pDimensionSizes,
                                                       Chunks,
                                                       CompressionLevel);

  HDF5_File.WriteHyperSlab(HDF5_Dataset_id,
                           DimensionSizes(0, 0, 0),
                           pDimensionSizes,
                           pMatrixData);

  HDF5_File.CloseDataset(HDF5_Dataset_id);

  // Write data and domain type
  HDF5_File.WriteMatrixDataType  (HDF5_File.GetRootGroup(),
                                  MatrixName,
                                  Hdf5File::hdf5_mdt_float);
  HDF5_File.WriteMatrixDomainType(HDF5_File.GetRootGroup(),
                                  MatrixName,
                                  Hdf5File::hdf5_mdt_real);
}// end of WriteDataToHDF5File
//------------------------------------------------------------------------------

//----------------------------------------------------------------------------//
//                              Implementation                                //
//                             protected methods                              //
//----------------------------------------------------------------------------//

/**
 * Set necessary dimensions and auxiliary variables.
 * @param [in] DimensionSizes - 3D Dimension sizes
 */
void TRealMatrix::InitDimensions(const DimensionSizes & DimensionSizes)
{
  pDimensionSizes = DimensionSizes;

  pTotalElementCount = pDimensionSizes.nx *
                       pDimensionSizes.ny *
                       pDimensionSizes.nz;

  pTotalAllocatedElementCount = pTotalElementCount;

  pDataRowSize = pDimensionSizes.nx;

  p2DDataSliceSize = pDimensionSizes.nx *
                     pDimensionSizes.ny;
}// end of SetDimensions
//------------------------------------------------------------------------------/


//----------------------------------------------------------------------------//
//                              Implementation                                //
//                              private methods                               //
//----------------------------------------------------------------------------//
