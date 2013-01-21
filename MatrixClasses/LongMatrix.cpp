/**
 * @file        LongMatrix.cpp
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au
 * 
 * @brief       The implementation file containing the class for 64b integer matrices
 * 
 * @version     kspaceFirstOrder3D 2.13
 * @date        26 July 2011, 15:16             (created) \n
 *              17 September 2012, 15:55PM      (revised)
 * 
 * @section License
 * This file is part of the C++ extension of the k-Wave Toolbox (http://www.k-wave.org).\n
 * Copyright (C) 2012 Jiri Jaros and Bradley Treeby
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

#include <MatrixClasses/LongMatrix.h>

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
 * Constructor
 * @param [in] DimensionSizes - Dimension sizes 
 */
TLongMatrix::TLongMatrix(struct TDimensionSizes DimensionSizes) : 
                TBaseLongMatrix() 
{
     
    pDimensionSizes = DimensionSizes;
       
     
    pTotalElementCount = pDimensionSizes.X *
                         pDimensionSizes.Y *
                         pDimensionSizes.Z;

    pTotalAllocatedElementCount = pTotalElementCount;
    
    pDataRowSize       = pDimensionSizes.X;

    p2DDataSliceSize   = pDimensionSizes.X *
                         pDimensionSizes.Y;                        
            
    
    AllocateMemory();
}// end of TRealMatrixData
//-----------------------------------------------------------------------------



/**
 * Read data from HDF5 file 
 * @throw ios:failure if there's an error
 * 
 * @param HDF5_File - HDF5 file handle
 * @param MatrixName  - HDF5 dataset name
 */
void TLongMatrix::ReadDataFromHDF5File(THDF5_File & HDF5_File, const char * MatrixName){
    

    if (HDF5_File.ReadMatrixDataType(MatrixName) != THDF5_File::hdf5_mdt_long){        
        char ErrorMessage[256];
        sprintf(ErrorMessage,Matrix_ERR_FMT_MatrixNotLong,MatrixName);                        
        throw ios::failure(ErrorMessage);       
    }
    
    
    if (HDF5_File.ReadMatrixDomainType(MatrixName) != THDF5_File::hdf5_mdt_real){
        char ErrorMessage[256];
        sprintf(ErrorMessage,Matrix_ERR_FMT_MatrixNotReal,MatrixName);                        
        throw ios::failure(ErrorMessage);       
    }
    
    
    HDF5_File.ReadCompleteDataset(MatrixName,pDimensionSizes,pMatrixData);
    
}// end of LoadDataFromMatlabFile
//------------------------------------------------------------------------------
    



/**
 * Recompute indeces, MATLAB -> C++ 
 */
void TLongMatrix::RecomputeIndices()
{
    
    for (size_t i = 0; i< pTotalElementCount; i++){
        pMatrixData[i]--;
    }    
    
}// end of RecomputeIndices
//------------------------------------------------------------------------------



/**
 * Write data to HDF5 file 
 * @throw ios:failure 
 * 
 * @param [in] HDF5_File - HDF5 file handle
 * @param [in] MatrixName  - HDF5 Dataset name
 * @param [in] CompressionLevel - Compression level
 */
void TLongMatrix::WriteDataToHDF5File(THDF5_File & HDF5_File, const char * MatrixName, const int CompressionLevel){
    
    
    HDF5_File.WriteCompleteDataset(MatrixName,pDimensionSizes,pMatrixData);
    
    HDF5_File.WriteMatrixDataType  (MatrixName, THDF5_File::hdf5_mdt_long);        
    HDF5_File.WriteMatrixDomainType(MatrixName, THDF5_File::hdf5_mdt_real);
    
    
}// end of WriteDataToHDF5File
//---------------------------------------------------------------------------


//----------------------------------------------------------------------------//
//                              Implementation                                //
//                             protected methods                              //
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
//                              Implementation                                //
//                              private methods                               //
//----------------------------------------------------------------------------//
