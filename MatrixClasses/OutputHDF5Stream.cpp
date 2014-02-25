/**
 * @file        OutputHDF5Stream.cpp
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au
 * 
 * @brief       The implementation file of the class saving RealMatrix data into 
 *              the output HDF5 file
 * 
 * @version     kspaceFirstOrder3D 2.14
 * @date        11 July     2012, 10:30      (created) \n
 *              18 February 2014, 15:35 (revised)
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



#include <MatrixClasses/OutputHDF5Stream.h>


using namespace std;

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
 * This method initialize the output stream by creating a HDF5 dataset and writing
 * the parameters of the dataset
 * @param [in] HDF5_File          - Handle to HDF5 output file
 * @param [in] DatasetName        - Dataset name
 * @param [in] TotalSize          - Total size of the dataset (usually sensor_size * time_period)
 * @param [in] ChunkSize          - Chunk size (usually sensor_size)
 * @param [in] CompressionLevel   - Compression level
 */
void TOutputHDF5Stream::CreateStream(THDF5_File & HDF5_File, const char * DatasetName, 
                                     const TDimensionSizes & TotalSize, const TDimensionSizes & ChunkSize, 
                                     const int CompressionLevel){
    
    this->HDF5_File = &HDF5_File;
    
    HDF5_Dataset_id = HDF5_File.CreateFloatDataset(DatasetName,TotalSize,ChunkSize, CompressionLevel);
    HDF5_File.WriteMatrixDomainType(DatasetName,THDF5_File::hdf5_mdt_real);
    HDF5_File.WriteMatrixDataType(DatasetName,THDF5_File::hdf5_mdt_float);        
    
    Position = TDimensionSizes(0,0,0);
    
    
}// end CreateStream
//------------------------------------------------------------------------------

/**
 * Close the output stream and the dataset
 * 
 */
void TOutputHDF5Stream::CloseStream(){
    
    HDF5_File->CloseDataset(HDF5_Dataset_id);
    HDF5_Dataset_id  = H5I_BADID;
    
}// end of CloseStream
//------------------------------------------------------------------------------



/**
 * Add data into the stream (usually one time step sensor data) based on Index Mask
 * 
 * @param [in] SourceMatrix - Matrix from where to pick the values
 * @param [in] Index        - Index used to pick the values
 * @param [in, out] TempBuffer - Temp buffer to make the data block contiguous
 * 
 */
void TOutputHDF5Stream::AddDataIndex(TRealMatrix& SourceMatrix, TLongMatrix& Index, float * TempBuffer){
     
     
    #pragma omp parallel for if (Index.GetTotalElementCount() > 1e6)  
    for (size_t i = 0; i < Index.GetTotalElementCount(); i++) {
        TempBuffer[i] = SourceMatrix[Index[i]];        
    }
    
    
    TDimensionSizes BlockSize(Index.GetTotalElementCount(),1,1);
    
    HDF5_File->WriteHyperSlab(HDF5_Dataset_id,Position,BlockSize,TempBuffer);
    Position.Y++;
    
}// end of AddData
//------------------------------------------------------------------------------



/**
 * Add data into the stream (usually one time step sensor data) based on Cuboid Corners Mask
 * 
 * @param [in] SourceMatrix - Matrix from where to pick the values
 * @param [in] Corners      - Cuboid corners the interior is sampled
 * @param [in, out] TempBuffer - Temp buffer to make the data block contiguous
 * 
 */

void TOutputHDF5Stream::AddDataCorners(TRealMatrix& SourceMatrix, TLongMatrix& Corners, float * TempBuffer)
{
    ///@TODO Parallel section with some load balancing feature
    
    size_t BufferIndex = 0;
    const size_t XY_Size = SourceMatrix.GetDimensionSizes().Y * SourceMatrix.GetDimensionSizes().X;
    const size_t X_Size  = SourceMatrix.GetDimensionSizes().X;
    
    
    for (size_t CuboidIdx = 0; CuboidIdx < Corners.GetDimensionSizes().Y; CuboidIdx++)        
    {
        const TDimensionSizes TopLeftCorner     = Corners.GetTopLeftCorner(CuboidIdx);
        const TDimensionSizes BottomRightCorner = Corners.GetBottomRightCorner(CuboidIdx);
        
        for (size_t z = TopLeftCorner.Z; z <= BottomRightCorner.Z; z++) 
            for (size_t y = TopLeftCorner.Y; y <= BottomRightCorner.Y; y++) 
                for (size_t x = TopLeftCorner.X; x <= BottomRightCorner.X; x++) 
                {
                   TempBuffer[BufferIndex++] = SourceMatrix[z * XY_Size + y * X_Size + x];                           
                }
    }
    
    
    TDimensionSizes BlockSize(Corners.GetTotalNumberOfElementsInAllCuboids(),1,1);
    
    HDF5_File->WriteHyperSlab(HDF5_Dataset_id,Position,BlockSize,TempBuffer);
    Position.Y++;
    
    
    
}// end of AddDataCorners   
//------------------------------------------------------------------------------

/**
 * Destructor
 * 
 */
TOutputHDF5Stream::~TOutputHDF5Stream(){
    
    if (HDF5_Dataset_id != -1) HDF5_File->CloseDataset(HDF5_Dataset_id);
    
}// end of Destructor
//------------------------------------------------------------------------------


//----------------------------------------------------------------------------//
//                              Implementation                                //
//                             protected methods                              //
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
//                              Implementation                                //
//                              private methods                               //
//----------------------------------------------------------------------------//









