/**
 * @file        OutputHDF5Stream.cpp
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au
 * 
 * @brief       The implementation file of classes responsible for storing output 
 *              quantities into the output HDF5 file
 * 
 * @version     kspaceFirstOrder3D 2.14
 * @date        11 July     2012, 10:30      (created) \n
 *              02 July     2014, 14:38      (revised)
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


#include <string.h>
#include <cmath>
#include <malloc.h>

#include <MatrixClasses/OutputHDF5Stream.h>

#include <Parameters/Parameters.h>
#include <Utils/ErrorMessages.h>


using namespace std;

//----------------------------------------------------------------------------//
//                              Constants                                     //
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
//                              Definitions                                   //
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
//                  TBaseOutputHDF5Stream implementation                      //
//                              public methods                                //
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//                  TBaseOutputHDF5Stream implementation                      //
//                              public methods                                //
//----------------------------------------------------------------------------//


/**
 * Allocate memory using a proper memory alignment.
 * @warning - This can routine is not used in the base class (should be used in 
 *            derived ones
 */
void TBaseOutputHDF5Stream::AllocateMemory()
{
  StoreBuffer = (float *) memalign(DATA_ALIGNMENT, BufferSize * sizeof (float));

  if (!StoreBuffer)
  {
    fprintf(stderr, Matrix_ERR_FMT_NotEnoughMemory, "TBaseOutputHDF5Stream");
    throw bad_alloc();
  }

  // zero the matrix
  #pragma omp parallel for if (BufferSize > 1e6)  
  for (size_t i = 0; i < BufferSize; i++)
  {
    StoreBuffer[i] = 0.0f;
  }
  
}// end of AllocateMemory
//------------------------------------------------------------------------------

/**
 * Free memory.
 * @warning - This can routine is not used in the base class (should be used in 
 *            derived ones
 */
void TBaseOutputHDF5Stream::FreeMemory()
{
  if (StoreBuffer)
  {
    free(StoreBuffer);
    StoreBuffer = NULL;
  }
}// end of FreeMemory
//------------------------------------------------------------------------------




/**
 * Apply post-processing on the buffer. It supposes the elements are independent
 * @warning - This can routine is not used in the base class (should be used in 
 *            derived ones
 */
void TBaseOutputHDF5Stream::ApplyPostProcessing()
{     
  switch (ReductionOp)
  { 
    case roNONE :
    {
      // do nothing
      break;
    }
    
    case roRMS  : 
    {
      const float ScalingCoeff = 1.0f / (TParameters::GetInstance()->Get_Nt() - TParameters::GetInstance()->GetStartTimeIndex());
      
      #pragma omp parallel for if (BufferSize > 1e6)
      for (size_t i = 0; i < BufferSize; i++)
      {
        StoreBuffer[i] = sqrt(StoreBuffer[i] * ScalingCoeff);
      }   
      break;
    }
    
    case roMAX  :
    {
      // do nothing
      break;
    }
    
    case roMIN  : 
    {
      // do nothing
      break;
    }    
  }// switch  
    
}// end of ApplyPostProcessing
//------------------------------------------------------------------------------




//----------------------------------------------------------------------------//
//                 TIndexOutputHDF5Stream implementation                      //
//                              public methods                                //
//----------------------------------------------------------------------------//

/**
 * Constructor - links the HDF5 dataset, source (sampled matrix), Sensor mask 
 * and the reduction operator together. The constructor DOES NOT allocate memory 
 * because the size of the sensor mask is not known at the time the instance of 
 * the class is being created.
 * 
 * @param [in] HDF5_File       - Handle to the HDF5 (output) file 
 * @param [in] HDF5_ObjectName - The dataset's name (index based sensor data 
 *                                is store in a single dataset)
 * @param [in] SourceMatrix    - The source matrix (only real matrices are 
 *                               supported)
 * @param [in] SensorMask      - Index based sensor mask
 * @param [in] ReductionOp     - Reduction operator
 * @param [in] BufferToReuse   - An external buffer can be used to line up 
 *                               the grid points
 */
TIndexOutputHDF5Stream::TIndexOutputHDF5Stream(THDF5_File &             HDF5_File,
                                               const char *             HDF5_ObjectName,
                                               const TRealMatrix &      SourceMatrix, 
                                               const TLongMatrix &      SensorMask,
                                               const TReductionOperator ReductionOp,
                                               float *                  BufferToReuse)
        : TBaseOutputHDF5Stream(HDF5_File, HDF5_ObjectName, SourceMatrix, ReductionOp, BufferToReuse),
          SensorMask(SensorMask),
          HDF5_DatasetId(H5I_BADID),          
          Position(0,0,0)
{
  
}// end of TIndexOutputHDF5Stream
//------------------------------------------------------------------------------


/**
 * Destructor
 * if the file is still opened, it applies the post processing and flush the data.
 * Then, the object memory is freed and the object destroyed
 */
TIndexOutputHDF5Stream::~TIndexOutputHDF5Stream()
{
  
  Close();
  // free memory only if it was allocated
  if (!BufferReuse) FreeMemory();        
}// end of Destructor
//------------------------------------------------------------------------------



/**
 * Create a HDF5 stream, create a dataset, and allocate data for it
 * @param [in] NumberOfSampledElementsPerStep
 */
void TIndexOutputHDF5Stream::Create()
{
 
  size_t NumberOfSampledElementsPerStep = SensorMask.GetTotalElementCount();
  
  TParameters * Params = TParameters::GetInstance();
      
  // Derive dataset dimension sizes
  TDimensionSizes DatasetSize(NumberOfSampledElementsPerStep, 
                              (ReductionOp == roNONE) ?  Params->Get_Nt() - Params->GetStartTimeIndex() : 1,
                              1);
  
  // Set HDF5 chunk size
  TDimensionSizes ChunkSize(NumberOfSampledElementsPerStep, 1, 1);  
  // for chunks bigger than 32 MB 
  if (NumberOfSampledElementsPerStep > (ChunkSize_4MB * 8))
  { 
    ChunkSize.X = ChunkSize_4MB; // set chunk size to MB
  }
  
  // Create a dataset under the root group
  HDF5_DatasetId = HDF5_File.CreateFloatDataset(HDF5_File.GetRootGroup(),
                                                HDF5_RootObjectName,
                                                DatasetSize,
                                                ChunkSize, 
                                                Params->GetCompressionLevel());
  
  // Write dataset parameters
  HDF5_File.WriteMatrixDomainType(HDF5_File.GetRootGroup(), 
                                  HDF5_RootObjectName,
                                  THDF5_File::hdf5_mdt_real);
  HDF5_File.WriteMatrixDataType  (HDF5_File.GetRootGroup(), 
                                  HDF5_RootObjectName,
                                  THDF5_File::hdf5_mdt_float);        
  
    
  // Set positions in the dataset
  Position = TDimensionSizes(0,0,0);
      
  // Set buffer size
  BufferSize = NumberOfSampledElementsPerStep;
  
  // Allocate memory if needed
  if (!BufferReuse) AllocateMemory();
    
}// end of Create
//------------------------------------------------------------------------------

/**
 * Sample grid points, line them up in the buffer an flush to the disk unless a
 * reduction operator is applied
 */
void TIndexOutputHDF5Stream::Sample()
{  
  //@TODO: Do not use SensorMask as an object but as an array
  switch (ReductionOp)
  {
    case roNONE :
    {
      #pragma omp parallel for if (BufferSize > 1e6)  
      for (size_t i = 0; i < BufferSize; i++)
      {
        StoreBuffer[i] = SourceMatrix[SensorMask[i]];
      }    
      // only raw time series are flushed down to the disk every time step
      FlushToFile();
      break;
    }
    
    case roRMS  : 
    {
      #pragma omp parallel for if (BufferSize > 1e6)  
      for (size_t i = 0; i < BufferSize; i++)
      {
        StoreBuffer[i] += (SourceMatrix[SensorMask[i]] * SourceMatrix[SensorMask[i]]);
      }
      break;
    }
    
    case roMAX  :
    {
      #pragma omp parallel for if (BufferSize > 1e6)    
      for (size_t i = 0; i < BufferSize; i++)
      {
        if (StoreBuffer[i] < SourceMatrix[SensorMask[i]])
          StoreBuffer[i] = SourceMatrix[SensorMask[i]];
      }
      break;
    }
    
    case roMIN  : 
    {
      #pragma omp parallel for if (BufferSize > 1e6)  
      for (size_t i = 0; i < BufferSize; i++)
      {
        if (StoreBuffer[i] > SourceMatrix[SensorMask[i]])
          StoreBuffer[i] = SourceMatrix[SensorMask[i]];
      }
      break;
    }    
  }// switch  
}// end of Sample
//------------------------------------------------------------------------------

/**
 * Close stream (apply post-processing if necessary, flush data and close)
 */
void TIndexOutputHDF5Stream::Close()
{
  // the dataset is still opened
  if (HDF5_DatasetId != H5I_BADID) 
  {
    ApplyPostProcessing();
    // When no reduction operator is applied, the data is flushed after every time step
    if (ReductionOp != roNONE) FlushToFile();
    
    HDF5_File.CloseDataset(HDF5_DatasetId);
  }
    
  HDF5_DatasetId = H5I_BADID;
}// end of Close
//------------------------------------------------------------------------------

//----------------------------------------------------------------------------//
//                 TIndexOutputHDF5Stream implementation                      //
//                            protected methods                               //
//----------------------------------------------------------------------------//


/**
 * Flush the buffer down to the file at the actual position
 */
void TIndexOutputHDF5Stream::FlushToFile()
{
  HDF5_File.WriteHyperSlab(HDF5_DatasetId, 
                           Position, 
                           TDimensionSizes(BufferSize,1,1), 
                           StoreBuffer);
  Position.Y++;
}// end of FlushToFile
//------------------------------------------------------------------------------






//----------------------------------------------------------------------------//
//                TCuboidOutputHDF5Stream implementation                      //
//                              public methods                                //
//----------------------------------------------------------------------------//


/**
 * Constructor - links the HDF5 dataset, SourceMatrix, and SensorMask together
 * @param HDF5_File       - HDF5 file to write the output to 
 * @param HDF5_GroupName  - The name of the HDF5 group. This group contains datasets for particular cuboids
 * @param SourceMatrix    - Source matrix to be sampled      
 * @param SensorMask      - Sensor mask with the cuboid coordinates
 * @param ReductionOp     - Reduction operator
 * @param BufferToReuse   - If there is a memory space to be reused, provide a pointer
 */
TCuboidOutputHDF5Stream::TCuboidOutputHDF5Stream(THDF5_File &             HDF5_File,
                                                 const char *             HDF5_GroupName,
                                                 const TRealMatrix &      SourceMatrix,
                                                 const TLongMatrix &      SensorMask,
                                                 const TReductionOperator ReductionOp,
                                                 float *                  BufferToReuse)
        : TBaseOutputHDF5Stream(HDF5_File, HDF5_GroupName, SourceMatrix, ReductionOp, BufferToReuse),
          SensorMask(SensorMask),
          HDF5_GroupId(H5I_BADID),
          SampledTimeStep(0)
{
  
}// end of TCubodidOutputHDF5Stream
//------------------------------------------------------------------------------


/**
 * Destructor
 * if the file is still opened, it applies the post processing and flush the data.
 * Then, the object memory is freed and the object destroyed
 */
TCuboidOutputHDF5Stream::~TCuboidOutputHDF5Stream()
{
  Close();
  
  // free memory only if it was allocated
  if (!BufferReuse) FreeMemory();        
    
}// end ~TCubodidOutputHDF5Stream
//------------------------------------------------------------------------------




/**
 * Create a HDF5 stream and allocate data for it. It also creates a HDF5 group
 * with particular datasets (one per cuboid) 
 */
void TCuboidOutputHDF5Stream::Create()
{     
  // Create the HDF5 group and open it    
  HDF5_GroupId = HDF5_File.CreateGroup(HDF5_File.GetRootGroup(), HDF5_RootObjectName);
                      
  
  // Create all datasets (sizes, chunks, and attributes)
  size_t NumberOfCuboids        = SensorMask.GetDimensionSizes().Y;
  CuboidsInfo.reserve(NumberOfCuboids);  
  size_t ActualPositionInBuffer = 0;
  
  for (size_t CuboidIndex = 0; CuboidIndex < NumberOfCuboids; CuboidIndex++)
  {
    TCuboidInfo CuboidInfo;
    
    CuboidInfo.HDF5_CuboidId = CreateCuboidDataset(CuboidIndex);
    CuboidInfo.StartingPossitionInBuffer = ActualPositionInBuffer;
    
    ActualPositionInBuffer += (SensorMask.GetBottomRightCorner(CuboidIndex) - SensorMask.GetTopLeftCorner(CuboidIndex)).GetElementCount();
    CuboidsInfo.push_back(CuboidInfo);
  }
  
  //we're at the beginning
  SampledTimeStep = 0;
  
  // Create the memory buffer if necessary and set starting address  
  BufferSize = SensorMask.GetTotalNumberOfElementsInAllCuboids();
  
  // Allocate memory if needed
  if (!BufferReuse) AllocateMemory();
    
}// end of Create
//------------------------------------------------------------------------------



/**
 * Sample data into buffer and apply reduction, or flush to disk - based on a sensor mask
 */
void TCuboidOutputHDF5Stream::Sample()
{
  
  const size_t XY_Size = SourceMatrix.GetDimensionSizes().Y * SourceMatrix.GetDimensionSizes().X;
  const size_t X_Size  = SourceMatrix.GetDimensionSizes().X;

  switch (ReductionOp)
  {
    case roNONE :
    {
      TDimensionSizes DatasetPosition(0,0,0,0); //4D position in the dataset
      TDimensionSizes CuboidSize(0,0,0,0);       // Size of the cuboid

      DatasetPosition.T = SampledTimeStep;
      const float * MatrixData = SourceMatrix.GetRawData();

      // iterate over all cuboid to be sampled
      for (size_t CuboidIndex = 0; CuboidIndex < CuboidsInfo.size(); CuboidIndex++)
      {
        CuboidSize = SensorMask.GetBottomRightCorner(CuboidIndex) - SensorMask.GetTopLeftCorner(CuboidIndex);
        CuboidSize.T = 1; 

        HDF5_File.WriteCuboidToHyperSlab(CuboidsInfo[CuboidIndex].HDF5_CuboidId,
                                         DatasetPosition,
                                         SensorMask.GetTopLeftCorner(CuboidIndex), // position in the SourceMatrix
                                         CuboidSize,
                                         SourceMatrix.GetDimensionSizes(),
                                         MatrixData
        );
      }
      SampledTimeStep++;   // Move forward in time

      break;
    }
    
    case roRMS  : 
    {
      if (ReductionOp == roRMS)
      {    
        size_t BufferIndex = 0;

        // @TODO - figure out how to do this with OpenMP
        for (size_t CuboidIdx = 0; CuboidIdx < SensorMask.GetDimensionSizes().Y; CuboidIdx++)
        {
          const TDimensionSizes TopLeftCorner     = SensorMask.GetTopLeftCorner(CuboidIdx);
          const TDimensionSizes BottomRightCorner = SensorMask.GetBottomRightCorner(CuboidIdx);

          for (size_t z = TopLeftCorner.Z; z <= BottomRightCorner.Z; z++)
            for (size_t y = TopLeftCorner.Y; y <= BottomRightCorner.Y; y++)
              for (size_t x = TopLeftCorner.X; x <= BottomRightCorner.X; x++)
              {
                StoreBuffer[BufferIndex++] += (SourceMatrix[z * XY_Size + y * X_Size + x] * 
                                                 SourceMatrix[z * XY_Size + y * X_Size + x]);                         
              }
        }
      }
  
      break;
    }
    
    case roMAX  :
    {
      size_t BufferIndex = 0;
      ///@TODO Parallel section with some load balancing feature
      for (size_t CuboidIdx = 0; CuboidIdx < SensorMask.GetDimensionSizes().Y; CuboidIdx++)
      {
        const TDimensionSizes TopLeftCorner = SensorMask.GetTopLeftCorner(CuboidIdx);
        const TDimensionSizes BottomRightCorner = SensorMask.GetBottomRightCorner(CuboidIdx);

        for (size_t z = TopLeftCorner.Z; z <= BottomRightCorner.Z; z++)
          for (size_t y = TopLeftCorner.Y; y <= BottomRightCorner.Y; y++)
            for (size_t x = TopLeftCorner.X; x <= BottomRightCorner.X; x++)
            {
              if (StoreBuffer[BufferIndex] < SourceMatrix[z * XY_Size + y * X_Size + x])
              {
                StoreBuffer[BufferIndex] = SourceMatrix[z * XY_Size + y * X_Size + x];             
              }
              BufferIndex++;              
            }
      }
      break;
    }
    
    case roMIN  : 
    {
      size_t BufferIndex = 0;
      ///@TODO Parallel section with some load balancing feature
      for (size_t CuboidIdx = 0; CuboidIdx < SensorMask.GetDimensionSizes().Y; CuboidIdx++)
      {
        const TDimensionSizes TopLeftCorner = SensorMask.GetTopLeftCorner(CuboidIdx);
        const TDimensionSizes BottomRightCorner = SensorMask.GetBottomRightCorner(CuboidIdx);

        for (size_t z = TopLeftCorner.Z; z <= BottomRightCorner.Z; z++)
          for (size_t y = TopLeftCorner.Y; y <= BottomRightCorner.Y; y++)
            for (size_t x = TopLeftCorner.X; x <= BottomRightCorner.X; x++)
            {
              if (StoreBuffer[BufferIndex] > SourceMatrix[z * XY_Size + y * X_Size + x])
              {
                StoreBuffer[BufferIndex] = SourceMatrix[z * XY_Size + y * X_Size + x];             
              }
              BufferIndex++;              
            }
      }
       
      break;
    }    
  }// switch  
   
}// end of Sample
//------------------------------------------------------------------------------


/**
 * Close stream (apply post-processing if necessary, flush data, close datasets 
 * and the group )
 */
void TCuboidOutputHDF5Stream::Close()
{
  // the group is still open
  if (HDF5_GroupId != H5I_BADID) 
  {    
    ApplyPostProcessing();
    
    // When no reduction operator is applied, the data is flushed after every time step
    if (ReductionOp != roNONE) FlushBufferToFile();
    
    // Close all datasets and the group    
    for (size_t CuboidIndex = 0; CuboidIndex < CuboidsInfo.size(); CuboidIndex++)
    {
      HDF5_File.CloseDataset(CuboidsInfo[CuboidIndex].HDF5_CuboidId);
    }        
    CuboidsInfo.clear();    
    
    HDF5_File.CloseGroup(HDF5_GroupId);
    HDF5_GroupId = H5I_BADID;  
  }// if opened    
}// end of Close
//------------------------------------------------------------------------------



//----------------------------------------------------------------------------//
//                 TIndexOutputHDF5Stream implementation                      //
//                            protected methods                               //
//----------------------------------------------------------------------------//

/**
 *  Create a new dataset for a given cuboid specified by index (order)
 * @param Index - Index of the cuboid in the sensor mask
 * @return HDF5 handle to the dataset
 */
hid_t TCuboidOutputHDF5Stream::CreateCuboidDataset(const size_t Index)
{
  TParameters * Params = TParameters::GetInstance();
  
  // if time series then Number of steps else 1  
  size_t NumberOfSampledTimeSteps = (ReductionOp == roNONE)
                                      ? Params->Get_Nt() - Params->GetStartTimeIndex() 
                                      : 0; // will be a 3D dataset
  // Set cuboid dimensions   
  TDimensionSizes CuboidSize((SensorMask.GetBottomRightCorner(Index) - SensorMask.GetTopLeftCorner(Index)).X,
                             (SensorMask.GetBottomRightCorner(Index) - SensorMask.GetTopLeftCorner(Index)).Y,
                             (SensorMask.GetBottomRightCorner(Index) - SensorMask.GetTopLeftCorner(Index)).Z,
                             NumberOfSampledTimeSteps
                            );
  
  // Set chunk size   
  // If the size of the cuboid is bigger than 32 MB, set the chunk to approx 4MB
  size_t NumberOfSlabs = 0;
  TDimensionSizes CuboidChunkSize(CuboidSize.X, CuboidSize.Y, CuboidSize.Z, (ReductionOp == roNONE) ? 1 : 0);
  
  if (CuboidSize.GetElementCount() > (ChunkSize_4MB * 8))
  {
    while (NumberOfSlabs * CuboidSize.X * CuboidSize.Y > ChunkSize_4MB) NumberOfSlabs++;
    CuboidChunkSize.Z = NumberOfSlabs;
  }
    
  // @TODO: Can be done easily with std::to_string and c++0x or c++-11
  char HDF5_DatasetName[32] = "";
  // Indexed from 1
  sprintf(HDF5_DatasetName, "%ld",Index+1);
  hid_t HDF5_DatasetId = HDF5_File.CreateFloatDataset(HDF5_GroupId,
                                                      HDF5_DatasetName,
                                                      CuboidSize,
                                                      CuboidChunkSize, 
                                                      Params->GetCompressionLevel()
                                                     );
  
  // Write dataset parameters
  HDF5_File.WriteMatrixDomainType(HDF5_GroupId, 
                                  HDF5_DatasetName,
                                  THDF5_File::hdf5_mdt_real);
  HDF5_File.WriteMatrixDataType  (HDF5_GroupId, 
                                  HDF5_DatasetName,
                                  THDF5_File::hdf5_mdt_float);        
  
  
  return HDF5_DatasetId;
  
}//end of CreateCuboidDatasets
//------------------------------------------------------------------------------


    
/**
 * Flush the buffer to the file (to multiple datasets if necessary)
 */
void TCuboidOutputHDF5Stream::FlushBufferToFile()
{
  
  TDimensionSizes Position (0,0,0,0);
  TDimensionSizes BlockSize(0,0,0,0);
  
  if (ReductionOp == roNONE) Position.T = SampledTimeStep;
  
  for (size_t CuboidIndex = 0; CuboidIndex < CuboidsInfo.size(); CuboidIndex++)
  {       
    
    BlockSize = SensorMask.GetBottomRightCorner(CuboidIndex) - SensorMask.GetTopLeftCorner(CuboidIndex);
    BlockSize.T = 1;
    
    HDF5_File.WriteHyperSlab(CuboidsInfo[CuboidIndex].HDF5_CuboidId,
                             Position,
                             BlockSize,
                             StoreBuffer + CuboidsInfo[CuboidIndex].StartingPossitionInBuffer
                             );

  }
  
  SampledTimeStep++;    
}// end of FlushBufferToFile
//------------------------------------------------------------------------------



//----------------------------------------------------------------------------//
//             TWholeDomainOutputHDF5Stream implementation                    //
//                              public methods                                //
//----------------------------------------------------------------------------//

/**
 * Constructor - links the HDF5 dataset and SourceMatrix
 * @param HDF5_File        - HDF5 file to write the output to 
 * @param HDF5_DatasetName - The name of the HDF5 group. This group contains datasets for particular cuboids
 * @param SourceMatrix     - Source matrix to be sampled      
 * @param ReductionOp      - Reduction operator
 * @param BufferToReuse    - If there is a memory space to be reused, provide a pointer
 */
TWholeDomainOutputHDF5Stream::TWholeDomainOutputHDF5Stream(THDF5_File &             HDF5_File,
                                                           const char *             HDF5_DatasetName,
                                                           const TRealMatrix &      SourceMatrix,
                                                           const TReductionOperator ReductionOp,
                                                           float *                  BufferToReuse)
        : TBaseOutputHDF5Stream(HDF5_File, HDF5_DatasetName, SourceMatrix, ReductionOp, BufferToReuse),          
          HDF5_DatasetId(H5I_BADID),
          SampledTimeStep(0)
{
  
}// end of TWholeDomainOutputHDF5Stream
//------------------------------------------------------------------------------


/**
 * Destructor
 * if the file is still opened, it applies the post processing and flush the data.
 * Then, the object memory is freed and the object destroyed
 */
TWholeDomainOutputHDF5Stream::~TWholeDomainOutputHDF5Stream()
{  
  Close();
  // free memory only if it was allocated
  if (!BufferReuse) FreeMemory();        
}// end of Destructor
//------------------------------------------------------------------------------



/**
 * Create a HDF5 stream for the whole domain and allocate data for it
 */
void TWholeDomainOutputHDF5Stream::Create()
{  
  
  TDimensionSizes ChunkSize(SourceMatrix.GetDimensionSizes().X, SourceMatrix.GetDimensionSizes().Y, 1);
  
  // Create a dataset under the root group
  HDF5_DatasetId = HDF5_File.CreateFloatDataset(HDF5_File.GetRootGroup(),
                                                HDF5_RootObjectName,
                                                SourceMatrix.GetDimensionSizes(),
                                                ChunkSize, 
                                                TParameters::GetInstance()->GetCompressionLevel());
  
  // Write dataset parameters
  HDF5_File.WriteMatrixDomainType(HDF5_File.GetRootGroup(), 
                                  HDF5_RootObjectName,
                                  THDF5_File::hdf5_mdt_real);
  HDF5_File.WriteMatrixDataType  (HDF5_File.GetRootGroup(), 
                                  HDF5_RootObjectName,
                                  THDF5_File::hdf5_mdt_float);        
            
  // Set buffer size
  BufferSize = SourceMatrix.GetTotalElementCount();
  
  // Allocate memory if needed
  if (!BufferReuse) AllocateMemory();
      
}//end of Create
//------------------------------------------------------------------------------


/**
 * Sample all grid points, line them up in the buffer an flush to the disk unless 
 * a reduction operator is applied
 */
void TWholeDomainOutputHDF5Stream::Sample()
{  
  switch (ReductionOp)
  {
    case roNONE :
    { //@TODO: Use HDF5 with File space and Memory space
      #pragma omp parallel for if (BufferSize > 1e6)  
      for (size_t i = 0; i < BufferSize; i++)
      {
        StoreBuffer[i] = SourceMatrix[i];
      }    
      // only raw time series are flushed down to the disk every time step
      FlushToFile();
      break;
    }
    
    case roRMS  : 
    {
      #pragma omp parallel for if (BufferSize > 1e6)  
      for (size_t i = 0; i < BufferSize; i++)
      {
        StoreBuffer[i] += (SourceMatrix[i] * SourceMatrix[i]);
      }
      break;
    }
    
    case roMAX  :
    {
      #pragma omp parallel for if (BufferSize > 1e6)    
      for (size_t i = 0; i < BufferSize; i++)
      {
        if (StoreBuffer[i] < SourceMatrix[i])
          StoreBuffer[i] = SourceMatrix[i];
      }
      break;
    }
    
    case roMIN  : 
    {
      #pragma omp parallel for if (BufferSize > 1e6)  
      for (size_t i = 0; i < BufferSize; i++)
      {
        if (StoreBuffer[i] > SourceMatrix[i])
          StoreBuffer[i] = SourceMatrix[i];
      }
      break;
    }    
  }// switch  
}// end of Sample
//------------------------------------------------------------------------------

/**
 * Close stream (apply post-processing if necessary, flush data and close)
 */
void TWholeDomainOutputHDF5Stream::Close()
{
  // the dataset is still opened
  if (HDF5_DatasetId != H5I_BADID) 
  {
    ApplyPostProcessing();
    // When no reduction operator is applied, the data is flushed after every time step
    if (ReductionOp != roNONE) FlushToFile();
    
    HDF5_File.CloseDataset(HDF5_DatasetId);
  }
    
  HDF5_DatasetId = H5I_BADID;
}// end of Close
//------------------------------------------------------------------------------


//----------------------------------------------------------------------------//
//                TWholeDomainOutputHDF5Stream implementation                 //
//                            protected methods                               //
//----------------------------------------------------------------------------//


/**
 * Flush the buffer down to the file at the actual position
 */
void TWholeDomainOutputHDF5Stream::FlushToFile()
{
  TDimensionSizes Size = SourceMatrix.GetDimensionSizes();
  TDimensionSizes Position(0,0,0);    
  
  
  if (ReductionOp == roNONE) 
  {
    Position.T = SampledTimeStep;
    Size.T = SampledTimeStep;
  }
          
  HDF5_File.WriteHyperSlab(HDF5_DatasetId, 
                           Position, 
                           Size, 
                           StoreBuffer);
  SampledTimeStep++;
}// end of FlushToFile
//------------------------------------------------------------------------------


