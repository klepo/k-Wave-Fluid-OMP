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
 *              21 March    2014, 15:20      (revised)
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
  StoringBuffer = (float *) memalign(DATA_ALIGNMENT, BufferSize * sizeof (float));

  if (!StoringBuffer)
  {
    fprintf(stderr, Matrix_ERR_FMT_NotEnoughMemory, "TBaseOutputHDF5Stream");
    throw bad_alloc();
  }

  // zero the matrix
  #pragma omp parallel for if (BufferSize > 1e6)  
  for (size_t i = 0; i < BufferSize; i++)
  {
    StoringBuffer[i] = 0.0f;
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
  if (StoringBuffer)
  {
    free(StoringBuffer);
    StoringBuffer = NULL;
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
    roNONE :
    {
      // do nothing
      break;
    }
    
    roRMS  : 
    {
      const float ScalingCoeff = TParameters::GetInstance()->Get_Nt() - TParameters::GetInstance()->GetStartTimeIndex();
      
      #pragma omp parallel for if (BufferSize > 1e6)
      for (size_t i = 0; i < BufferSize; i++)
      {
        StoringBuffer[i] = sqrt(StoringBuffer[i] * ScalingCoeff);
      }   
      break;
    }
    
    roMAX  :
    {
      // do nothing
      break;
    }
    
    roMIN  : 
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
 * @param [in] HDF5_File           - Handle to the HDF5 (output) file 
 * @param [in] HDF5_RootObjectName - The dataset's name (index based sensor data 
 *                                   is store in a single dataset)
 * @param [in] SourceMatrix        - The source matrix (only real matrices are 
 *                                   supported)
 * @param [in] SensorMask          - Index based sensor mask
 * @param [in] ReductionOp         - Reduction operator
 * @param [in] BufferToReuse       - An external buffer can be used to line up 
 *                                   the grid points
 */
TIndexOutputHDF5Stream::TIndexOutputHDF5Stream(THDF5_File &             HDF5_File,
                                               const char *             HDF5_RootObjectName,
                                               const TRealMatrix &      SourceMatrix, 
                                               const TLongMatrix        SensorMask,
                                               const TReductionOperator ReductionOp,
                                               float *                  BufferToReuse)
        : TBaseOutputHDF5Stream(HDF5_File, HDF5_RootObjectName, SourceMatrix, ReductionOp, BufferToReuse),
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
void TIndexOutputHDF5Stream::Create(const size_t NumberOfSampledElementsPerStep)
{
 
  TParameters * Params = TParameters::GetInstance();
      
  // Derive dataset dimension sizes
  TDimensionSizes DatasetSize(NumberOfSampledElementsPerStep, 
                              Params->GetNumberOfThreads() - Params->GetStartTimeIndex(),
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
  switch (ReductionOp)
  {
    roNONE :
    {
      #pragma omp parallel for if (BufferSize > 1e6)  
      for (size_t i = 0; i < BufferSize; i++)
      {
        StoringBuffer[i] = SourceMatrix[SensorMask[i]];
      }    
      // only raw time series are flushed down to the disk every time step
      FlushToFile();
      break;
    }
    
    roRMS  : 
    {
      #pragma omp parallel for if (BufferSize > 1e6)  
      for (size_t i = 0; i < BufferSize; i++)
      {
        StoringBuffer[i] += (SourceMatrix[SensorMask[i]] * SourceMatrix[SensorMask[i]]);
      }
      break;
    }
    
    roMAX  :
    {
      #pragma omp parallel for if (BufferSize > 1e6)    
      for (size_t i = 0; i < BufferSize; i++)
      {
        if (StoringBuffer[i] < SourceMatrix[SensorMask[i]])
          StoringBuffer[i] = SourceMatrix[SensorMask[i]];
      }
      break;
    }
    
    roMIN  : 
    {
      #pragma omp parallel for if (BufferSize > 1e6)  
      for (size_t i = 0; i < BufferSize; i++)
      {
        if (StoringBuffer[i] > SourceMatrix[SensorMask[i]])
          StoringBuffer[i] = SourceMatrix[SensorMask[i]];
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
                           StoringBuffer);
  Position.Y++;
}// end of FlushToFile
//------------------------------------------------------------------------------






//----------------------------------------------------------------------------//
//                 TIndexOutputHDF5Stream implementation                      //
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




//----------------------------------------------------------------------------//
//                              Implementation                                //
//                              public methods                                //
//----------------------------------------------------------------------------//
/**
 * This method creates the time series output stream by creating either a HDF5 
 * group of datasets dataset or a dataset and writing the parameters of the dataset
 * 
 * @param [in,out] HDF5_File                - Handle to HDF5 output file
 * @param [in]     DatasetName              - Dataset name
 * @param [in]     NumberOfSensorElements   - Number of elements sampled 
 * @param [in]     NumberOfTimeSteps        - Number of time steps to sample
 * @param [in,out] BufferToReuse            - if NULL the class will create its own memory
 *                                            if !NULL the class will reuse specified memory 
 *                                            to linearise the sampled data
 */
void TOutputHDF5Stream::CreateStream(THDF5_File & HDF5_File,
                                     const char * DatasetName,
                                     const size_t NumberOfSensorElements,
                                     const size_t NumberOfTimeSteps,
                                     float *      BufferToReuse)
{
  
  WholeDomainStream = false;
  
  // if reuse buffer, set the pointer
  BufferReuse = (BufferToReuse != NULL);
  StoringBuffer = BufferToReuse;
          
  // copy the dataset name
  this->DatasetName = new char[strlen(DatasetName)];
  strcpy(this->DatasetName, DatasetName);

  // Set the HDF5 file
  this->HDF5_File = &HDF5_File;
    
  // set HDF5 domain size
  TotalDatasetSize.X = NumberOfSensorElements;
  TotalDatasetSize.Y = NumberOfTimeSteps;
  TotalDatasetSize.Z = 1;
  
  // Set HDF5 chunk size
  TDimensionSizes ChunkSize(NumberOfSensorElements, 1, 1);  
  // for chunks bigger than 32 MB 
  if (NumberOfSensorElements > (ChunkSize_4MB * 8))
  { 
    ChunkSize.X = ChunkSize_4MB; // set chunk size to MB
  }
  
  // Set positions in the dataset
  Position = TDimensionSizes(0,0,0);
    
  hid_t GroupId =  HDF5_File.CreateGroup(DatasetName);
  // Create float dataset in the file
  HDF5_Dataset_id = HDF5_File.CreateFloatDataset(GroupId, "Test",TotalDatasetSize,
                                                 ChunkSize, TParameters::GetInstance()->GetCompressionLevel());
  // Write dataset parameters
  HDF5_File.WriteMatrixDomainType(DatasetName,THDF5_File::hdf5_mdt_real);
  HDF5_File.WriteMatrixDataType(DatasetName,THDF5_File::hdf5_mdt_float);        
        
  // Set buffer size
  BufferSize = TDimensionSizes(NumberOfSensorElements,1,1);
  
  // Allocate memory if needed
  if (!BufferReuse) AllocateMemory();
  
    
}// end CreateStream
//------------------------------------------------------------------------------

/**
 * Create stream operating on the whole domain (no sensor_mask), only for aggregated quantities
 * @param [in,out] HDF5_File
 * @param [in] DatasetName
 */
void TOutputHDF5Stream::CreateWholeDomainAggregatedStream(THDF5_File & HDF5_File,
                                                          const char * DatasetName)
{

  WholeDomainStream = true;
  // Never reuse buffer for aggregated quantities
  BufferReuse = false;

  // copy the dataset name
  this->DatasetName = new char[strlen(DatasetName)];
  strcpy(this->DatasetName, DatasetName);

  // Set the HDF5 file
  this->HDF5_File = &HDF5_File;

  // set HDF5 domain size
  TotalDatasetSize = TParameters::GetInstance()->GetFullDimensionSizes();

  // Set HDF5 chunk size to be a 1D slab
  TDimensionSizes ChunkSize(TotalDatasetSize.X, TotalDatasetSize.Y, 1);

  // Set positions in the dataset
  Position = TDimensionSizes(0, 0, 0);


  // Create float dataset in the file  
  HDF5_Dataset_id = HDF5_File.CreateFloatDataset(DatasetName, TotalDatasetSize, 
                                                 ChunkSize, TParameters::GetInstance()->GetCompressionLevel());

  // Write dataset parameters
  HDF5_File.WriteMatrixDomainType(DatasetName, THDF5_File::hdf5_mdt_real);
  HDF5_File.WriteMatrixDataType(DatasetName, THDF5_File::hdf5_mdt_float);

  // Set buffer size
  BufferSize = TotalDatasetSize;
  
  // Allocate memory - always
  AllocateMemory();


}// end of CreateWholeDomainAggregatedStream
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
 * Write data directly into the dataset based on the sensor mask (linear or cuboid)
 * @param [in] SourceMatrix  - Matrix to be sampled and stored
 * @param [in] SensorMask - Sensor mask 
 */
void TOutputHDF5Stream::WriteDataDirectly( const TRealMatrix& SourceMatrix,
                                           const TLongMatrix& SensorMask)
{

  // Sample points based on linear mask index
  if (TParameters::GetInstance()->Get_sensor_mask_type() == TParameters::smt_index)
  {
    #pragma omp parallel for if (SensorMask.GetTotalElementCount() > 1e6)  
    for (size_t i = 0; i < SensorMask.GetTotalElementCount(); i++)
    {
      StoringBuffer[i] = SourceMatrix[SensorMask[i]];
    }    
  }
  
  // Sample points based on cuboid corners sensor mask
  if (TParameters::GetInstance()->Get_sensor_mask_type() == TParameters::smt_corners)
  {
    
    size_t BufferIndex = 0;
    const size_t XY_Size = SourceMatrix.GetDimensionSizes().Y * SourceMatrix.GetDimensionSizes().X;
    const size_t X_Size = SourceMatrix.GetDimensionSizes().X;

    ///@TODO Parallel section with some load balancing feature
    for (size_t CuboidIdx = 0; CuboidIdx < SensorMask.GetDimensionSizes().Y; CuboidIdx++)
    {
      const TDimensionSizes TopLeftCorner     = SensorMask.GetTopLeftCorner(CuboidIdx);
      const TDimensionSizes BottomRightCorner = SensorMask.GetBottomRightCorner(CuboidIdx);

      for (size_t z = TopLeftCorner.Z; z <= BottomRightCorner.Z; z++)
        for (size_t y = TopLeftCorner.Y; y <= BottomRightCorner.Y; y++)
          for (size_t x = TopLeftCorner.X; x <= BottomRightCorner.X; x++)
          {
            StoringBuffer[BufferIndex++] = SourceMatrix[z * XY_Size + y * X_Size + x];
          }
    }

  } 
      
  // Write storing buffer into the dataset
  WriteBuffer(StoringBuffer);
  
}// end of WriteDataDirectly
//------------------------------------------------------------------------------    


/**
 * Reduce data into buffer - add data from matrix,
 * @warning Distribution given by indices must be compatible 
 *      with dimension sizes provide when the stream was opened!
 * @param [in] SourceMatrix - Matrix to be sampled
 * @param [in] SensorMask   - Sensor mask (corner or cuboid)
 * @param [in] ReductionOp  - Reduction operator to be applied
 */
void TOutputHDF5Stream::ReduceDataIntoBuffer(const TRealMatrix& SourceMatrix,
                                             const TLongMatrix& SensorMask,
                                             const TReductionOperator ReductionOp)
{
  // Sample points based on linear mask index
  if (TParameters::GetInstance()->Get_sensor_mask_type() == TParameters::smt_index)
  {
    if (ReductionOp == roRMS)
    {
      #pragma omp parallel for if (SensorMask.GetTotalElementCount() > 1e6)  
      for (size_t i = 0; i < BufferSize.GetElementCount(); i++)
      {
        StoringBuffer[i] += (SourceMatrix[SensorMask[i]] * SourceMatrix[SensorMask[i]]);
      }
    }
    else if (ReductionOp == roMAX)
    {
      #pragma omp parallel for if (SensorMask.GetTotalElementCount() > 1e6)    
      for (size_t i = 0; i < BufferSize.GetElementCount(); i++)
      {
        if (StoringBuffer[i] < SourceMatrix[SensorMask[i]])
          StoringBuffer[i] = SourceMatrix[SensorMask[i]];
      }
    }
    else if (ReductionOp == roMIN)
    {
      #pragma omp parallel for if (SensorMask.GetTotalElementCount() > 1e6)  
      for (size_t i = 0; i < BufferSize.GetElementCount(); i++)
      {
        if (StoringBuffer[i] > SourceMatrix[SensorMask[i]])
          StoringBuffer[i] = SourceMatrix[SensorMask[i]];
      }
    }    
  }// smt_index
    
  // Sample points based on cuboid corners sensor mask
  if (TParameters::GetInstance()->Get_sensor_mask_type() == TParameters::smt_corners)
  {    
    const size_t XY_Size = SourceMatrix.GetDimensionSizes().Y * SourceMatrix.GetDimensionSizes().X;
    const size_t X_Size = SourceMatrix.GetDimensionSizes().X;

    if (ReductionOp == roRMS)
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
              StoringBuffer[BufferIndex++] += (SourceMatrix[z * XY_Size + y * X_Size + x] * SourceMatrix[z * XY_Size + y * X_Size + x]);                         
            }
      }
    }
    
    if (ReductionOp == roMAX)
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
              if (StoringBuffer[BufferIndex] < SourceMatrix[z * XY_Size + y * X_Size + x]) StoringBuffer[BufferIndex] = SourceMatrix[z * XY_Size + y * X_Size + x];             
              BufferIndex++;              
            }
      }
    }
    
    if (ReductionOp == roMIN)
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
              if (StoringBuffer[BufferIndex] > SourceMatrix[z * XY_Size + y * X_Size + x]) StoringBuffer[BufferIndex] = SourceMatrix[z * XY_Size + y * X_Size + x];             
              BufferIndex++;              
            }
      }
    }    
  }   
}// end of ReduceDataIntoBuffer
//------------------------------------------------------------------------------


/**
 * Add data into reduce buffer (rms, max, min), all domain
 * Processes the full domain
 * @param [in] SourceMatrix
 * @param [in] ReductionOp
 */
void TOutputHDF5Stream::ReduceDataIntoBuffer(const TRealMatrix& SourceMatrix,
                                             const TReductionOperator ReductionOp)
{

  if (ReductionOp == roRMS)
  {
    #pragma omp parallel for if (BufferSize.GetElementCount() > 1e6)  
    for (size_t i = 0; i < BufferSize.GetElementCount(); i++)
    {
      StoringBuffer[i] += (SourceMatrix[i] * SourceMatrix[i]);
    }
  }
  
  if (ReductionOp == roMAX)
  {
    #pragma omp parallel for if (BufferSize.GetElementCount() > 1e6)  
    for (size_t i = 0; i < BufferSize.GetElementCount(); i++)
    {
      if (StoringBuffer[i] < SourceMatrix[i]) StoringBuffer[i] = SourceMatrix[i];

    }
  }
  
  if (ReductionOp == roMIN)
  {
    #pragma omp parallel for if (BufferSize.GetElementCount() > 1e6)  
    for (size_t i = 0; i < BufferSize.GetElementCount(); i++)
    {
      if (StoringBuffer[i] > SourceMatrix[i]) StoringBuffer[i] = SourceMatrix[i];

    }
  }
}// end of ReduceDataIntoBuffer
//------------------------------------------------------------------------------    

/**
 * RMS post-processing
 * @param [in] ScalingCoeff - scaling coefficient
 */
void TOutputHDF5Stream::RMSPostprocessing(const float ScalingCoeff)
{
  #pragma omp parallel for if (BufferSize.GetElementCount() > 1e6)
  for (size_t i = 0; i < BufferSize.GetElementCount(); i++)
  {
    StoringBuffer[i] = sqrt(StoringBuffer[i] * ScalingCoeff);
  }        
}// end of RMSPostprocessing
//------------------------------------------------------------------------------    

/**
 * Write reduced buffer (after post-processing)
 */
void TOutputHDF5Stream::WriteRudcedBuffer()
{
  WriteBuffer(StoringBuffer);
}// end of WriteRudcedBuffer
//------------------------------------------------------------------------------


/**
 * Destructor
 * 
 */
TOutputHDF5Stream::~TOutputHDF5Stream(){
    
    if (HDF5_Dataset_id != -1) HDF5_File->CloseDataset(HDF5_Dataset_id);
    // free memory only if it was allocated
    if (!BufferReuse) FreeMemory();    
    
}// end of Destructor
//------------------------------------------------------------------------------


//----------------------------------------------------------------------------//
//                              Implementation                                //
//                             protected methods                              //
//----------------------------------------------------------------------------//

/**
 * Allocate memory. 
 */
void TOutputHDF5Stream::AllocateMemory()
{
  StoringBuffer = (float *) memalign(DATA_ALIGNMENT, BufferSize.GetElementCount() * sizeof (float));

  if (!StoringBuffer)
  {
    fprintf(stderr, Matrix_ERR_FMT_NotEnoughMemory, "TOutputHDF5Stream");
    throw bad_alloc();
  }

  #pragma omp parallel for if (BufferSize.GetElementCount() > 1e6)  
  for (size_t i = 0; i < BufferSize.GetElementCount(); i++)
  {
    StoringBuffer[i] = 0.0f;
  }
}// end of AllocateMemory
//------------------------------------------------------------------------------

void TOutputHDF5Stream::FreeMemory()
{
  if (StoringBuffer)
  {
    free(StoringBuffer);
    StoringBuffer = NULL;
  }

  if (DatasetName)
  {
    delete [] DatasetName;
  }
}// end of FreeMemory
//------------------------------------------------------------------------------


/**
 * Write buffer into the dataset 
 * @param [in] Buffer - buffer to be written
 */
void TOutputHDF5Stream::WriteBuffer(const float * Buffer)
{
  HDF5_File->WriteHyperSlab(HDF5_Dataset_id, Position, BufferSize, StoringBuffer);
  Position.Y++;
}// end of WriteBuffer
//------------------------------------------------------------------------------


//----------------------------------------------------------------------------//
//                              Implementation                                //
//                              private methods                               //
//----------------------------------------------------------------------------//
