/**
 * @file        OutputHDF5Stream.cpp
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The implementation file of classes responsible for storing output
 *              quantities into the output HDF5 file
 *
 * @version     kspaceFirstOrder3D 2.16
 *
 * @date        11 July      2012, 10:30 (created) \n
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


#include <string.h>
#include <cmath>
#include <immintrin.h>
#include <MatrixClasses/OutputHDF5Stream.h>

#include <Parameters/Parameters.h>
#include <Utils/ErrorMessages.h>
#include <limits>


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


/**
 * Apply post-processing on the buffer. It supposes the elements are independent.
 *
 */
void TBaseOutputHDF5Stream::PostProcess()
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

      #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
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
//                  TBaseOutputHDF5Stream implementation                      //
//                              public methods                                //
//----------------------------------------------------------------------------//


/**
 * Allocate memory using a proper memory alignment.
 * @warning This can routine is not used in the base class (should be used in
 *          derived ones).
 */
void TBaseOutputHDF5Stream::AllocateMemory()
{
  StoreBuffer = (float *) _mm_malloc(BufferSize * sizeof (float), kDataAlignment);

  if (!StoreBuffer)
  {
    fprintf(stderr, kErrFmtNotEnoughMemory, "TBaseOutputHDF5Stream");
    throw bad_alloc();
  }

  // we need different initialization for different reduction ops
  switch (ReductionOp)
  {
    case roNONE :
    {
      // zero the matrix
      #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
      for (size_t i = 0; i < BufferSize; i++)
      {
        StoreBuffer[i] = 0.0f;
      }
      break;
    }

    case roRMS  :
    {
      // zero the matrix
      #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
      for (size_t i = 0; i < BufferSize; i++)
      {
        StoreBuffer[i] = 0.0f;
      }
      break;
    }

    case roMAX  :
    {
      // set the values to the highest negative float value
      #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
      for (size_t i = 0; i < BufferSize; i++)
      {
        StoreBuffer[i] = -1 * std::numeric_limits<float>::max();
      }
      break;
    }

    case roMIN  :
    {
      // set the values to the highest float value
      #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
      for (size_t i = 0; i < BufferSize; i++)
      {
        StoreBuffer[i] = std::numeric_limits<float>::max();
      }
      break;
    }
  }// switch
}// end of AllocateMemory
//------------------------------------------------------------------------------

/**
 * Free memory.
 * @warning This can routine is not used in the base class (should be used in
 *          derived ones).
 */
void TBaseOutputHDF5Stream::FreeMemory()
{
  if (StoreBuffer)
  {
    _mm_free(StoreBuffer);
    StoreBuffer = NULL;
  }
}// end of FreeMemory
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
TIndexOutputHDF5Stream::TIndexOutputHDF5Stream(Hdf5File &             HDF5_File,
                                               const char *             HDF5_ObjectName,
                                               const TRealMatrix &      SourceMatrix,
                                               const TIndexMatrix &     SensorMask,
                                               const TReductionOperator ReductionOp,
                                               float *                  BufferToReuse)
        : TBaseOutputHDF5Stream(HDF5_File, HDF5_ObjectName, SourceMatrix, ReductionOp, BufferToReuse),
          SensorMask(SensorMask),
          HDF5_DatasetId(H5I_BADID),
          SampledTimeStep(0)
{

}// end of TIndexOutputHDF5Stream
//------------------------------------------------------------------------------


/**
 * Destructor.
 * If the file is still opened, it applies the post processing and flush the data.
 * Then, the object memory is freed and the object destroyed.
 */
TIndexOutputHDF5Stream::~TIndexOutputHDF5Stream()
{
  Close();
  // free memory only if it was allocated
  if (!BufferReuse) FreeMemory();
}// end of Destructor
//------------------------------------------------------------------------------



/**
 * Create a HDF5 stream, create a dataset, and allocate data for it.
 */
void TIndexOutputHDF5Stream::Create()
{
  size_t NumberOfSampledElementsPerStep = SensorMask.GetTotalElementCount();

  TParameters * Params = TParameters::GetInstance();

  // Derive dataset dimension sizes
  DimensionSizes DatasetSize(NumberOfSampledElementsPerStep,
                              (ReductionOp == roNONE) ?  Params->Get_Nt() - Params->GetStartTimeIndex() : 1,
                              1);

  // Set HDF5 chunk size
  DimensionSizes ChunkSize(NumberOfSampledElementsPerStep, 1, 1);
  // for chunks bigger than 32 MB
  if (NumberOfSampledElementsPerStep > (ChunkSize_4MB * 8))
  {
    ChunkSize.nx = ChunkSize_4MB; // set chunk size to MB
  }

  // Create a dataset under the root group
  HDF5_DatasetId = HDF5_File.createDataset(HDF5_File.getRootGroup(),
                                           HDF5_RootObjectName,
                                           DatasetSize,
                                           ChunkSize,
                                           Hdf5File::MatrixDataType::kFloat,
                                           Params->GetCompressionLevel());

  // Write dataset parameters
  HDF5_File.writeMatrixDomainType(HDF5_File.getRootGroup(),
                                  HDF5_RootObjectName,
                                  Hdf5File::MatrixDomainType::kReal);
  HDF5_File.writeMatrixDataType  (HDF5_File.getRootGroup(),
                                  HDF5_RootObjectName,
                                  Hdf5File::MatrixDataType::kFloat);


  // Sampled time step
  SampledTimeStep = 0;

  // Set buffer size
  BufferSize = NumberOfSampledElementsPerStep;

  // Allocate memory if needed
  if (!BufferReuse) AllocateMemory();
}// end of Create
//------------------------------------------------------------------------------

/**
 * Reopen the output stream after restart.
 */
void TIndexOutputHDF5Stream::Reopen()
{
  // Get parameters
  TParameters * Params = TParameters::GetInstance();

  // Set buffer size
  BufferSize = SensorMask.GetTotalElementCount();

  // Allocate memory if needed
  if (!BufferReuse) AllocateMemory();

  // Reopen the dataset
  HDF5_DatasetId = HDF5_File.openDataset(HDF5_File.getRootGroup(),
                                         HDF5_RootObjectName);


  if (ReductionOp == roNONE)
  { // raw time series - just seek to the right place in the dataset
    SampledTimeStep = (Params->Get_t_index() < Params->GetStartTimeIndex()) ?
                              0 : (Params->Get_t_index() - Params->GetStartTimeIndex());

  }
  else
  { // aggregated quantities - reload data
    SampledTimeStep = 0;
    // read only if it is necessary (it is anything to read).
    if (Params->Get_t_index() > Params->GetStartTimeIndex())
    {
      // Since there is only a single timestep in the dataset, I can read the whole dataset
      HDF5_File.readCompleteDataset(HDF5_File.getRootGroup(),
                                    HDF5_RootObjectName,
                                    DimensionSizes(BufferSize, 1, 1),
                                    StoreBuffer);
    }
  }
}// end of Reopen
//------------------------------------------------------------------------------


/**
 * Sample grid points, line them up in the buffer an flush to the disk unless a
 * reduction operator is applied.
 */
void TIndexOutputHDF5Stream::Sample()
{
  const float  * SourceData = SourceMatrix.GetRawData();
  const size_t * SensorData = SensorMask.GetRawData();

  switch (ReductionOp)
  {
    case roNONE :
    {
      #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
      for (size_t i = 0; i < BufferSize; i++)
      {
        StoreBuffer[i] = SourceData[SensorData[i]];
      }
      // only raw time series are flushed down to the disk every time step
      FlushBufferToFile();
      /* - for future use when offloading the sampling work to HDF5 - now it seem to be slower
      HDF5_File.WriteSensorbyMaskToHyperSlab(HDF5_DatasetId,
                                             Position,        // position in the dataset
                                             BufferSize,      // number of elements sampled
                                             SensorMask.GetRawData(), // Sensor
                                             SourceMatrix.GetDimensionSizes(), // Matrix dims
                                             SourceMatrix.GetRawData());
      Position.Y++;
       */
      break;
    }// case roNONE

    case roRMS  :
    {
      #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
      for (size_t i = 0; i < BufferSize; i++)
      {
        StoreBuffer[i] += (SourceData[SensorData[i]] * SourceData[SensorData[i]]);
      }
      break;
    }// case roRMS

    case roMAX  :
    {
      #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
      for (size_t i = 0; i < BufferSize; i++)
      {
        if (StoreBuffer[i] < SourceData[SensorData[i]])
          StoreBuffer[i] = SourceData[SensorData[i]];
      }
      break;
    }// case roMAX

    case roMIN  :
    {
      #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
      for (size_t i = 0; i < BufferSize; i++)
      {
        if (StoreBuffer[i] > SourceData[SensorData[i]])
          StoreBuffer[i] = SourceData[SensorData[i]];
      }
      break;
    } //case roMIN
  }// switch
}// end of Sample
//------------------------------------------------------------------------------

/**
 * Apply post-processing on the buffer and flush it to the file.
 */
void TIndexOutputHDF5Stream::PostProcess()
{
  // run inherited method
  TBaseOutputHDF5Stream::PostProcess();
  // When no reduction operator is applied, the data is flushed after every time step
  if (ReductionOp != roNONE) FlushBufferToFile();
}// end of PostProcessing
//------------------------------------------------------------------------------


/**
 * Checkpoint the stream and close.
 */
void TIndexOutputHDF5Stream::Checkpoint()
{
  // raw data has already been flushed, others has to be fushed here
  if (ReductionOp != roNONE) FlushBufferToFile();

}// end of Checkpoint
//------------------------------------------------------------------------------

/**
 * Close stream (apply post-processing if necessary, flush data and close).
 */
void TIndexOutputHDF5Stream::Close()
{
  // the dataset is still opened
  if (HDF5_DatasetId != H5I_BADID)
  {
    HDF5_File.closeDataset(HDF5_DatasetId);
  }

  HDF5_DatasetId = H5I_BADID;
}// end of Close
//------------------------------------------------------------------------------

//----------------------------------------------------------------------------//
//                 TIndexOutputHDF5Stream implementation                      //
//                            protected methods                               //
//----------------------------------------------------------------------------//


/**
 * Flush the buffer down to the file at the actual position.
 */
void TIndexOutputHDF5Stream::FlushBufferToFile()
{
  HDF5_File.writeHyperSlab(HDF5_DatasetId,
                           DimensionSizes(0,SampledTimeStep,0),
                           DimensionSizes(BufferSize,1,1),
                           StoreBuffer);
  SampledTimeStep++;
}// end of FlushToFile
//------------------------------------------------------------------------------






//----------------------------------------------------------------------------//
//                TCuboidOutputHDF5Stream implementation                      //
//                              public methods                                //
//----------------------------------------------------------------------------//


/**
 * Constructor - links the HDF5 dataset, SourceMatrix, and SensorMask together.
 * @param [in] HDF5_File       - HDF5 file to write the output to
 * @param [in] HDF5_GroupName  - The name of the HDF5 group. This group contains datasets for particular cuboids
 * @param [in] SourceMatrix    - Source matrix to be sampled
 * @param [in] SensorMask      - Sensor mask with the cuboid coordinates
 * @param [in] ReductionOp     - Reduction operator
 * @param [in] BufferToReuse   - If there is a memory space to be reused, provide a pointer
 */
TCuboidOutputHDF5Stream::TCuboidOutputHDF5Stream(Hdf5File &             HDF5_File,
                                                 const char *             HDF5_GroupName,
                                                 const TRealMatrix &      SourceMatrix,
                                                 const TIndexMatrix &     SensorMask,
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
 * Destructor.
 * if the file is still opened, it applies the post processing and flush the data.
 * Then, the object memory is freed and the object destroyed.
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
 * with particular datasets (one per cuboid).
 */
void TCuboidOutputHDF5Stream::Create()
{
  // Create the HDF5 group and open it
  HDF5_GroupId = HDF5_File.createGroup(HDF5_File.getRootGroup(), HDF5_RootObjectName);

  // Create all datasets (sizes, chunks, and attributes)
  size_t NumberOfCuboids        = SensorMask.GetDimensionSizes().ny;
  CuboidsInfo.reserve(NumberOfCuboids);
  size_t ActualPositionInBuffer = 0;

  for (size_t CuboidIndex = 0; CuboidIndex < NumberOfCuboids; CuboidIndex++)
  {
    TCuboidInfo CuboidInfo;

    CuboidInfo.HDF5_CuboidId = CreateCuboidDataset(CuboidIndex);
    CuboidInfo.StartingPossitionInBuffer = ActualPositionInBuffer;
    CuboidsInfo.push_back(CuboidInfo);

    ActualPositionInBuffer += (SensorMask.GetBottomRightCorner(CuboidIndex) - SensorMask.GetTopLeftCorner(CuboidIndex)).nElements();
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
 * Reopen the output stream after restart and reload data.
 */
void TCuboidOutputHDF5Stream::Reopen()
{
  // Get parameters
  TParameters * Params = TParameters::GetInstance();

  SampledTimeStep = 0;
  if (ReductionOp == roNONE) // set correct sampled timestep for raw data series
  {
    SampledTimeStep = (Params->Get_t_index() < Params->GetStartTimeIndex()) ?
                        0 : (Params->Get_t_index() - Params->GetStartTimeIndex());
  }

  // Create the memory buffer if necessary and set starting address
  BufferSize = SensorMask.GetTotalNumberOfElementsInAllCuboids();

  // Allocate memory if needed
  if (!BufferReuse) AllocateMemory();


  // Open all datasets (sizes, chunks, and attributes)
  size_t NumberOfCuboids        = SensorMask.GetDimensionSizes().ny;
  CuboidsInfo.reserve(NumberOfCuboids);
  size_t ActualPositionInBuffer = 0;

  // Open the HDF5 group
  HDF5_GroupId = HDF5_File.openGroup(HDF5_File.getRootGroup(), HDF5_RootObjectName);

  for (size_t CuboidIndex = 0; CuboidIndex < NumberOfCuboids; CuboidIndex++)
  {
    TCuboidInfo CuboidInfo;

    // @todo: Can be done easily with std::to_string and c++0x or c++-11
    char HDF5_DatasetName[32] = "";
    // Indexed from 1
    sprintf(HDF5_DatasetName, "%ld",CuboidIndex + 1);

    // open the dataset
    CuboidInfo.HDF5_CuboidId = HDF5_File.openDataset(HDF5_GroupId,
                                                     HDF5_DatasetName);
    CuboidInfo.StartingPossitionInBuffer = ActualPositionInBuffer;
    CuboidsInfo.push_back(CuboidInfo);

    // read only if there is anything to read
    if (Params->Get_t_index() > Params->GetStartTimeIndex())
    {
      if (ReductionOp != roNONE)
      { // Reload data
        DimensionSizes CuboidSize((SensorMask.GetBottomRightCorner(CuboidIndex) - SensorMask.GetTopLeftCorner(CuboidIndex)).nx,
                                   (SensorMask.GetBottomRightCorner(CuboidIndex) - SensorMask.GetTopLeftCorner(CuboidIndex)).ny,
                                   (SensorMask.GetBottomRightCorner(CuboidIndex) - SensorMask.GetTopLeftCorner(CuboidIndex)).nz);

        HDF5_File.readCompleteDataset(HDF5_GroupId,
                                      HDF5_DatasetName,
                                      CuboidSize,
                                      StoreBuffer + ActualPositionInBuffer);
      }
    }
    // move the pointer for the next cuboid beginning (this inits the locations)
    ActualPositionInBuffer += (SensorMask.GetBottomRightCorner(CuboidIndex) -
                               SensorMask.GetTopLeftCorner(CuboidIndex)).nElements();
  }
}// end of Reopen
//------------------------------------------------------------------------------

/**
 * Sample data into buffer and apply reduction, or flush to disk - based on a sensor mask.
 */
void TCuboidOutputHDF5Stream::Sample()
{
  const size_t XY_Size = SourceMatrix.GetDimensionSizes().ny * SourceMatrix.GetDimensionSizes().nx;
  const size_t X_Size  = SourceMatrix.GetDimensionSizes().nx;

  const float * SourceData = SourceMatrix.GetRawData();

  switch (ReductionOp)
  {
    case roNONE :
    {
      /* We use here direct HDF5 offload using MEMSPACE - seems to be faster for bigger datasets*/
      DimensionSizes DatasetPosition(0,0,0,0); //4D position in the dataset
      DimensionSizes CuboidSize(0,0,0,0);      // Size of the cuboid

      DatasetPosition.nt = SampledTimeStep;
      const float * MatrixData = SourceMatrix.GetRawData();

      // iterate over all cuboid to be sampled
      for (size_t CuboidIndex = 0; CuboidIndex < CuboidsInfo.size(); CuboidIndex++)
      {
        CuboidSize = SensorMask.GetBottomRightCorner(CuboidIndex) - SensorMask.GetTopLeftCorner(CuboidIndex);
        CuboidSize.nt = 1;

        HDF5_File.writeCuboidToHyperSlab(CuboidsInfo[CuboidIndex].HDF5_CuboidId,
                                         DatasetPosition,
                                         SensorMask.GetTopLeftCorner(CuboidIndex), // position in the SourceMatrix
                                         CuboidSize,
                                         SourceMatrix.GetDimensionSizes(),
                                         MatrixData);
      }
      SampledTimeStep++;   // Move forward in time

      break;

    /* At the time being, this version using manual data lining up seems to be slower, and is not used
      size_t BufferStart = 0;

      for (size_t CuboidIdx = 0; CuboidIdx < SensorMask.GetDimensionSizes().Y; CuboidIdx++)
      {
        const TDimensionSizes TopLeftCorner     = SensorMask.GetTopLeftCorner(CuboidIdx);
        const TDimensionSizes BottomRightCorner = SensorMask.GetBottomRightCorner(CuboidIdx);

        size_t cuboid_XY_plane_size = ( BottomRightCorner.Y - TopLeftCorner.Y) * ( BottomRightCorner.X - TopLeftCorner.X);
        size_t cuboid_X_plane_size  = ( BottomRightCorner.X - TopLeftCorner.X);

        #pragma omp parallel for collapse(3)
        for (size_t z = TopLeftCorner.Z; z <= BottomRightCorner.Z; z++)
          for (size_t y = TopLeftCorner.Y; y <= BottomRightCorner.Y; y++)
            for (size_t x = TopLeftCorner.X; x <= BottomRightCorner.X; x++)
            {
              StoreBuffer[BufferStart +
                            (z - TopLeftCorner.Z) * cuboid_XY_plane_size +
                            (y - TopLeftCorner.Y) * cuboid_X_plane_size  +
                            (x - TopLeftCorner.X) ]

                    = SourceMatrix[z * XY_Size + y * X_Size + x];
            }
        BufferStart += (BottomRightCorner - TopLeftCorner).GetElementCount();
      }
      FlushBufferToFile();
      break;
        */
    }// case roNONE

    case roRMS  :
    {
      size_t CuboidInBufferStart = 0;

      // Parallelise within the cuboid - Since a typical number of cuboids is 1, than we have to paralelise inside
      for (size_t CuboidIdx = 0; CuboidIdx < SensorMask.GetDimensionSizes().ny; CuboidIdx++)
      {
        const DimensionSizes TopLeftCorner     = SensorMask.GetTopLeftCorner(CuboidIdx);
        const DimensionSizes BottomRightCorner = SensorMask.GetBottomRightCorner(CuboidIdx);

        size_t cuboid_XY_plane_size = (BottomRightCorner.ny - TopLeftCorner.ny + 1) *
                                      (BottomRightCorner.nx - TopLeftCorner.nx + 1);
        size_t cuboid_X_plane_size  = (BottomRightCorner.nx - TopLeftCorner.nx + 1);

        #pragma omp parallel for collapse(3) \
                if ((BottomRightCorner - TopLeftCorner).nElements() > MinGridpointsToSampleInParallel)
        for (size_t z = TopLeftCorner.nz; z <= BottomRightCorner.nz; z++)
          for (size_t y = TopLeftCorner.ny; y <= BottomRightCorner.ny; y++)
            for (size_t x = TopLeftCorner.nx; x <= BottomRightCorner.nx; x++)
            {
              const size_t StoreBufferIndex = CuboidInBufferStart +
                                              (z - TopLeftCorner.nz) * cuboid_XY_plane_size +
                                              (y - TopLeftCorner.ny) * cuboid_X_plane_size  +
                                              (x - TopLeftCorner.nx);

              const size_t SourceIndex = z * XY_Size + y * X_Size + x;

              // aggregating the sum over timesteps
              StoreBuffer[StoreBufferIndex] += (SourceData[SourceIndex] * SourceData[SourceIndex]);
            }

        CuboidInBufferStart += (BottomRightCorner - TopLeftCorner).nElements();
      }

      break;
    }// case roRMS

    case roMAX  :
    {
      size_t CuboidInBufferStart = 0;

      // Parallelise within the cuboid - Since a typical number of cuboids is 1, than we have to paralelise inside.
      for (size_t CuboidIdx = 0; CuboidIdx < SensorMask.GetDimensionSizes().ny; CuboidIdx++)
      {
        const DimensionSizes TopLeftCorner     = SensorMask.GetTopLeftCorner(CuboidIdx);
        const DimensionSizes BottomRightCorner = SensorMask.GetBottomRightCorner(CuboidIdx);

        size_t cuboid_XY_plane_size = (BottomRightCorner.ny - TopLeftCorner.ny + 1) *
                                      (BottomRightCorner.nx - TopLeftCorner.nx + 1);
        size_t cuboid_X_plane_size  = (BottomRightCorner.nx - TopLeftCorner.nx + 1);

        #pragma omp parallel for collapse(3) \
                if ((BottomRightCorner - TopLeftCorner).nElements() > MinGridpointsToSampleInParallel)
        for (size_t z = TopLeftCorner.nz; z <= BottomRightCorner.nz; z++)
          for (size_t y = TopLeftCorner.ny; y <= BottomRightCorner.ny; y++)
            for (size_t x = TopLeftCorner.nx; x <= BottomRightCorner.nx; x++)
            {
              const size_t StoreBufferIndex = CuboidInBufferStart +
                                              (z - TopLeftCorner.nz) * cuboid_XY_plane_size +
                                              (y - TopLeftCorner.ny) * cuboid_X_plane_size  +
                                              (x - TopLeftCorner.nx);

              const size_t SourceIndex = z * XY_Size + y * X_Size + x;

              // finding max
              if (StoreBuffer[StoreBufferIndex] < SourceData[SourceIndex])
              {
                StoreBuffer[StoreBufferIndex] = SourceData[SourceIndex];
              }
            }
        CuboidInBufferStart += (BottomRightCorner - TopLeftCorner).nElements();
      }
      break;
    }// case roMAX

    case roMIN  :
    {
      size_t CuboidInBufferStart = 0;

      // Parallelise within the cuboid - Since a typical number of cuboids is 1, than we have to paralelise inside.
      for (size_t CuboidIdx = 0; CuboidIdx < SensorMask.GetDimensionSizes().ny; CuboidIdx++)
      {
        const DimensionSizes TopLeftCorner     = SensorMask.GetTopLeftCorner(CuboidIdx);
        const DimensionSizes BottomRightCorner = SensorMask.GetBottomRightCorner(CuboidIdx);

        size_t cuboid_XY_plane_size = (BottomRightCorner.ny - TopLeftCorner.ny + 1) *
                                      (BottomRightCorner.nx - TopLeftCorner.nx + 1);
        size_t cuboid_X_plane_size  = (BottomRightCorner.nx - TopLeftCorner.nx + 1);

        #pragma omp parallel for collapse(3) \
                if ((BottomRightCorner - TopLeftCorner).nElements() > MinGridpointsToSampleInParallel)
        for (size_t z = TopLeftCorner.nz; z <= BottomRightCorner.nz; z++)
          for (size_t y = TopLeftCorner.ny; y <= BottomRightCorner.ny; y++)
            for (size_t x = TopLeftCorner.nx; x <= BottomRightCorner.nx; x++)
            {
              const size_t StoreBufferIndex = CuboidInBufferStart +
                                              (z - TopLeftCorner.nz) * cuboid_XY_plane_size +
                                              (y - TopLeftCorner.ny) * cuboid_X_plane_size  +
                                              (x - TopLeftCorner.nx);
              const size_t SourceIndex = z * XY_Size + y * X_Size + x;

              // finding min
              if (StoreBuffer[StoreBufferIndex] > SourceData[SourceIndex])
              {
                StoreBuffer[StoreBufferIndex] = SourceData[SourceIndex];
              }
            }
        CuboidInBufferStart += (BottomRightCorner - TopLeftCorner).nElements();
      }

      break;
    }// case roMIN
  }// switch
}// end of Sample
//------------------------------------------------------------------------------


/**
 * Apply post-processing on the buffer and flush it to the file.
 */
void TCuboidOutputHDF5Stream::PostProcess()
{
  // run inherited method
  TBaseOutputHDF5Stream::PostProcess();
  // When no reduction operator is applied, the data is flushed after every time step
  if (ReductionOp != roNONE) FlushBufferToFile();
}// end of PostProcessing
//------------------------------------------------------------------------------

/**
 * Checkpoint the stream.
 */
void TCuboidOutputHDF5Stream::Checkpoint()
{
  // raw data has already been flushed, others has to be fushed here
  if (ReductionOp != roNONE) FlushBufferToFile();
}// end of Checkpoint
//------------------------------------------------------------------------------


/**
 * Close stream (apply post-processing if necessary, flush data, close datasets
 * and the group).
 */
void TCuboidOutputHDF5Stream::Close()
{
  // the group is still open
  if (HDF5_GroupId != H5I_BADID)
  {
    // Close all datasets and the group
    for (size_t CuboidIndex = 0; CuboidIndex < CuboidsInfo.size(); CuboidIndex++)
    {
      HDF5_File.closeDataset(CuboidsInfo[CuboidIndex].HDF5_CuboidId);
    }
    CuboidsInfo.clear();

    HDF5_File.closeGroup(HDF5_GroupId);
    HDF5_GroupId = H5I_BADID;
  }// if opened
}// end of Close
//------------------------------------------------------------------------------



//----------------------------------------------------------------------------//
//                 TIndexOutputHDF5Stream implementation                      //
//                            protected methods                               //
//----------------------------------------------------------------------------//

/**
 *  Create a new dataset for a given cuboid specified by index (order).
 * @param [in] Index - Index of the cuboid in the sensor mask
 * @return HDF5 handle to the dataset.
 */
hid_t TCuboidOutputHDF5Stream::CreateCuboidDataset(const size_t Index)
{
  TParameters * Params = TParameters::GetInstance();

  // if time series then Number of steps else 1
  size_t NumberOfSampledTimeSteps = (ReductionOp == roNONE)
                                      ? Params->Get_Nt() - Params->GetStartTimeIndex()
                                      : 0; // will be a 3D dataset
  // Set cuboid dimensions (subtract two corners (add 1) and use the appropriate component)
  DimensionSizes CuboidSize((SensorMask.GetBottomRightCorner(Index) - SensorMask.GetTopLeftCorner(Index)).nx,
                             (SensorMask.GetBottomRightCorner(Index) - SensorMask.GetTopLeftCorner(Index)).ny,
                             (SensorMask.GetBottomRightCorner(Index) - SensorMask.GetTopLeftCorner(Index)).nz,
                             NumberOfSampledTimeSteps
                            );

  // Set chunk size
  // If the size of the cuboid is bigger than 32 MB per timestep, set the chunk to approx 4MB
  size_t NumberOfSlabs = 1; //at least one slab
  DimensionSizes CuboidChunkSize(CuboidSize.nx, CuboidSize.ny, CuboidSize.nz, (ReductionOp == roNONE) ? 1 : 0);

  if (CuboidChunkSize.nElements() > (ChunkSize_4MB * 8))
  {
    while (NumberOfSlabs * CuboidSize.nx * CuboidSize.ny < ChunkSize_4MB) NumberOfSlabs++;
    CuboidChunkSize.nz = NumberOfSlabs;
  }

  // @todo: Can be done easily with std::to_string and c++0x or c++-11
  char HDF5_DatasetName[32] = "";
  // Indexed from 1
  sprintf(HDF5_DatasetName, "%ld",Index+1);
  hid_t HDF5_DatasetId = HDF5_File.createDataset(HDF5_GroupId,
                                                 HDF5_DatasetName,
                                                 CuboidSize,
                                                 CuboidChunkSize,
                                                 Hdf5File::MatrixDataType::kFloat,
                                                 Params->GetCompressionLevel());

  // Write dataset parameters
  HDF5_File.writeMatrixDomainType(HDF5_GroupId,
                                  HDF5_DatasetName,
                                  Hdf5File::MatrixDomainType::kReal);
  HDF5_File.writeMatrixDataType  (HDF5_GroupId,
                                  HDF5_DatasetName,
                                  Hdf5File::MatrixDataType::kFloat);

  return HDF5_DatasetId;
}//end of CreateCuboidDatasets
//------------------------------------------------------------------------------



/**
 * Flush the buffer to the file (to multiple datasets if necessary).
 */
void TCuboidOutputHDF5Stream::FlushBufferToFile()
{
  DimensionSizes Position (0,0,0,0);
  DimensionSizes BlockSize(0,0,0,0);

  if (ReductionOp == roNONE) Position.nt = SampledTimeStep;

  for (size_t CuboidIndex = 0; CuboidIndex < CuboidsInfo.size(); CuboidIndex++)
  {

    BlockSize = SensorMask.GetBottomRightCorner(CuboidIndex) - SensorMask.GetTopLeftCorner(CuboidIndex);
    BlockSize.nt = 1;

    HDF5_File.writeHyperSlab(CuboidsInfo[CuboidIndex].HDF5_CuboidId,
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
 * Constructor - links the HDF5 dataset and SourceMatrix.
 * @param [in] HDF5_File        - HDF5 file to write the output to
 * @param [in] HDF5_DatasetName - The name of the HDF5 group. This group contains datasets for particular cuboids
 * @param [in] SourceMatrix     - Source matrix to be sampled
 * @param [in] ReductionOp      - Reduction operator
 * @param [in] BufferToReuse    - If there is a memory space to be reused, provide a pointer
 */
TWholeDomainOutputHDF5Stream::TWholeDomainOutputHDF5Stream(Hdf5File &             HDF5_File,
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
 * Destructor.
 * if the file is still opened, it applies the post processing and flush the data.
 * Then, the object memory is freed and the object destroyed.
 */
TWholeDomainOutputHDF5Stream::~TWholeDomainOutputHDF5Stream()
{
  Close();
  // free memory only if it was allocated
  if (!BufferReuse) FreeMemory();
}// end of Destructor
//------------------------------------------------------------------------------

/**
 * Create a HDF5 stream for the whole domain and allocate data for it.
 */
void TWholeDomainOutputHDF5Stream::Create()
{
  DimensionSizes ChunkSize(SourceMatrix.GetDimensionSizes().nx, SourceMatrix.GetDimensionSizes().ny, 1);

  // Create a dataset under the root group
  HDF5_DatasetId = HDF5_File.createDataset(HDF5_File.getRootGroup(),
                                           HDF5_RootObjectName,
                                           SourceMatrix.GetDimensionSizes(),
                                           ChunkSize,
                                           Hdf5File::MatrixDataType::kFloat,
                                           TParameters::GetInstance()->GetCompressionLevel());

  // Write dataset parameters
  HDF5_File.writeMatrixDomainType(HDF5_File.getRootGroup(),
                                  HDF5_RootObjectName,
                                  Hdf5File::MatrixDomainType::kReal);
  HDF5_File.writeMatrixDataType  (HDF5_File.getRootGroup(),
                                  HDF5_RootObjectName,
                                  Hdf5File::MatrixDataType::kFloat);

  // Set buffer size
  BufferSize = SourceMatrix.GetTotalElementCount();

  // Allocate memory if needed
  if (!BufferReuse) AllocateMemory();
}//end of Create
//------------------------------------------------------------------------------


/**
 * Reopen the output stream after restart and reload data.
 */
void TWholeDomainOutputHDF5Stream::Reopen()
{
  TParameters * Params = TParameters::GetInstance();

  // Set buffer size
  BufferSize = SourceMatrix.GetTotalElementCount();

  // Allocate memory if needed
  if (!BufferReuse) AllocateMemory();

  // Open the dataset under the root group
  HDF5_DatasetId = HDF5_File.openDataset(HDF5_File.getRootGroup(),
                                         HDF5_RootObjectName);

  SampledTimeStep = 0;
  if (ReductionOp == roNONE)
  { // seek in the dataset
    SampledTimeStep = (Params->Get_t_index() < Params->GetStartTimeIndex()) ?
                        0 : (Params->Get_t_index() - Params->GetStartTimeIndex());
  }
  else
  { // reload data
    if (Params->Get_t_index() > Params->GetStartTimeIndex())
    {
      HDF5_File.readCompleteDataset(HDF5_File.getRootGroup(),
                                    HDF5_RootObjectName,
                                    SourceMatrix.GetDimensionSizes(),
                                    StoreBuffer);
    }
  }
}// end of Reopen
//------------------------------------------------------------------------------


/**
 * Sample all grid points, line them up in the buffer an flush to the disk unless
 * a reduction operator is applied.
 */
void TWholeDomainOutputHDF5Stream::Sample()
{
  const float * SourceData = SourceMatrix.GetRawData();

  switch (ReductionOp)
  {
    case roNONE :
    {

       /* We use here direct HDF5 offload using MEMSPACE - seems to be faster for bigger datasets*/
      const DimensionSizes DatasetPosition(0,0,0,SampledTimeStep); //4D position in the dataset

      DimensionSizes CuboidSize(SourceMatrix.GetDimensionSizes());// Size of the cuboid
      CuboidSize.nt = 1;

      // iterate over all cuboid to be sampled
      HDF5_File.writeCuboidToHyperSlab(HDF5_DatasetId,
                                       DatasetPosition,
                                       DimensionSizes(0,0,0,0), // position in the SourceMatrix
                                       CuboidSize,
                                       SourceMatrix.GetDimensionSizes(),
                                       SourceMatrix.GetRawData());

      SampledTimeStep++;   // Move forward in time

      // This version is obsolete (needs one more data movement)
      /*
       #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
      for (size_t i = 0; i < BufferSize; i++)
      {
        StoreBuffer[i] = SourceData[i];
      }
      // only raw time series are flushed down to the disk every time step
      FlushBufferToFile();
      */

      break;
    }// case roNONE

    case roRMS  :
    {
      #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
      for (size_t i = 0; i < BufferSize; i++)
      {
        StoreBuffer[i] += (SourceData[i] * SourceData[i]);
      }
      break;
    }// case roRMS

    case roMAX  :
    {
      #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
      for (size_t i = 0; i < BufferSize; i++)
      {
        if (StoreBuffer[i] < SourceData[i])  StoreBuffer[i] = SourceData[i];
      }
      break;
    }//case roMAX

    case roMIN  :
    {
      #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
      for (size_t i = 0; i < BufferSize; i++)
      {
        if (StoreBuffer[i] > SourceData[i]) StoreBuffer[i] = SourceData[i];
      }
      break;
    } //case roMIN
  }// switch
}// end of Sample
//------------------------------------------------------------------------------


/**
 * Apply post-processing on the buffer and flush it to the file.
 */
void TWholeDomainOutputHDF5Stream::PostProcess()
{
  // run inherited method
  TBaseOutputHDF5Stream::PostProcess();
  // When no reduction operator is applied, the data is flushed after every time step
  if (ReductionOp != roNONE) FlushBufferToFile();
}// end of PostProcessing
//------------------------------------------------------------------------------

/**
 * Checkpoint the stream
 */
void TWholeDomainOutputHDF5Stream::Checkpoint()
{
  // raw data has already been flushed, others has to be flushed here.
  if (ReductionOp != roNONE) FlushBufferToFile();
}// end of Checkpoint
//------------------------------------------------------------------------------

/**
 * Close stream (apply post-processing if necessary, flush data and close).
 */
void TWholeDomainOutputHDF5Stream::Close()
{
  // the dataset is still opened
  if (HDF5_DatasetId != H5I_BADID)
  {
    HDF5_File.closeDataset(HDF5_DatasetId);
  }

  HDF5_DatasetId = H5I_BADID;
}// end of Close
//------------------------------------------------------------------------------


//----------------------------------------------------------------------------//
//                TWholeDomainOutputHDF5Stream implementation                 //
//                            protected methods                               //
//----------------------------------------------------------------------------//


/**
 * Flush the buffer down to the file at the actual position.
 */
void TWholeDomainOutputHDF5Stream::FlushBufferToFile()
{
  DimensionSizes Size = SourceMatrix.GetDimensionSizes();
  DimensionSizes Position(0,0,0);

  // Not used for roNONE now!
  if (ReductionOp == roNONE)
  {
    Position.nt = SampledTimeStep;
    Size.nt = SampledTimeStep;
  }

  HDF5_File.writeHyperSlab(HDF5_DatasetId,
                           Position,
                           Size,
                           StoreBuffer);
  SampledTimeStep++;
}// end of FlushToFile
//------------------------------------------------------------------------------

