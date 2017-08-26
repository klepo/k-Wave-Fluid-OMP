/**
 * @file        CuboidOutputStream.cpp
 * @author      Jiri Jaros \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The implementation file of classes responsible for storing output quantities based
 *              on the cuboid sensor mask into the output HDF5 file.
 *
 * @version     kspaceFirstOrder3D 2.16
 *
 * @date        26 August    2017, 17:03 (created) \n
 *              26 August    2017, 17:03 (revised)
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
#include <OutputStreams/CuboidOutputStream.h>
#include <Parameters/Parameters.h>



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
                                                 const RealMatrix &      SourceMatrix,
                                                 const IndexMatrix &     SensorMask,
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
  size_t NumberOfCuboids        = SensorMask.getDimensionSizes().ny;
  CuboidsInfo.reserve(NumberOfCuboids);
  size_t ActualPositionInBuffer = 0;

  for (size_t CuboidIndex = 0; CuboidIndex < NumberOfCuboids; CuboidIndex++)
  {
    TCuboidInfo CuboidInfo;

    CuboidInfo.HDF5_CuboidId = CreateCuboidDataset(CuboidIndex);
    CuboidInfo.StartingPossitionInBuffer = ActualPositionInBuffer;
    CuboidsInfo.push_back(CuboidInfo);

    ActualPositionInBuffer += (SensorMask.getBottomRightCorner(CuboidIndex) - SensorMask.getTopLeftCorner(CuboidIndex)).nElements();
  }

  //we're at the beginning
  SampledTimeStep = 0;

  // Create the memory buffer if necessary and set starting address
  BufferSize = SensorMask.getSizeOfAllCuboids();

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
  Parameters& Params = Parameters::getInstance();

  SampledTimeStep = 0;
  if (ReductionOp == roNONE) // set correct sampled timestep for raw data series
  {
    SampledTimeStep = (Params.getTimeIndex() < Params.getSamplingStartTimeIndex()) ?
                        0 : (Params.getTimeIndex() - Params.getSamplingStartTimeIndex());
  }

  // Create the memory buffer if necessary and set starting address
  BufferSize = SensorMask.getSizeOfAllCuboids();

  // Allocate memory if needed
  if (!BufferReuse) AllocateMemory();


  // Open all datasets (sizes, chunks, and attributes)
  size_t NumberOfCuboids        = SensorMask.getDimensionSizes().ny;
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
    if (Params.getTimeIndex() > Params.getSamplingStartTimeIndex())
    {
      if (ReductionOp != roNONE)
      { // Reload data
        DimensionSizes CuboidSize((SensorMask.getBottomRightCorner(CuboidIndex) - SensorMask.getTopLeftCorner(CuboidIndex)).nx,
                                   (SensorMask.getBottomRightCorner(CuboidIndex) - SensorMask.getTopLeftCorner(CuboidIndex)).ny,
                                   (SensorMask.getBottomRightCorner(CuboidIndex) - SensorMask.getTopLeftCorner(CuboidIndex)).nz);

        HDF5_File.readCompleteDataset(HDF5_GroupId,
                                      HDF5_DatasetName,
                                      CuboidSize,
                                      StoreBuffer + ActualPositionInBuffer);
      }
    }
    // move the pointer for the next cuboid beginning (this inits the locations)
    ActualPositionInBuffer += (SensorMask.getBottomRightCorner(CuboidIndex) -
                               SensorMask.getTopLeftCorner(CuboidIndex)).nElements();
  }
}// end of Reopen
//------------------------------------------------------------------------------

/**
 * Sample data into buffer and apply reduction, or flush to disk - based on a sensor mask.
 */
void TCuboidOutputHDF5Stream::Sample()
{
  const size_t XY_Size = SourceMatrix.getDimensionSizes().ny * SourceMatrix.getDimensionSizes().nx;
  const size_t X_Size  = SourceMatrix.getDimensionSizes().nx;

  const float * SourceData = SourceMatrix.getData();

  switch (ReductionOp)
  {
    case roNONE :
    {
      /* We use here direct HDF5 offload using MEMSPACE - seems to be faster for bigger datasets*/
      DimensionSizes DatasetPosition(0,0,0,0); //4D position in the dataset
      DimensionSizes CuboidSize(0,0,0,0);      // Size of the cuboid

      DatasetPosition.nt = SampledTimeStep;
      const float * MatrixData = SourceMatrix.getData();

      // iterate over all cuboid to be sampled
      for (size_t CuboidIndex = 0; CuboidIndex < CuboidsInfo.size(); CuboidIndex++)
      {
        CuboidSize = SensorMask.getBottomRightCorner(CuboidIndex) - SensorMask.getTopLeftCorner(CuboidIndex);
        CuboidSize.nt = 1;

        HDF5_File.writeCuboidToHyperSlab(CuboidsInfo[CuboidIndex].HDF5_CuboidId,
                                         DatasetPosition,
                                         SensorMask.getTopLeftCorner(CuboidIndex), // position in the SourceMatrix
                                         CuboidSize,
                                         SourceMatrix.getDimensionSizes(),
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
      for (size_t CuboidIdx = 0; CuboidIdx < SensorMask.getDimensionSizes().ny; CuboidIdx++)
      {
        const DimensionSizes TopLeftCorner     = SensorMask.getTopLeftCorner(CuboidIdx);
        const DimensionSizes BottomRightCorner = SensorMask.getBottomRightCorner(CuboidIdx);

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
      for (size_t CuboidIdx = 0; CuboidIdx < SensorMask.getDimensionSizes().ny; CuboidIdx++)
      {
        const DimensionSizes TopLeftCorner     = SensorMask.getTopLeftCorner(CuboidIdx);
        const DimensionSizes BottomRightCorner = SensorMask.getBottomRightCorner(CuboidIdx);

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
      for (size_t CuboidIdx = 0; CuboidIdx < SensorMask.getDimensionSizes().ny; CuboidIdx++)
      {
        const DimensionSizes TopLeftCorner     = SensorMask.getTopLeftCorner(CuboidIdx);
        const DimensionSizes BottomRightCorner = SensorMask.getBottomRightCorner(CuboidIdx);

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
  Parameters& Params = Parameters::getInstance();

  // if time series then Number of steps else 1
  size_t NumberOfSampledTimeSteps = (ReductionOp == roNONE)
                                      ? Params.getNt() - Params.getSamplingStartTimeIndex()
                                      : 0; // will be a 3D dataset
  // Set cuboid dimensions (subtract two corners (add 1) and use the appropriate component)
  DimensionSizes CuboidSize((SensorMask.getBottomRightCorner(Index) - SensorMask.getTopLeftCorner(Index)).nx,
                             (SensorMask.getBottomRightCorner(Index) - SensorMask.getTopLeftCorner(Index)).ny,
                             (SensorMask.getBottomRightCorner(Index) - SensorMask.getTopLeftCorner(Index)).nz,
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
                                                 Params.getCompressionLevel());

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

    BlockSize = SensorMask.getBottomRightCorner(CuboidIndex) - SensorMask.getTopLeftCorner(CuboidIndex);
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

