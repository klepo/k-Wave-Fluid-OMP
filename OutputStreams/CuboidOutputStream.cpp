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
 *              31 August    2017, 15:25 (revised)
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

#include <algorithm>

#include <OutputStreams/CuboidOutputStream.h>
#include <Parameters/Parameters.h>

//--------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------- Constants -----------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//


//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Public methods ---------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Constructor - links the HDF5 dataset, SourceMatrix, and SensorMask together.
 */
CuboidOutputStream::CuboidOutputStream(Hdf5File&            file,
                                       MatrixName&          groupName,
                                       const RealMatrix&    sourceMatrix,
                                       const IndexMatrix&   sensorMask,
                                       const ReduceOperator reduceOp,
                                       float*               bufferToReuse)
  : BaseOutputStream(file, groupName, sourceMatrix, reduceOp, bufferToReuse),
    mSensorMask(sensorMask),
    mGroup(H5I_BADID),
    mSampledTimeStep(0)
{

}// end of CuboidOutputStream
//----------------------------------------------------------------------------------------------------------------------

/**
 * Destructor.
 */
CuboidOutputStream::~CuboidOutputStream()
{
  close();

  // free memory only if it was allocated
  if (!mBufferReuse) freeMemory();
}// end ~CuboidOutputStream
//----------------------------------------------------------------------------------------------------------------------

/**
 * Create a HDF5 stream and allocate data for it. It also creates a HDF5 group with particular datasets
 * (one per cuboid).
 */
void CuboidOutputStream::create()
{
  // Create the HDF5 group and open it
  mGroup = mFile.createGroup(mFile.getRootGroup(), mRootObjectName);

  // Create all datasets (sizes, chunks, and attributes)
  size_t nCuboids = mSensorMask.getDimensionSizes().ny;
  mCuboidsInfo.reserve(nCuboids);

  size_t actualPositionInBuffer = 0;

  for (size_t cuboidIdx = 0; cuboidIdx < nCuboids; cuboidIdx++)
  {
    CuboidInfo cuboidInfo;

    cuboidInfo.cuboidId = createCuboidDataset(cuboidIdx);
    cuboidInfo.startingPossitionInBuffer = actualPositionInBuffer;
    mCuboidsInfo.push_back(cuboidInfo);

    actualPositionInBuffer += (mSensorMask.getBottomRightCorner(cuboidIdx) -
                               mSensorMask.getTopLeftCorner(cuboidIdx)
                              ).nElements();
  }

  //we're at the beginning
  mSampledTimeStep = 0;

  // Create the memory buffer if necessary and set starting address
  mBufferSize = mSensorMask.getSizeOfAllCuboids();

  // Allocate memory if needed
  if (!mBufferReuse) allocateMemory();
}// end of create
//----------------------------------------------------------------------------------------------------------------------

/**
 * Reopen the output stream after restart and reload data.
 */
void CuboidOutputStream::reopen()
{
  // Get parameters
  const Parameters& params = Parameters::getInstance();

  mSampledTimeStep = 0;
  if (mReduceOp == ReduceOperator::kNone) // set correct sampled timestep for raw data series
  {
    mSampledTimeStep = (params.getTimeIndex() < params.getSamplingStartTimeIndex()) ?
                        0 : (params.getTimeIndex() - params.getSamplingStartTimeIndex());
  }

  // Create the memory buffer if necessary and set starting address
  mBufferSize = mSensorMask.getSizeOfAllCuboids();

  // Allocate memory if needed
  if (!mBufferReuse) allocateMemory();


  // Open all datasets (sizes, chunks, and attributes)
  size_t nCuboids = mSensorMask.getDimensionSizes().ny;
  mCuboidsInfo.reserve(nCuboids);
  size_t actualPositionInBuffer = 0;

  // Open the HDF5 group
  mGroup = mFile.openGroup(mFile.getRootGroup(), mRootObjectName);

  for (size_t cuboidIdx = 0; cuboidIdx < nCuboids; cuboidIdx++)
  {
    CuboidInfo cuboidInfo;

    // Indexed from 1
    const std::string datasetName = std::to_string(cuboidIdx + 1);

    // open the dataset
    cuboidInfo.cuboidId = mFile.openDataset(mGroup, datasetName);
    cuboidInfo.startingPossitionInBuffer = actualPositionInBuffer;
    mCuboidsInfo.push_back(cuboidInfo);

    // read only if there is anything to read
    if (params.getTimeIndex() > params.getSamplingStartTimeIndex())
    {
      if (mReduceOp != ReduceOperator::kNone)
      { // Reload data
        DimensionSizes cuboidSize((mSensorMask.getBottomRightCorner(cuboidIdx) -
                                   mSensorMask.getTopLeftCorner(cuboidIdx)).nx,
                                  (mSensorMask.getBottomRightCorner(cuboidIdx) -
                                   mSensorMask.getTopLeftCorner(cuboidIdx)).ny,
                                  (mSensorMask.getBottomRightCorner(cuboidIdx) -
                                   mSensorMask.getTopLeftCorner(cuboidIdx)).nz);

        mFile.readCompleteDataset(mGroup,
                                  datasetName,
                                  cuboidSize,
                                  mStoreBuffer + actualPositionInBuffer);
      }
    }
    // move the pointer for the next cuboid beginning (this inits the locations)
    actualPositionInBuffer += (mSensorMask.getBottomRightCorner(cuboidIdx) -
                               mSensorMask.getTopLeftCorner(cuboidIdx)).nElements();
  }
}// end of reopen
//----------------------------------------------------------------------------------------------------------------------

/**
 * Sample data into buffer and apply reduction, or flush to disk - based on a sensor mask.
 */
void CuboidOutputStream::sample()
{
  switch (mReduceOp)
  {
    case ReduceOperator::kNone :
    {
      /* We use here direct HDF5 offload using MEMSPACE - seems to be faster for bigger datasets*/
      DimensionSizes datasetPosition(0, 0, 0, 0); //4D position in the dataset
      DimensionSizes cuboidSize(0, 0, 0, 0);      // Size of the cuboid

      datasetPosition.nt = mSampledTimeStep;

      // iterate over all cuboid to be sampled
      for (size_t cuboidIdx = 0; cuboidIdx < mCuboidsInfo.size(); cuboidIdx++)
      {
        cuboidSize    = mSensorMask.getBottomRightCorner(cuboidIdx) - mSensorMask.getTopLeftCorner(cuboidIdx);
        cuboidSize.nt = 1;

        mFile.writeCuboidToHyperSlab(mCuboidsInfo[cuboidIdx].cuboidId,
                                     datasetPosition,
                                     mSensorMask.getTopLeftCorner(cuboidIdx), // position in the SourceMatrix
                                     cuboidSize,
                                     mSourceMatrix.getDimensionSizes(),
                                     mSourceMatrix.getData());
      }
      mSampledTimeStep++;   // Move forward in time

      break;

     // At the time being, this version using manual data lining up seems to be slower, and is not used.
     // sampleAggregated<ReduceOperator::kNone>();
    }// case kNone

    case ReduceOperator::kRms:
    {
      sampleAggregated<ReduceOperator::kRms>();
      break;
    }
    case ReduceOperator::kMax:
    {
      sampleAggregated<ReduceOperator::kMax>();
      break;
    }
    case ReduceOperator::kMin:
    {
      sampleAggregated<ReduceOperator::kMin>();
      break;
    }
  }// switch
}// end of sample
//----------------------------------------------------------------------------------------------------------------------

/**
 * Apply post-processing on the buffer and flush it to the file.
 */
void CuboidOutputStream::postProcess()
{
  // run inherited method
  BaseOutputStream::postProcess();
  // When no reduction operator is applied, the data is flushed after every time step
  if (mReduceOp != ReduceOperator::kNone) flushBufferToFile();
}// end of postProcess
//----------------------------------------------------------------------------------------------------------------------

/**
 * Checkpoint the stream.
 */
void CuboidOutputStream::checkpoint()
{
  // raw data has already been flushed, others has to be flushed here
  if (mReduceOp != ReduceOperator::kNone) flushBufferToFile();
}// end of checkpoint
//----------------------------------------------------------------------------------------------------------------------

/**
 * Close stream (apply post-processing if necessary, flush data, close datasets and the group).
 */
void CuboidOutputStream::close()
{
  // the group is still open
  if (mGroup != H5I_BADID)
  {
    // Close all datasets and the group
    for (size_t cuboidIdx = 0; cuboidIdx < mCuboidsInfo.size(); cuboidIdx++)
    {
      mFile.closeDataset(mCuboidsInfo[cuboidIdx].cuboidId);
    }
    mCuboidsInfo.clear();

    mFile.closeGroup(mGroup);
    mGroup = H5I_BADID;
  }// if opened
}// end of close
//----------------------------------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------ Protected methods -------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Create a new dataset for a given cuboid specified by index (order).
 */
hid_t CuboidOutputStream::createCuboidDataset(const size_t cuboidIdx)
{
  Parameters& params = Parameters::getInstance();

  // if time series then Number of steps else 1
  size_t nSampledTimeSteps = (mReduceOp == ReduceOperator::kNone)
                                ? params.getNt() - params.getSamplingStartTimeIndex()
                                : 0; // will be a 3D dataset
  // Set cuboid dimensions (subtract two corners (add 1) and use the appropriate component)
  DimensionSizes cuboidSize((mSensorMask.getBottomRightCorner(cuboidIdx) - mSensorMask.getTopLeftCorner(cuboidIdx)).nx,
                            (mSensorMask.getBottomRightCorner(cuboidIdx) - mSensorMask.getTopLeftCorner(cuboidIdx)).ny,
                            (mSensorMask.getBottomRightCorner(cuboidIdx) - mSensorMask.getTopLeftCorner(cuboidIdx)).nz,
                             nSampledTimeSteps);

  // Set chunk size
  // If the size of the cuboid is bigger than 32 MB per timestep, set the chunk to approx 4MB
  size_t nSlabs = 1; //at least one slab
  DimensionSizes cuboidChunkSize(cuboidSize.nx,
                                 cuboidSize.ny,
                                 cuboidSize.nz,
                                 (mReduceOp == ReduceOperator::kNone) ? 1 : 0);

  if (cuboidChunkSize.nElements() > (kChunkSize4MB * 8))
  {
    while (nSlabs * cuboidSize.nx * cuboidSize.ny < kChunkSize4MB) nSlabs++;
    cuboidChunkSize.nz = nSlabs;
  }

  // Indexed from 1
  const std::string datasetName = std::to_string(cuboidIdx + 1);

  hid_t dataset = mFile.createDataset(mGroup,
                                     datasetName,
                                     cuboidSize,
                                     cuboidChunkSize,
                                     Hdf5File::MatrixDataType::kFloat,
                                     params.getCompressionLevel());

  // Write dataset parameters
  mFile.writeMatrixDomainType(mGroup, datasetName, Hdf5File::MatrixDomainType::kReal);
  mFile.writeMatrixDataType  (mGroup, datasetName, Hdf5File::MatrixDataType::kFloat);

  return dataset;
}//end of createCuboidDataset
//----------------------------------------------------------------------------------------------------------------------

/**
 * Sample aggregated values.
 */
template<BaseOutputStream::ReduceOperator reduceOp>
void CuboidOutputStream::sampleAggregated()
{
  const size_t slabSize = mSourceMatrix.getDimensionSizes().ny * mSourceMatrix.getDimensionSizes().nx;
  const size_t rowSize  = mSourceMatrix.getDimensionSizes().nx;

  const float* sourceData = mSourceMatrix.getData();

  size_t cuboidInBufferStart = 0;

  // Parallelise within the cuboid - Since a typical number of cuboids is 1, than we have to paralelise inside
  #pragma omp parallel
  for (size_t cuboidIdx = 0; cuboidIdx < mSensorMask.getDimensionSizes().ny; cuboidIdx++)
  {
    const DimensionSizes topLeftCorner     = mSensorMask.getTopLeftCorner(cuboidIdx);
    const DimensionSizes bottomRightCorner = mSensorMask.getBottomRightCorner(cuboidIdx);

    size_t cuboidSlabSize = (bottomRightCorner.ny - topLeftCorner.ny + 1) *
                            (bottomRightCorner.nx - topLeftCorner.nx + 1);
    size_t cuboidRowSize  = (bottomRightCorner.nx - topLeftCorner.nx + 1);

    #pragma omp for collapse(3)
    for (size_t z = topLeftCorner.nz; z <= bottomRightCorner.nz; z++)
      for (size_t y = topLeftCorner.ny; y <= bottomRightCorner.ny; y++)
        for (size_t x = topLeftCorner.nx; x <= bottomRightCorner.nx; x++)
        {
          const size_t storeBufferIndex = cuboidInBufferStart +
                                          (z - topLeftCorner.nz) * cuboidSlabSize +
                                          (y - topLeftCorner.ny) * cuboidRowSize  +
                                          (x - topLeftCorner.nx);

          const size_t sourceIndex = z * slabSize + y * rowSize + x;

          // based on template parameter
          switch (reduceOp)
          {
            case ReduceOperator::kNone:
            {
              mStoreBuffer[storeBufferIndex] = sourceData[sourceIndex];
              break;
            }
            case ReduceOperator::kRms:
            {
              mStoreBuffer[storeBufferIndex] += (sourceData[sourceIndex] * sourceData[sourceIndex]);
              break;
            }
            case ReduceOperator::kMax:
            {
              mStoreBuffer[storeBufferIndex] = std::max(mStoreBuffer[storeBufferIndex], sourceData[sourceIndex]);
              break;
            }
            case ReduceOperator::kMin:
            {
              mStoreBuffer[storeBufferIndex] = std::min(mStoreBuffer[storeBufferIndex], sourceData[sourceIndex]);
              break;
            }
          }// switch
        }
     // must be done only once
     #pragma omp single
     {
       cuboidInBufferStart += (bottomRightCorner - topLeftCorner).nElements();
     }
  }
}// end of sampleAggregated
//----------------------------------------------------------------------------------------------------------------------

/**
 * Flush the buffer to the file (to multiple datasets if necessary).
 */
void CuboidOutputStream::flushBufferToFile()
{
  DimensionSizes position (0, 0, 0, 0);
  DimensionSizes blockSize(0, 0, 0, 0);

  if (mReduceOp == ReduceOperator::kNone) position.nt = mSampledTimeStep;

  for (size_t cuboidIdx = 0; cuboidIdx < mCuboidsInfo.size(); cuboidIdx++)
  {
    blockSize    = mSensorMask.getBottomRightCorner(cuboidIdx) - mSensorMask.getTopLeftCorner(cuboidIdx);
    blockSize.nt = 1;

    mFile.writeHyperSlab(mCuboidsInfo[cuboidIdx].cuboidId,
                         position,
                         blockSize,
                         mStoreBuffer + mCuboidsInfo[cuboidIdx].startingPossitionInBuffer);
  }

  mSampledTimeStep++;
}// end of flushBufferToFile
//----------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Private methods --------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//
