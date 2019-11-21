/**
 * @file      CuboidOutputStream.cpp
 *
 * @author    Jiri Jaros \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            jarosjir@fit.vutbr.cz
 *
 * @brief     The implementation file of classes responsible for storing output quantities based
 *            on the cuboid sensor mask into the output HDF5 file.
 *
 * @version   kspaceFirstOrder 2.17
 *
 * @date      26 August    2017, 17:03 (created) \n
 *            20 February  2019, 14:45 (revised)
 *
 * @copyright Copyright (C) 2019 Jiri Jaros and Bradley Treeby.
 *
 * This file is part of the C++ extension of the [k-Wave Toolbox](http://www.k-wave.org).
 *
 * This file is part of the k-Wave. k-Wave is free software: you can redistribute it and/or modify it under the terms
 * of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with k-Wave.
 * If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).
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
  allocateMinMaxMemory(mSensorMask.getDimensionSizes().ny);
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
  // Create and open or only open the HDF5 group
  if (mFile.groupExists(mFile.getRootGroup(), mRootObjectName))
  {
    mGroup = mFile.openGroup(mFile.getRootGroup(), mRootObjectName);
  }
  else
  {
    mGroup = mFile.createGroup(mFile.getRootGroup(), mRootObjectName);
  }

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

    if (mReduceOp == ReduceOperator::kC)
    {
      actualPositionInBuffer += (mSensorMask.getBottomRightCorner(cuboidIdx) -
                                 mSensorMask.getTopLeftCorner(cuboidIdx)
                                ).nElements() * mCompressHelper->getHarmonics() * 2;
    }
    else
    {
      actualPositionInBuffer += (mSensorMask.getBottomRightCorner(cuboidIdx) -
                                 mSensorMask.getTopLeftCorner(cuboidIdx)
                                ).nElements();
    }
  }

  // We're at the beginning
  mSampledTimeStep = 0;

  // Set buffer size
  if (mReduceOp == ReduceOperator::kC)
  {
    mBufferSize = mSensorMask.getSizeOfAllCuboids() * mCompressHelper->getHarmonics() * 2;
  }
  else
  {
    mBufferSize = mSensorMask.getSizeOfAllCuboids();
  }

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
  if (mReduceOp == ReduceOperator::kNone || mReduceOp == ReduceOperator::kC) // set correct sampled timestep for raw data series
  {
    mSampledTimeStep = (params.getTimeIndex() < params.getSamplingStartTimeIndex()) ?
                        0 : (params.getTimeIndex() - params.getSamplingStartTimeIndex());
    if (mReduceOp == ReduceOperator::kC)
    {
      mCompressedTimeStep = size_t(std::max(float(floor(float(mSampledTimeStep) / mCompressHelper->getOSize())), 0.0f));
    }
  }

  // Create the memory buffer if necessary and set starting address
  if (mReduceOp == ReduceOperator::kC)
  {
    mBufferSize = mSensorMask.getSizeOfAllCuboids() * mCompressHelper->getHarmonics() * 2;
  }
  else
  {
    mBufferSize = mSensorMask.getSizeOfAllCuboids();
  }

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
    const std::string datasetName = (mReduceOp == ReduceOperator::kC) ? std::to_string(cuboidIdx + 1) + "_c" : std::to_string(cuboidIdx + 1);

    // open the dataset
    cuboidInfo.cuboidId = mFile.openDataset(mGroup, datasetName);
    cuboidInfo.startingPossitionInBuffer = actualPositionInBuffer;
    mCuboidsInfo.push_back(cuboidInfo);

    // read only if there is anything to read
    if (params.getTimeIndex() > params.getSamplingStartTimeIndex())
    {
      if (mReduceOp != ReduceOperator::kNone && mReduceOp != ReduceOperator::kC)
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
    if (mReduceOp == ReduceOperator::kC)
    {
      actualPositionInBuffer += (mSensorMask.getBottomRightCorner(cuboidIdx) -
                                 mSensorMask.getTopLeftCorner(cuboidIdx)
                                ).nElements() * mCompressHelper->getHarmonics() * 2;
    }
    else
    {
      actualPositionInBuffer += (mSensorMask.getBottomRightCorner(cuboidIdx) -
                                 mSensorMask.getTopLeftCorner(cuboidIdx)
                                ).nElements();
    }
  }

  if (params.getTimeIndex() > params.getSamplingStartTimeIndex())
  {
    // Reload temp coefficients from checkpoint file
    loadCheckpointCompressionCoefficients();

    // Reload min and max values
    for (size_t cuboidIdx = 0; cuboidIdx < mCuboidsInfo.size(); cuboidIdx++)
    {
      const std::string datasetName = (mReduceOp == ReduceOperator::kC) ? std::to_string(cuboidIdx + 1) + "_c" : std::to_string(cuboidIdx + 1);
      //Logger::log(Logger::LogLevel::kBasic, datasetName + " ");
      //Logger::log(Logger::LogLevel::kBasic, std::to_string(minValue[cuboidIdx]));
      loadMinMaxValues(mFile, mGroup, datasetName, cuboidIdx);
    }
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
      DimensionSizes datasetPosition(0, 0, 0, 0); // 4D position in the dataset
      DimensionSizes cuboidSize(0, 0, 0, 0);      // Size of the cuboid

      datasetPosition.nt = mSampledTimeStep;

      const size_t slabSize = mSourceMatrix.getDimensionSizes().ny * mSourceMatrix.getDimensionSizes().nx;
      const size_t rowSize  = mSourceMatrix.getDimensionSizes().nx;
      const float* sourceData = mSourceMatrix.getData();

      // iterate over all cuboid to be sampled
      for (size_t cuboidIdx = 0; cuboidIdx < mCuboidsInfo.size(); cuboidIdx++)
      {
        const DimensionSizes topLeftCorner     = mSensorMask.getTopLeftCorner(cuboidIdx);
        const DimensionSizes bottomRightCorner = mSensorMask.getBottomRightCorner(cuboidIdx);

        cuboidSize    = bottomRightCorner - topLeftCorner;
        cuboidSize.nt = 1;

        mFile.writeCuboidToHyperSlab(mCuboidsInfo[cuboidIdx].cuboidId,
                                     datasetPosition,
                                     topLeftCorner, // position in the SourceMatrix
                                     cuboidSize,
                                     mSourceMatrix.getDimensionSizes(),
                                     mSourceMatrix.getData());

        size_t storeBufferIndex = 0;
        for (size_t z = topLeftCorner.nz; z <= bottomRightCorner.nz; z++)
        {
          for (size_t y = topLeftCorner.ny; y <= bottomRightCorner.ny; y++)
          {
            for (size_t x = topLeftCorner.nx; x <= bottomRightCorner.nx; x++)
            {
              const size_t sourceIndex = z * slabSize + y * rowSize + x;
              checkOrSetMinMaxValue(minValue[cuboidIdx], maxValue[cuboidIdx], sourceData[sourceIndex], minValueIndex[cuboidIdx], maxValueIndex[cuboidIdx], cuboidSize.nElements() * mSampledTimeStep + storeBufferIndex);
              storeBufferIndex++;
            }
          }
        }

      }
      mSampledTimeStep++;   // Move forward in time

      break;

     // At the time being, this version using manual data lining up seems to be slower, and is not used.
     // sampleAggregated<ReduceOperator::kNone>();
    }// case kNone
    case ReduceOperator::kC:
    {
      sampleAggregated<ReduceOperator::kC>();
      break;
    }
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
  if (mReduceOp != ReduceOperator::kNone && mReduceOp != ReduceOperator::kC)
  {
    flushBufferToFile();
  }

  // Store min and max values
  for (size_t cuboidIdx = 0; cuboidIdx < mCuboidsInfo.size(); cuboidIdx++)
  {
    const std::string datasetName = (mReduceOp == ReduceOperator::kC) ? std::to_string(cuboidIdx + 1) + "_c" : std::to_string(cuboidIdx + 1);
    //Logger::log(Logger::LogLevel::kBasic, datasetName + " ");
    //Logger::log(Logger::LogLevel::kBasic, std::to_string(minValue[cuboidIdx]));
    storeMinMaxValues(mFile, mGroup, datasetName, cuboidIdx);
  }
}// end of postProcess
//----------------------------------------------------------------------------------------------------------------------

/**
 * Checkpoint the stream.
 */
void CuboidOutputStream::checkpoint()
{
  // raw data has already been flushed, others has to be flushed here
  if (mReduceOp != ReduceOperator::kNone && mReduceOp != ReduceOperator::kC)
  {
    flushBufferToFile();
  }
  storeCheckpointCompressionCoefficients();

  // Store min and max values
  for (size_t cuboidIdx = 0; cuboidIdx < mCuboidsInfo.size(); cuboidIdx++)
  {
    const std::string datasetName = (mReduceOp == ReduceOperator::kC) ? std::to_string(cuboidIdx + 1) + "_c" : std::to_string(cuboidIdx + 1);
    //Logger::log(Logger::LogLevel::kBasic, datasetName + " ");
    //Logger::log(Logger::LogLevel::kBasic, std::to_string(minValue[cuboidIdx]));
    storeMinMaxValues(mFile, mGroup, datasetName, cuboidIdx);
  }
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
  const Parameters& params = Parameters::getInstance();

  // if time series then Number of steps else 1
  size_t nSampledTimeSteps = (mReduceOp == ReduceOperator::kNone)
                             ? params.getNt() - params.getSamplingStartTimeIndex()
                             : 0; // will be a 3D dataset

  // Set cuboid dimensions (subtract two corners (add 1) and use the appropriate component)
  DimensionSizes cuboidSize((mSensorMask.getBottomRightCorner(cuboidIdx) - mSensorMask.getTopLeftCorner(cuboidIdx)).nx,
                            (mSensorMask.getBottomRightCorner(cuboidIdx) - mSensorMask.getTopLeftCorner(cuboidIdx)).ny,
                            (mSensorMask.getBottomRightCorner(cuboidIdx) - mSensorMask.getTopLeftCorner(cuboidIdx)).nz,
                             nSampledTimeSteps);

  if (mReduceOp == ReduceOperator::kC)
  {
    size_t steps = params.getNt() - params.getSamplingStartTimeIndex();
    nSampledTimeSteps = size_t(std::max(float(floor(float(steps) / mCompressHelper->getOSize())) - 1, 1.0f));

    cuboidSize.nt = nSampledTimeSteps;
    cuboidSize.nx *= mCompressHelper->getHarmonics() * 2;
  }

  // Set chunk size
  // If the size of the cuboid is bigger than 32 MB per timestep, set the chunk to approx 4MB
  size_t nSlabs = 1; //at least one slab
  DimensionSizes cuboidChunkSize(cuboidSize.nx,
                                 cuboidSize.ny,
                                 cuboidSize.nz,
                                 (mReduceOp == ReduceOperator::kNone || mReduceOp == ReduceOperator::kC) ? 1 : 0);

  if (cuboidChunkSize.nElements() > (kChunkSize4MB * 8))
  {
    while (nSlabs * cuboidSize.nx * cuboidSize.ny < kChunkSize4MB) nSlabs++;
    cuboidChunkSize.nz = nSlabs;
  }

  // Indexed from 1
  const std::string datasetName = (mReduceOp == ReduceOperator::kC) ? std::to_string(cuboidIdx + 1) + "_c" : std::to_string(cuboidIdx + 1);

  hid_t dataset = mFile.createDataset(mGroup,
                                     datasetName,
                                     cuboidSize,
                                     cuboidChunkSize,
                                     Hdf5File::MatrixDataType::kFloat,
                                     params.getCompressionLevel());

  // Write dataset parameters
  mFile.writeMatrixDomainType(mGroup, datasetName, Hdf5File::MatrixDomainType::kReal);
  mFile.writeMatrixDataType  (mGroup, datasetName, Hdf5File::MatrixDataType::kFloat);

  // Write compression parameters as attributes
  if (mReduceOp == ReduceOperator::kC)
  {
    mFile.writeLongLongAttribute(mGroup, datasetName, "c_harmonics", ssize_t(mCompressHelper->getHarmonics()));
    mFile.writeStringAttribute(mGroup, datasetName, "c_type", "c");
    mFile.writeFloatAttribute(mGroup, datasetName, "c_period", mCompressHelper->getPeriod());
    mFile.writeLongLongAttribute(mGroup, datasetName, "c_mos", ssize_t(mCompressHelper->getMos()));
    mFile.writeStringAttribute(mGroup, datasetName, "src_dataset_name", "/" + mRootObjectName + "/" + std::to_string(cuboidIdx + 1));
  }

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

  DimensionSizes cuboidSize(0, 0, 0, 0);      // Size of the cuboid

  const float* sourceData = mSourceMatrix.getData();

  size_t cuboidInBufferStart = 0;

  if (reduceOp == ReduceOperator::kC)
  {
    // Compression
    // Compute local index and flags
    mStepLocal = mSampledTimeStep % (mCompressHelper->getBSize() - 1);
    mSavingFlag = ((mStepLocal + 1) % mCompressHelper->getOSize() == 0) ? true : false;
    mOddFrameFlag = ((mCompressedTimeStep + 1) % 2 == 0) ? true : false;
  }

  // Parallelise within the cuboid - Since a typical number of cuboids is 1, than we have to paralelise inside
  #pragma omp parallel
  for (size_t cuboidIdx = 0; cuboidIdx < mSensorMask.getDimensionSizes().ny; cuboidIdx++)
  {
    const DimensionSizes topLeftCorner     = mSensorMask.getTopLeftCorner(cuboidIdx);
    const DimensionSizes bottomRightCorner = mSensorMask.getBottomRightCorner(cuboidIdx);

    size_t cuboidSlabSize = (bottomRightCorner.ny - topLeftCorner.ny + 1) *
                            ((reduceOp == ReduceOperator::kC) ? ((bottomRightCorner.nx - topLeftCorner.nx + 1) * mCompressHelper->getHarmonics()) : (bottomRightCorner.nx - topLeftCorner.nx + 1));
    size_t cuboidRowSize  = ((reduceOp == ReduceOperator::kC) ? ((bottomRightCorner.nx - topLeftCorner.nx + 1) * mCompressHelper->getHarmonics()) : (bottomRightCorner.nx - topLeftCorner.nx + 1));

    cuboidSize    = bottomRightCorner - topLeftCorner;
    cuboidSize.nt = 1;

    #pragma omp for collapse(3)
    for (size_t z = topLeftCorner.nz; z <= bottomRightCorner.nz; z++)
      for (size_t y = topLeftCorner.ny; y <= bottomRightCorner.ny; y++)
        for (size_t x = topLeftCorner.nx; x <= bottomRightCorner.nx; x++)
        {
          const size_t storeBufferIndex = cuboidInBufferStart +
                                          (z - topLeftCorner.nz) * cuboidSlabSize +
                                          (y - topLeftCorner.ny) * cuboidRowSize  +
                                          ((reduceOp == ReduceOperator::kC) ? ((x - topLeftCorner.nx) * mCompressHelper->getHarmonics()) : (x - topLeftCorner.nx));

          const size_t sourceIndex = z * slabSize + y * rowSize + x;

          // based on template parameter
          switch (reduceOp)
          {
            case ReduceOperator::kNone:
            {
              mStoreBuffer[storeBufferIndex] = sourceData[sourceIndex];
              checkOrSetMinMaxValue(minValue[cuboidIdx], maxValue[cuboidIdx], mStoreBuffer[storeBufferIndex], minValueIndex[cuboidIdx], maxValueIndex[cuboidIdx], cuboidSize.nElements() * mSampledTimeStep + storeBufferIndex);

              break;
            }
            case ReduceOperator::kC:
            {
              checkOrSetMinMaxValue(minValue[cuboidIdx], maxValue[cuboidIdx], sourceData[sourceIndex], minValueIndex[cuboidIdx], maxValueIndex[cuboidIdx], cuboidSize.nElements() * mSampledTimeStep + (storeBufferIndex / mCompressHelper->getHarmonics()));

              //For every harmonics
              for (size_t ih = 0; ih < mCompressHelper->getHarmonics(); ih++)
              {
                size_t pH = storeBufferIndex + ih;
                size_t bIndex = ih * mCompressHelper->getBSize() + mStepLocal;

                // Correlation step
                reinterpret_cast<FloatComplex*>(mStoreBuffer)[pH] += mCompressHelper->getBE()[bIndex] * sourceData[sourceIndex];
                reinterpret_cast<FloatComplex*>(mStoreBuffer2)[pH] += mCompressHelper->getBE_1()[bIndex] * sourceData[sourceIndex];
              }
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
       if (mReduceOp == ReduceOperator::kC)
         cuboidInBufferStart += (bottomRightCorner - topLeftCorner).nElements() * mCompressHelper->getHarmonics();
       else
         cuboidInBufferStart += (bottomRightCorner - topLeftCorner).nElements();
     }
  }

  if (reduceOp == ReduceOperator::kC)
  {
    if (mSavingFlag)
    {
      // Select accumulated value
      float* data = mOddFrameFlag ? mStoreBuffer : mStoreBuffer2;

      // Store selected buffer
      if (mCompressedTimeStep > 0)
      {
        flushBufferToFile(data);
      }

      // Set zeros for next accumulation
      //memset(data, 0, mBufferSize * sizeof(float));
      #pragma omp parallel for simd schedule(static)
      for (size_t i = 0; i < mBufferSize; i++)
      {
        data[i] = 0.0f;
      }
      mCompressedTimeStep++;
    }
    mSampledTimeStep++;
  }

}// end of sampleAggregated
//----------------------------------------------------------------------------------------------------------------------

/**
 * Flush the buffer to the file (to multiple datasets if necessary).
 */
void CuboidOutputStream::flushBufferToFile(float* bufferToFlush)
{
  DimensionSizes position (0, 0, 0, 0);
  DimensionSizes blockSize(0, 0, 0, 0);

  if (mReduceOp == ReduceOperator::kNone) position.nt = mSampledTimeStep;
  if (mReduceOp == ReduceOperator::kC) position.nt = mCompressedTimeStep - 1;

  for (size_t cuboidIdx = 0; cuboidIdx < mCuboidsInfo.size(); cuboidIdx++)
  {
    blockSize = mSensorMask.getBottomRightCorner(cuboidIdx) - mSensorMask.getTopLeftCorner(cuboidIdx);
    if (mReduceOp == ReduceOperator::kC) blockSize.nx *= mCompressHelper->getHarmonics() * 2;
    blockSize.nt = 1;

    mFile.writeHyperSlab(mCuboidsInfo[cuboidIdx].cuboidId,
                         position,
                         blockSize,
                         ((bufferToFlush != nullptr) ? bufferToFlush : mStoreBuffer) + mCuboidsInfo[cuboidIdx].startingPossitionInBuffer);
  }

  if (mReduceOp != ReduceOperator::kC) mSampledTimeStep++;
}// end of flushBufferToFile
//----------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Private methods --------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//
