/**
 * @file      IndexOutputStream.cpp
 *
 * @author    Jiri Jaros \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            jarosjir@fit.vutbr.cz
 *
 * @brief     The implementation file of the class saving data based on index senor mask into
 *            the output HDF5 file.
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

#include <OutputStreams/IndexOutputStream.h>
#include <Parameters/Parameters.h>
#include <Containers/OutputStreamContainer.h>

//--------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------- Constants -----------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//


//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Public methods ---------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Constructor - links the HDF5 dataset, SourceMatrix, and SensorMask together.
 */
IndexOutputStream::IndexOutputStream(Hdf5File&              file,
                                     MatrixName&            datasetName,
                                     const RealMatrix&      sourceMatrix,
                                     const IndexMatrix&     sensorMask,
                                     const ReduceOperator   reduceOp,
                                     float*                 bufferToReuse,
                                     OutputStreamContainer* outputStreamContainer,
                                     bool                   doNotSaveFlag)
  : BaseOutputStream(file, datasetName, sourceMatrix, reduceOp, bufferToReuse, outputStreamContainer, doNotSaveFlag),
    sensorMask(sensorMask),
    mDataset(H5I_BADID),
    mSampledTimeStep(0)
{
  mMinValue.value = std::numeric_limits<float>::max();
  mMaxValue.value = std::numeric_limits<float>::min();
  mMinValue.index = 0;
  mMaxValue.index = 0;
}// end of IndexOutputStream
//----------------------------------------------------------------------------------------------------------------------

/**
 * Destructor.
 */
IndexOutputStream::~IndexOutputStream()
{
  close();
  // free memory only if it was allocated
  if (!mBufferReuse) freeMemory();
}// end of Destructor
//----------------------------------------------------------------------------------------------------------------------

/**
 * Create a HDF5 stream, create a dataset, and allocate data for it.
 */
void IndexOutputStream::create()
{
  // Set buffer size
  // Extend "x" dimension for compression coefficients
  mBufferSize = (mReduceOp == ReduceOperator::kC) ? sensorMask.size() * mCompressHelper->getHarmonics() * 2 : sensorMask.size();

  // Dont create dataset for compression coefficients if only kIAvgC or kQTermC should be stored
  if (!mDoNotSaveFlag)
  {
    const Parameters& params = Parameters::getInstance();

    // Derive dataset dimension sizes
    DimensionSizes datasetSize(mBufferSize, 1, 1);

    if (mReduceOp == ReduceOperator::kNone)
    {
      datasetSize = DimensionSizes(mBufferSize,
                                   params.getNt() - params.getSamplingStartTimeIndex(),
                                   1);
    }
    else if (mReduceOp == ReduceOperator::kC)
    {
      // TODO minimal number of steps for compression (1 period)
      size_t steps = params.getNt() - params.getSamplingStartTimeIndex();
      size_t compressedSteps = size_t(std::max(float(floor(float(steps) / mCompressHelper->getOSize())), 1.0f));
      datasetSize = DimensionSizes(mBufferSize, compressedSteps, 1);
    }

    // Set HDF5 chunk size
    DimensionSizes chunkSize(mBufferSize, 1, 1);
    // for data bigger than 32 MB
    if (mBufferSize > (kChunkSize4MB * 8))
    {
      chunkSize.nx = kChunkSize4MB; // set chunk size to MB
    }

    const std::string objectName = (mReduceOp == ReduceOperator::kC) ? mRootObjectName + kCompressSuffix : mRootObjectName;

    if (mFile.datasetExists(mFile.getRootGroup(), objectName))
    {
      mDataset = mFile.openDataset(mFile.getRootGroup(), objectName);
    }
    else
    {
      // Create a dataset under the root group
      mDataset = mFile.createDataset(mFile.getRootGroup(),
                                     objectName,
                                     datasetSize,
                                     chunkSize,
                                     Hdf5File::MatrixDataType::kFloat,
                                     params.getCompressionLevel());

      // Write dataset parameters
      mFile.writeMatrixDomainType(mFile.getRootGroup(), objectName, Hdf5File::MatrixDomainType::kReal);
      mFile.writeMatrixDataType  (mFile.getRootGroup(), objectName, Hdf5File::MatrixDataType::kFloat);
    }

    // Write compression parameters as attributes
    if (mReduceOp == ReduceOperator::kC)
    {
      mFile.writeLongLongAttribute(mFile.getRootGroup(), objectName, "c_harmonics", ssize_t(mCompressHelper->getHarmonics()));
      mFile.writeStringAttribute  (mFile.getRootGroup(), objectName, "c_type", "c");
      mFile.writeFloatAttribute   (mFile.getRootGroup(), objectName, "c_period", ssize_t(mCompressHelper->getPeriod()));
      mFile.writeLongLongAttribute(mFile.getRootGroup(), objectName, "c_mos", ssize_t(mCompressHelper->getMos()));
      mFile.writeStringAttribute  (mFile.getRootGroup(), objectName, "src_dataset_name", mRootObjectName);
    }
  }

  // Sampled time step
  mSampledTimeStep = 0;

  // Allocate memory if needed
  if (!mBufferReuse) allocateMemory();
  mCurrentStoreBuffer = mStoreBuffer;

  if (mReduceOp == ReduceOperator::kC)
  {
    mCurrentStoreBuffer = nullptr;
  }
}// end of create
//----------------------------------------------------------------------------------------------------------------------

/**
 * Reopen the output stream after restart.
 */
void IndexOutputStream::reopen()
{
  // Get parameters
  const Parameters& params = Parameters::getInstance();

  // Set buffer size
  mBufferSize = sensorMask.size();
  std::string objectName = mRootObjectName;

  if (mReduceOp == ReduceOperator::kC)
  {
    mBufferSize *= mCompressHelper->getHarmonics() * 2;
    objectName = objectName + kCompressSuffix;
  }

  // Allocate memory if needed
  if (!mBufferReuse) allocateMemory();
  mCurrentStoreBuffer = mStoreBuffer;

  // Reopen the dataset
  if (!mDoNotSaveFlag)
  {
    mDataset = mFile.openDataset(mFile.getRootGroup(), objectName);
  }

  if (mReduceOp == ReduceOperator::kNone ||
      mReduceOp == ReduceOperator::kC ||
      mReduceOp == ReduceOperator::kIAvgC)
  { // raw time series - just seek to the right place in the dataset
    mSampledTimeStep = (params.getTimeIndex() < params.getSamplingStartTimeIndex()) ?
                          0 : (params.getTimeIndex() - params.getSamplingStartTimeIndex());
    if (mReduceOp == ReduceOperator::kC || mReduceOp == ReduceOperator::kIAvgC)
    {
      mCompressedTimeStep = size_t(std::max(float(floor(float(mSampledTimeStep) / mCompressHelper->getOSize())), 0.0f));
    }
  }
  else if (mReduceOp != ReduceOperator::kIAvg &&
           mReduceOp != ReduceOperator::kQTerm &&
           mReduceOp != ReduceOperator::kQTermC &&
           !mDoNotSaveFlag)
  { // aggregated quantities - reload data
    mSampledTimeStep = 0;
    // read only if it is necessary (it is anything to read).
    if (params.getTimeIndex() > params.getSamplingStartTimeIndex())
    {
      // Since there is only a single timestep in the dataset, I can read the whole dataset
      mFile.readCompleteDataset(mFile.getRootGroup(),
                                objectName,
                                DimensionSizes(mBufferSize, 1, 1),
                                mStoreBuffer);
    }
  }

  if (params.getTimeIndex() > params.getSamplingStartTimeIndex())
  {
    // Reload temp coefficients from checkpoint file
    loadCheckpointCompressionCoefficients();

    if (!mDoNotSaveFlag)
    {
      // Reload min and max values
      loadMinMaxValues(mFile, mFile.getRootGroup(), objectName, mMinValue, mMaxValue);
    }
  }
}// end of reopen
//----------------------------------------------------------------------------------------------------------------------

/**
 * Sample grid points, line them up in the buffer an flush to the disk unless a reduction operator is applied.
 */
void IndexOutputStream::sample()
{
  if (mReduceOp == ReduceOperator::kIAvg ||
      mReduceOp == ReduceOperator::kIAvgC ||
      mReduceOp == ReduceOperator::kQTerm ||
      mReduceOp == ReduceOperator::kQTermC)
  {
    return;
  }

  const float*  sourceData = mSourceMatrix.getData();
  const size_t* sensorData = sensorMask.getData();

  switch (mReduceOp)
  {
    case ReduceOperator::kNone :
    {
      #pragma omp parallel
      {
        ReducedValue minValueLocal = mMinValue;
        ReducedValue maxValueLocal = mMaxValue;

        // For every point
        #pragma omp for nowait
        for (size_t i = 0; i < mBufferSize; i++)
        {
          mStoreBuffer[i] = sourceData[sensorData[i]];
          checkOrSetMinMaxValue(minValueLocal, maxValueLocal, mStoreBuffer[i], mBufferSize * mSampledTimeStep + i);
        }
        checkOrSetMinMaxValueGlobal(mMinValue, mMaxValue, minValueLocal, maxValueLocal);
      }
      // only raw time series are flushed down to the disk every time step
      flushBufferToFile();
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
    }// case kNone

    case ReduceOperator::kC:
    {
      // Compression
      // Compute local index and flags
      mStepLocal = mSampledTimeStep % (mCompressHelper->getBSize() - 1);
      mSavingFlag = ((mStepLocal + 1) % mCompressHelper->getOSize() == 0) ? true : false;
      mOddFrameFlag = ((mCompressedTimeStep + 1) % 2 == 0) ? true : false;

      #pragma omp parallel
      {
        ReducedValue minValueLocal = mMinValue;
        ReducedValue maxValueLocal = mMaxValue;

        // For every point
        #pragma omp for nowait
        for (size_t i = 0; i < sensorMask.size(); i++)
        {
          checkOrSetMinMaxValue(minValueLocal, maxValueLocal, sourceData[sensorData[i]], sensorMask.size() * mSampledTimeStep + i);
          size_t offset = mCompressHelper->getHarmonics() * i;

          //For every harmonics
          for (size_t ih = 0; ih < mCompressHelper->getHarmonics(); ih++)
          {
            size_t pH = offset + ih;
            size_t bIndex = ih * mCompressHelper->getBSize() + mStepLocal;

            // Correlation step
            // Time shift of velocity
            if (mRootObjectName == kUxNonStaggeredName
                || mRootObjectName == kUyNonStaggeredName
                || mRootObjectName == kUzNonStaggeredName)
            {
              reinterpret_cast<FloatComplex*>(mStoreBuffer)[pH] += mCompressHelper->getBEShifted()[bIndex] * sourceData[sensorData[i]];
              reinterpret_cast<FloatComplex*>(mStoreBuffer2)[pH] += mCompressHelper->getBE_1Shifted()[bIndex] * sourceData[sensorData[i]];
            }
            else
            {
              reinterpret_cast<FloatComplex*>(mStoreBuffer)[pH] += mCompressHelper->getBE()[bIndex] * sourceData[sensorData[i]];
              reinterpret_cast<FloatComplex*>(mStoreBuffer2)[pH] += mCompressHelper->getBE_1()[bIndex] * sourceData[sensorData[i]];
            }

            if (mCompressedTimeStep == 0 && mSavingFlag)
            {
              reinterpret_cast<FloatComplex*>(mStoreBuffer2)[pH] += reinterpret_cast<FloatComplex*>(mStoreBuffer)[pH];
            }
          }
        }
        checkOrSetMinMaxValueGlobal(mMinValue, mMaxValue, minValueLocal, maxValueLocal);
      }

      if (mSavingFlag)
      {
        // Select accumulated value
        mCurrentStoreBuffer = mOddFrameFlag ? mStoreBuffer : mStoreBuffer2;

        // Store selected buffer
        //if (mCompressedTimeStep > 0)
        {
          if (!mDoNotSaveFlag)
          {
            flushBufferToFile(mCurrentStoreBuffer);
          }
        }
        mCompressedTimeStep++;
      }
      mSampledTimeStep++;
      break;
    }// case kC

    case ReduceOperator::kRms:
    {
      #pragma omp parallel for
      for (size_t i = 0; i < mBufferSize; i++)
      {
        mStoreBuffer[i] += (sourceData[sensorData[i]] * sourceData[sensorData[i]]);
      }
      break;
    }// case kRms

    case ReduceOperator::kMax:
    {
      #pragma omp parallel for
      for (size_t i = 0; i < mBufferSize; i++)
      {
        mStoreBuffer[i] = std::max(mStoreBuffer[i], sourceData[sensorData[i]]);
      }
      break;
    }// case kMax

    case ReduceOperator::kMin:
    {
      #pragma omp parallel for
      for (size_t i = 0; i < mBufferSize; i++)
      {
        mStoreBuffer[i] = std::min(mStoreBuffer[i], sourceData[sensorData[i]]);
      }
      break;
    } //case kMin
  }// switch
}// end of sample
//----------------------------------------------------------------------------------------------------------------------

/**
 * Post sampling step, can work with other filled stream buffers
 */
void IndexOutputStream::postSample()
{
  if (mReduceOp == ReduceOperator::kIAvgC)
  {
    FloatComplex* bufferP = reinterpret_cast<FloatComplex*>((*mOutputStreamContainer)[OutputStreamContainer::OutputStreamIdx::kPressureC].getCurrentStoreBuffer());
    FloatComplex* bufferI = nullptr;

    if (mRootObjectName == kIxAvgName + kCompressSuffix)
    {
      bufferI = reinterpret_cast<FloatComplex*>((*mOutputStreamContainer)[OutputStreamContainer::OutputStreamIdx::kVelocityXNonStaggeredC].getCurrentStoreBuffer());
    }
    else if (mRootObjectName == kIyAvgName + kCompressSuffix)
    {
      bufferI = reinterpret_cast<FloatComplex*>((*mOutputStreamContainer)[OutputStreamContainer::OutputStreamIdx::kVelocityYNonStaggeredC].getCurrentStoreBuffer());
    }
    else
    {
      bufferI = reinterpret_cast<FloatComplex*>((*mOutputStreamContainer)[OutputStreamContainer::OutputStreamIdx::kVelocityZNonStaggeredC].getCurrentStoreBuffer());
    }
    // TODO check the length of bufferP == the length of bufferI

    if (bufferP && bufferI)
    {
      #pragma omp parallel for
      for (size_t i = 0; i < mBufferSize; i++)
      {
        size_t offset = mCompressHelper->getHarmonics() * i;
        //For every harmonics
        for (size_t ih = 0; ih < mCompressHelper->getHarmonics(); ih++)
        {
          size_t pH = offset + ih;
          mStoreBuffer[i] += real(bufferP[pH] * conj(bufferI[pH])) / 2.0f;
        }
      }
      mCompressedTimeStep++;
    }
  }
}// end of postSample
//----------------------------------------------------------------------------------------------------------------------

/**
 * Apply post-processing on the buffer and flush it to the file.
 */
void IndexOutputStream::postProcess()
{
  // run inherited method
  BaseOutputStream::postProcess();

  if (mReduceOp == ReduceOperator::kIAvgC)
  {
    #pragma omp parallel for simd
    for (size_t i = 0; i < mBufferSize; i++)
    {
      mStoreBuffer[i] = mStoreBuffer[i] / (mCompressedTimeStep);
    }
    mCompressedTimeStep = 0;
  }

  // When no reduction operator is applied, the data is flushed after every time step
  if (mReduceOp != ReduceOperator::kNone &&
      mReduceOp != ReduceOperator::kC &&
      mReduceOp != ReduceOperator::kQTerm &&
      mReduceOp != ReduceOperator::kQTermC &&
      !mDoNotSaveFlag)
  {
    flushBufferToFile();
  }

  if (!mDoNotSaveFlag && !Parameters::getInstance().getOnlyPostProcessingFlag())
  {
    // Store min and max values
    const std::string datasetName = (mReduceOp == ReduceOperator::kC) ? mRootObjectName + kCompressSuffix : mRootObjectName;
    storeMinMaxValues(mFile, mFile.getRootGroup(), datasetName, mMinValue, mMaxValue);
  }
}// end of postProcessing
//----------------------------------------------------------------------------------------------------------------------

/**
 * Apply post-processing 2 on the buffer and flush it to the file.
 */
void IndexOutputStream::postProcess2()
{
  // run inherited method
  BaseOutputStream::postProcess2();

  if (mReduceOp == ReduceOperator::kQTerm || mReduceOp == ReduceOperator::kQTermC)
  {
    flushBufferToFile();
  }

}// end of postProcess2
//----------------------------------------------------------------------------------------------------------------------

/**
 * Checkpoint the stream and close.
 */
void IndexOutputStream::checkpoint()
{
  // raw data has already been flushed, others has to be flushed here
  if (mReduceOp != ReduceOperator::kNone &&
      mReduceOp != ReduceOperator::kC &&
      mReduceOp != ReduceOperator::kQTerm &&
      mReduceOp != ReduceOperator::kQTermC &&
      !mDoNotSaveFlag)
  {
    flushBufferToFile();
  }
  storeCheckpointCompressionCoefficients();

  if (!mDoNotSaveFlag)
  {
    // Store min and max values
    const std::string datasetName = (mReduceOp == ReduceOperator::kC) ? mRootObjectName + kCompressSuffix : mRootObjectName;
    storeMinMaxValues(mFile, mFile.getRootGroup(), datasetName, mMinValue, mMaxValue);
  }
}// end of checkpoint
//----------------------------------------------------------------------------------------------------------------------

/**
 * Close stream (apply post-processing if necessary, flush data and close).
 */
void IndexOutputStream::close()
{
  // the dataset is still opened
  if (mDataset != H5I_BADID)
  {
    if (!mDoNotSaveFlag)
    {
      mFile.closeDataset(mDataset);
    }
  }

  mDataset = H5I_BADID;
}// end of close
//----------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------ Protected methods -------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Flush the buffer down to the file at the actual position.
 */
void IndexOutputStream::flushBufferToFile(float* bufferToFlush)
{
  mFile.writeHyperSlab(mDataset,
                       DimensionSizes(0, (mReduceOp == ReduceOperator::kC) ? mCompressedTimeStep : mSampledTimeStep, 0),
                       DimensionSizes(mBufferSize, 1, 1),
                       (bufferToFlush != nullptr) ? bufferToFlush : mStoreBuffer);
  if (mReduceOp != ReduceOperator::kC) mSampledTimeStep++;
}// end of flushToFile
//----------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Private methods --------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

