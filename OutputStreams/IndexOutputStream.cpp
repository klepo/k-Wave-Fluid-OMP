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
IndexOutputStream::IndexOutputStream(Hdf5File& file,
                                     MatrixName& datasetName,
                                     const RealMatrix& sourceMatrix,
                                     const IndexMatrix& sensorMask,
                                     const ReduceOperator reduceOp,
                                     float* bufferToReuse,
                                     OutputStreamContainer* outputStreamContainer,
                                     bool doNotSaveFlag)
  : BaseOutputStream(file, datasetName, sourceMatrix, reduceOp, bufferToReuse, outputStreamContainer, doNotSaveFlag),
    mSensorMask(sensorMask),
    mDataset(H5I_BADID),
    mSampledTimeStep(0) {
  mMinValue.value = std::numeric_limits<float>::max();
  mMaxValue.value = std::numeric_limits<float>::min();
  mMinValue.index = 0;
  mMaxValue.index = 0;
} // end of IndexOutputStream
//----------------------------------------------------------------------------------------------------------------------

/**
 * Destructor.
 */
IndexOutputStream::~IndexOutputStream() {
  close();
  // free memory only if it was allocated
  if (!mBufferReuse)
    freeMemory();
} // end of Destructor
//----------------------------------------------------------------------------------------------------------------------

/**
 * Create a HDF5 stream, create a dataset, and allocate data for it.
 */
void IndexOutputStream::create() {
  // Set buffer size
  // Extend "x" dimension for compression coefficients
  mBufferSize = (mReduceOp == ReduceOperator::kC) ? ceilf(mSensorMask.size() * mComplexSize) * mCompressHelper->getHarmonics() : mSensorMask.size();

  // Don't create dataset for compression coefficients if only kIAvgC or kQTermC should be stored
  if (!mDoNotSaveFlag) {
    const Parameters& params = Parameters::getInstance();

    // Derive dataset dimension sizes
    DimensionSizes datasetSize(mBufferSize, 1, 1);

    if (mReduceOp == ReduceOperator::kNone) {
      datasetSize = DimensionSizes(mBufferSize,
                                   params.getNt() - params.getSamplingStartTimeIndex(),
                                   1);
    } else if (mReduceOp == ReduceOperator::kC) {
      // NOTE minimal useful number of steps for compression is 1 period.
      size_t steps = params.getNt() - params.getSamplingStartTimeIndex();
      size_t compressedSteps = size_t(std::max(float(floor(float(steps) / mCompressHelper->getOSize())), 1.0f));
      datasetSize = DimensionSizes(mBufferSize, compressedSteps, 1);
    }

    // Set HDF5 chunk size
    DimensionSizes chunkSize(mBufferSize, 1, 1);
    // for data bigger than 32 MB
    if (mBufferSize > (kChunkSize4MB * 8)) {
      chunkSize.nx = kChunkSize4MB; // set chunk size to MB
    }

    if (mFile.datasetExists(mFile.getRootGroup(), mRootObjectName)) {
      mDataset = mFile.openDataset(mFile.getRootGroup(), mRootObjectName);
    } else {
      // Create a dataset under the root group
      mDataset = mFile.createDataset(mFile.getRootGroup(),
                                     mRootObjectName,
                                     datasetSize,
                                     chunkSize,
                                     Hdf5File::MatrixDataType::kFloat,
                                     params.getCompressionLevel());

      // Write dataset parameters
      mFile.writeMatrixDomainType(mFile.getRootGroup(), mRootObjectName, Hdf5File::MatrixDomainType::kReal);
      mFile.writeMatrixDataType(mFile.getRootGroup(), mRootObjectName, Hdf5File::MatrixDataType::kFloat);
    }

    // Write compression parameters as attributes
    if (mReduceOp == ReduceOperator::kC) {
      mFile.writeLongLongAttribute(mFile.getRootGroup(), mRootObjectName, "c_harmonics", ssize_t(mCompressHelper->getHarmonics()));
      mFile.writeStringAttribute(mFile.getRootGroup(), mRootObjectName, "c_type", "c");
      mFile.writeFloatAttribute(mFile.getRootGroup(), mRootObjectName, "c_period", mCompressHelper->getPeriod());
      mFile.writeLongLongAttribute(mFile.getRootGroup(), mRootObjectName, "c_mos", ssize_t(mCompressHelper->getMos()));
      mFile.writeLongLongAttribute(mFile.getRootGroup(), mRootObjectName, "c_shift", ssize_t(mShiftFlag));
      mFile.writeFloatAttribute(mFile.getRootGroup(), mRootObjectName, "c_complex_size", mComplexSize);
      mFile.writeLongLongAttribute(mFile.getRootGroup(), mRootObjectName, "c_max_exp", mE);
    }
  }

  // Sampled time step
  mSampledTimeStep = 0;

  // Allocate memory if needed
  if (!mBufferReuse)
    allocateMemory();

  mCurrentStoreBuffer = mStoreBuffer;
  if (mReduceOp == ReduceOperator::kC) {
    mCurrentStoreBuffer = nullptr;
  }
} // end of create
//----------------------------------------------------------------------------------------------------------------------

/**
 * Reopen the output stream after restart.
 */
void IndexOutputStream::reopen() {
  // Get parameters
  const Parameters& params = Parameters::getInstance();

  // Set buffer size
  mBufferSize = mSensorMask.size();
  std::string objectName = mRootObjectName;

  if (mReduceOp == ReduceOperator::kC) {
    mBufferSize = size_t(ceilf(mSensorMask.size() * mComplexSize)) * mCompressHelper->getHarmonics();
  }

  // Allocate memory if needed
  if (!mBufferReuse)
    allocateMemory();
  mCurrentStoreBuffer = mStoreBuffer;

  // Reopen the dataset
  if (!mDoNotSaveFlag) {
    mDataset = mFile.openDataset(mFile.getRootGroup(), objectName);
  }

  if (mReduceOp == ReduceOperator::kNone || mReduceOp == ReduceOperator::kC || mReduceOp == ReduceOperator::kIAvgC) { // raw time series - just seek to the right place in the dataset
    mSampledTimeStep = (params.getTimeIndex() < params.getSamplingStartTimeIndex()) ? 0 : (params.getTimeIndex() - params.getSamplingStartTimeIndex());
    if (mReduceOp == ReduceOperator::kC || mReduceOp == ReduceOperator::kIAvgC) {
      mCompressedTimeStep = size_t(std::max(float(floor(float(mSampledTimeStep) / mCompressHelper->getOSize())), 0.0f));
    }
  } else if (mReduceOp != ReduceOperator::kIAvg && mReduceOp != ReduceOperator::kQTerm && mReduceOp != ReduceOperator::kQTermC && !mDoNotSaveFlag) { // aggregated quantities - reload data
    mSampledTimeStep = 0;
    // read only if it is necessary (it is anything to read).
    if (params.getTimeIndex() > params.getSamplingStartTimeIndex()) {
      // Since there is only a single timestep in the dataset, I can read the whole dataset
      mFile.readCompleteDataset(mFile.getRootGroup(),
                                objectName,
                                DimensionSizes(mBufferSize, 1, 1),
                                mStoreBuffer);
    }
  }

  if (params.getTimeIndex() > params.getSamplingStartTimeIndex()) {
    // Reload temp coefficients from checkpoint file
    loadCheckpointCompressionCoefficients();

    if (!mDoNotSaveFlag) {
      // Reload min and max values
      loadMinMaxValues(mFile, mFile.getRootGroup(), objectName, mMinValue, mMaxValue);
    }
  }
} // end of reopen
//----------------------------------------------------------------------------------------------------------------------

/**
 * Sample grid points, line them up in the buffer an flush to the disk unless a reduction operator is applied.
 */
void IndexOutputStream::sample() {
  const float* sourceData = mSourceMatrix.getData();
  const size_t* sensorData = mSensorMask.getData();

  switch (mReduceOp) {
    case ReduceOperator::kNone: {
#pragma omp parallel
      {
        // ReducedValue minValueLocal = mMinValue;
        // ReducedValue maxValueLocal = mMaxValue;

        // For every point
#pragma omp for nowait
        for (size_t i = 0; i < mBufferSize; i++) {
          mStoreBuffer[i] = sourceData[sensorData[i]];
          // checkOrSetMinMaxValue(minValueLocal, maxValueLocal, mStoreBuffer[i], mBufferSize * mSampledTimeStep + i);
        }
        //checkOrSetMinMaxValueGlobal(mMinValue, mMaxValue, minValueLocal, maxValueLocal);
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
    } // case kNone

    case ReduceOperator::kC: {
      // Compression
      // Compute local index and flags
      mStepLocal = mSampledTimeStep % (mCompressHelper->getBSize() - 1);
      mSavingFlag = ((mStepLocal + 1) % mCompressHelper->getOSize() == 0) ? true : false;
      mOddFrameFlag = ((mCompressedTimeStep + 1) % 2 == 0) ? true : false;
      const bool noCompressionOverlapFlag = Parameters::getInstance().getNoCompressionOverlapFlag();
      // noCompressionOverlapFlag -> mStoreBuffer == mStoreBuffer2
      mMirrorFirstHalfFrameFlag = (mCompressedTimeStep == 0 && mSavingFlag && !noCompressionOverlapFlag) ? true : false;

      FloatComplex* mStoreBufferFloatC = reinterpret_cast<FloatComplex*>(mStoreBuffer);
      FloatComplex* mStoreBuffer2FloatC = reinterpret_cast<FloatComplex*>(mStoreBuffer2);
      uint8_t* mStoreBufferInt8 = reinterpret_cast<uint8_t*>(mStoreBuffer);
      uint8_t* mStoreBuffer2Int8 = reinterpret_cast<uint8_t*>(mStoreBuffer2);

#pragma omp parallel
      {
        //ReducedValue minValueLocal = mMinValue;
        //ReducedValue maxValueLocal = mMaxValue;

        // For every point
#pragma omp for nowait
        for (size_t i = 0; i < mSensorMask.size(); i++) {
          //checkOrSetMinMaxValue(minValueLocal, maxValueLocal, sourceData[sensorData[i]], sensorMask.size() * mSampledTimeStep + i);
          const size_t storeBufferIndexC = mCompressHelper->getHarmonics() * i;
          const size_t sourceIndex = sensorData[i];

          //For every harmonics
          for (size_t ih = 0; ih < mCompressHelper->getHarmonics(); ih++) {
            size_t pH = storeBufferIndexC + ih;
            const size_t bIndex = ih * mCompressHelper->getBSize() + mStepLocal;

            // 40-bit complex float compression
            if (Parameters::getInstance().get40bitCompressionFlag()) {
              pH = pH * 5;
              FloatComplex cc1;
              FloatComplex cc2;
              if (Parameters::getInstance().getNoCompressionOverlapFlag()) {
                CompressHelper::convert40bToFloatC(&mStoreBufferInt8[pH], cc1, mE);
                cc1 += mBE[bIndex] * sourceData[sourceIndex] + mBE_1[bIndex] * sourceData[sourceIndex];
                CompressHelper::convertFloatCTo40b(cc1, &mStoreBufferInt8[pH], mE);
              } else {
                CompressHelper::convert40bToFloatC(&mStoreBufferInt8[pH], cc1, mE);
                CompressHelper::convert40bToFloatC(&mStoreBuffer2Int8[pH], cc2, mE);
                cc1 += mBE[bIndex] * sourceData[sourceIndex];
                cc2 += mBE_1[bIndex] * sourceData[sourceIndex];
                CompressHelper::convertFloatCTo40b(cc1, &mStoreBufferInt8[pH], mE);
                CompressHelper::convertFloatCTo40b(cc2, &mStoreBuffer2Int8[pH], mE);
                // Mirror first "half" frame
                if (mMirrorFirstHalfFrameFlag) {
                  cc2 += cc1;
                  CompressHelper::convertFloatCTo40b(cc2, &mStoreBuffer2Int8[pH], mE);
                }
              }
            } else {
              // Correlation step
              mStoreBufferFloatC[pH] += mBE[bIndex] * sourceData[sourceIndex];
              mStoreBuffer2FloatC[pH] += mBE_1[bIndex] * sourceData[sourceIndex];
              // Mirror first "half" frame
              if (mMirrorFirstHalfFrameFlag) {
                mStoreBuffer2FloatC[pH] += mStoreBufferFloatC[pH];
              }
            }
          }
        }
        //checkOrSetMinMaxValueGlobal(mMinValue, mMaxValue, minValueLocal, maxValueLocal);
      }

      const size_t steps = Parameters::getInstance().getNt() - Parameters::getInstance().getSamplingStartTimeIndex();
      const bool lastStep = ((steps - mSampledTimeStep == 1) && steps <= mCompressHelper->getOSize()) ? true : false;
      if (mSavingFlag || lastStep) {
        // Select accumulated value (mStoreBuffer2 is first)
        mCurrentStoreBuffer = mOddFrameFlag ? mStoreBuffer : mStoreBuffer2;

        // Store selected buffer
        if (!mDoNotSaveFlag) {
          flushBufferToFile(mCurrentStoreBuffer);
        }

        mCompressedTimeStep++;
      }
      mSampledTimeStep++;
      break;
    } // case kC

    case ReduceOperator::kRms: {
#pragma omp parallel for
      for (size_t i = 0; i < mBufferSize; i++) {
        mStoreBuffer[i] += (sourceData[sensorData[i]] * sourceData[sensorData[i]]);
      }
      break;
    } // case kRms

    case ReduceOperator::kMax: {
#pragma omp parallel for
      for (size_t i = 0; i < mBufferSize; i++) {
        mStoreBuffer[i] = std::max(mStoreBuffer[i], sourceData[sensorData[i]]);
      }
      break;
    } // case kMax

    case ReduceOperator::kMin: {
#pragma omp parallel for
      for (size_t i = 0; i < mBufferSize; i++) {
        mStoreBuffer[i] = std::min(mStoreBuffer[i], sourceData[sensorData[i]]);
      }
      break;
    } //case kMin

    default: {
      break;
    }
  } // switch
} // end of sample
//----------------------------------------------------------------------------------------------------------------------

/**
 * Post sampling step, can work with other filled stream buffers
 */
void IndexOutputStream::postSample() {
  if (mReduceOp == ReduceOperator::kIAvgC && !Parameters::getInstance().getOnlyPostProcessingFlag()) {
    float* bufferP = (*mOutputStreamContainer)[OutputStreamContainer::OutputStreamIdx::kPressureC].getCurrentStoreBuffer();
    float* bufferU = (*mOutputStreamContainer)[static_cast<OutputStreamContainer::OutputStreamIdx>(mVelocityOutputStreamIdx)].getCurrentStoreBuffer();

    uint8_t* mBufferPInt8 = reinterpret_cast<uint8_t*>(bufferP);
    uint8_t* mBufferUInt8 = reinterpret_cast<uint8_t*>(bufferU);

    // TODO check the length of bufferP == the length of bufferU
    if (bufferP && bufferU) {
#pragma omp parallel for
      for (size_t i = 0; i < mBufferSize; i++) {
        size_t offset = mCompressHelper->getHarmonics() * i;
        //For every harmonics
        for (size_t ih = 0; ih < mCompressHelper->getHarmonics(); ih++) {
          size_t pH = offset + ih;
          FloatComplex sCP;
          FloatComplex sCU;
          if (Parameters::getInstance().get40bitCompressionFlag()) {
            pH = pH * 5;
            CompressHelper::convert40bToFloatC(&mBufferPInt8[pH], sCP, CompressHelper::kMaxExpP);
            CompressHelper::convert40bToFloatC(&mBufferUInt8[pH], sCU, CompressHelper::kMaxExpU);
          } else {
            sCP = reinterpret_cast<FloatComplex*>(bufferP)[pH];
            sCU = reinterpret_cast<FloatComplex*>(bufferU)[pH];
          }
          mStoreBuffer[i] += real(sCP * conj(sCU)) / 2.0f;
        }
      }
      mCompressedTimeStep++;
    }
  }
} // end of postSample
//----------------------------------------------------------------------------------------------------------------------

/**
 * Apply post-processing on the buffer and flush it to the file.
 */
void IndexOutputStream::postProcess() {
  // run inherited method
  BaseOutputStream::postProcess();

  if (mReduceOp == ReduceOperator::kIAvgC &&
      !Parameters::getInstance().getOnlyPostProcessingFlag()) {
#pragma omp parallel for simd
    for (size_t i = 0; i < mBufferSize; i++) {
      mStoreBuffer[i] = mStoreBuffer[i] / (mCompressedTimeStep);
    }
    mCompressedTimeStep = 0;
  }

  // When no reduction operator is applied, the data is flushed after every time step
  if (mReduceOp != ReduceOperator::kNone &&
      mReduceOp != ReduceOperator::kC &&
      mReduceOp != ReduceOperator::kQTerm &&
      mReduceOp != ReduceOperator::kQTermC &&
      !mDoNotSaveFlag) {
    flushBufferToFile();
  }

  /*if (!mDoNotSaveFlag && !Parameters::getInstance().getOnlyPostProcessingFlag())
  {
    // Store min and max values
    storeMinMaxValues(mFile, mFile.getRootGroup(), mRootObjectName, mMinValue, mMaxValue);
  }*/
} // end of postProcessing
//----------------------------------------------------------------------------------------------------------------------

/**
 * Apply post-processing 2 on the buffer and flush it to the file.
 */
void IndexOutputStream::postProcess2() {
  // run inherited method
  BaseOutputStream::postProcess2();

  if (mReduceOp == ReduceOperator::kQTerm ||
      mReduceOp == ReduceOperator::kQTermC) {
    flushBufferToFile();
  }

} // end of postProcess2
//----------------------------------------------------------------------------------------------------------------------

/**
 * Checkpoint the stream and close.
 */
void IndexOutputStream::checkpoint() {
  // raw data has already been flushed, others has to be flushed here
  if (mReduceOp != ReduceOperator::kNone &&
      mReduceOp != ReduceOperator::kC &&
      mReduceOp != ReduceOperator::kQTerm &&
      mReduceOp != ReduceOperator::kQTermC &&
      !mDoNotSaveFlag) {
    flushBufferToFile();
  }
  storeCheckpointCompressionCoefficients();

  /*if (!mDoNotSaveFlag)
  {
    // Store min and max values
    storeMinMaxValues(mFile, mFile.getRootGroup(), mRootObjectName, mMinValue, mMaxValue);
  }*/
} // end of checkpoint
//----------------------------------------------------------------------------------------------------------------------

/**
 * Close stream (apply post-processing if necessary, flush data and close).
 */
void IndexOutputStream::close() {
  if (!mDoNotSaveFlag) { // the dataset is still opened
    if (mDataset != H5I_BADID) {
      mFile.closeDataset(mDataset);
    }
  }
  mDataset = H5I_BADID;
} // end of close
//----------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------ Protected methods -------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Flush the buffer down to the file at the actual position.
 */
void IndexOutputStream::flushBufferToFile(float* bufferToFlush) {
  mFile.writeHyperSlab(mDataset,
                       DimensionSizes(0, (mReduceOp == ReduceOperator::kC) ? mCompressedTimeStep : mSampledTimeStep, 0),
                       DimensionSizes(mBufferSize, 1, 1),
                       (bufferToFlush != nullptr) ? bufferToFlush : mStoreBuffer);

  if (mReduceOp != ReduceOperator::kC)
    mSampledTimeStep++;
} // end of flushToFile
//----------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Private methods --------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//
