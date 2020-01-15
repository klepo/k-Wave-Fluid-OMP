/**
 * @file      BaseOutputStream.cpp
 *
 * @author    Jiri Jaros \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            jarosjir@fit.vutbr.cz
 *
 * @brief     The implementation file of the class saving RealMatrix data into the output HDF5 file.
 *
 * @version   kspaceFirstOrder 2.17
 *
 * @date      11 July      2012, 10:30 (created) \n
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


#include <cmath>
#include <immintrin.h>
#include <limits>

#include <OutputStreams/BaseOutputStream.h>
#include <Parameters/Parameters.h>
#include <Containers/OutputStreamContainer.h>


//--------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------- Constants -----------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Public methods ---------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Constructor - there is no sensor mask by default!
 */
BaseOutputStream::BaseOutputStream(Hdf5File&              file,
                                   MatrixName&            rootObjectName,
                                   const RealMatrix&      sourceMatrix,
                                   const ReduceOperator   reduceOp,
                                   float*                 bufferToReuse,
                                   OutputStreamContainer* outputStreamContainer,
                                   bool                   doNotSaveFlag)
  : mFile(file),
    mRootObjectName(rootObjectName),
    mSourceMatrix(sourceMatrix),
    mReduceOp(reduceOp),
    mBufferReuse(bufferToReuse != nullptr),
    mBufferSize(0),
    mStoreBuffer(bufferToReuse),
    mOutputStreamContainer(outputStreamContainer),
    mDoNotSaveFlag(doNotSaveFlag)
{
  // Set compression variables
  if (mReduceOp == ReduceOperator::kC || mReduceOp == ReduceOperator::kIAvgC)
  {
    // Set compression helper
    mCompressHelper = &CompressHelper::getInstance();

    if (mRootObjectName == kUxNonStaggeredName + kCompressSuffix
        || mRootObjectName == kUyNonStaggeredName + kCompressSuffix
        || mRootObjectName == kUzNonStaggeredName + kCompressSuffix)
    {
      // Time shift of velocity
      mBE = mCompressHelper->getBEShifted();
      mBE_1 = mCompressHelper->getBE_1Shifted();
      mShiftFlag = true;
      mE = CompressHelper::kMaxExpU;
    }
    else
    {
      mBE = mCompressHelper->getBE();
      mBE_1 = mCompressHelper->getBE_1();
      mE = CompressHelper::kMaxExpP;
    }

    if (mRootObjectName == kIxAvgName + kCompressSuffix)
    {
      mVelocityOutputStreamIdx = static_cast<int>(OutputStreamContainer::OutputStreamIdx::kVelocityXNonStaggeredC);
    }
    else if (mRootObjectName == kIyAvgName + kCompressSuffix)
    {
      mVelocityOutputStreamIdx = static_cast<int>(OutputStreamContainer::OutputStreamIdx::kVelocityYNonStaggeredC);
    }
    else if (mRootObjectName == kIzAvgName + kCompressSuffix)
    {
      mVelocityOutputStreamIdx = static_cast<int>(OutputStreamContainer::OutputStreamIdx::kVelocityZNonStaggeredC);
    }

    if (Parameters::getInstance().get40bitCompressionFlag())
    {
      mComplexSize = 1.25f;
    }
  }
}

/**
 * Post sampling step, can work with other filled stream buffers.
 */
void BaseOutputStream::postSample()
{
}// end of postSample
//----------------------------------------------------------------------------------------------------------------------

/**
 * Post sampling step 2, can work with other filled stream buffers.
 */
void BaseOutputStream::postSample2()
{
  // Compression stuff
  if (mReduceOp == ReduceOperator::kC && mSavingFlag && mCurrentStoreBuffer)
  {
    // Set zeros for next accumulation
    {
      #pragma omp parallel for simd schedule(static)
      for (size_t i = 0; i < mBufferSize; i++)
      {
        mCurrentStoreBuffer[i] = 0.0f;
      }
    }
    mCurrentStoreBuffer = nullptr;
  }
}// end of postSample2
//----------------------------------------------------------------------------------------------------------------------

/**
 * Apply post-processing on the buffer. It supposes the elements are independent.
 */
void BaseOutputStream::postProcess()
{
  switch (mReduceOp)
  {
    case ReduceOperator::kNone:
    {
      // do nothing
      break;
    }
    case ReduceOperator::kC:
    {
      // do nothing
      break;
    }
    case ReduceOperator::kIAvg:
    {
      // do nothing
      break;
    }
    case ReduceOperator::kIAvgC:
    {
      // do nothing
      break;
    }
    case ReduceOperator::kQTerm:
    {
      // do nothing
      break;
    }
    case ReduceOperator::kQTermC:
    {
      // do nothing
      break;
    }
    case ReduceOperator::kRms:
    {
      const float scalingCoeff = 1.0f / (Parameters::getInstance().getNt() -
                                         Parameters::getInstance().getSamplingStartTimeIndex());

      #pragma omp parallel for simd
      for (size_t i = 0; i < mBufferSize; i++)
      {
        mStoreBuffer[i] = sqrt(mStoreBuffer[i] * scalingCoeff);
      }
      break;
    }
    case ReduceOperator::kMax:
    {
      // do nothing
      break;
    }
    case ReduceOperator::kMin:
    {
      // do nothing
      break;
    }
  }// switch

}// end of postProcess
//----------------------------------------------------------------------------------------------------------------------

/**
 * Apply post-processing 2 on the buffer.
 */
void BaseOutputStream::postProcess2()
{
}// end of postProcess2
//----------------------------------------------------------------------------------------------------------------------

/**
 * Get current store buffer.
 * @return Current store buffer.
 */
float *BaseOutputStream::getCurrentStoreBuffer()
{
  return mCurrentStoreBuffer;
}// end of getCurrentStoreBuffer
//----------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------ Protected methods -------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Allocate memory using a proper memory alignment.
 */
void BaseOutputStream::allocateMemory()
{
  if (mReduceOp == ReduceOperator::kC)
  {
    mStoreBuffer = (float*) _mm_malloc((mBufferSize) * sizeof(float), kDataAlignment);
    if (!mStoreBuffer)
    {
      throw std::bad_alloc();
    }
    if (Parameters::getInstance().getNoCompressionOverlapFlag())
    {
      mStoreBuffer2 = mStoreBuffer;
    }
    else
    {
      mStoreBuffer2 = (float*) _mm_malloc((mBufferSize) * sizeof(float), kDataAlignment);
      if (!mStoreBuffer2)
      {
        throw std::bad_alloc();
      }
    }
  }
  else
  {
    mStoreBuffer = (float*) _mm_malloc((mBufferSize) * sizeof(float), kDataAlignment);
    if (!mStoreBuffer)
    {
      throw std::bad_alloc();
    }
  }

  // we need different initialization for different reduction ops
  switch (mReduceOp)
  {
    case ReduceOperator::kNone:
    {
      // zero the matrix
      #pragma omp parallel for simd
      for (size_t i = 0; i < mBufferSize; i++)
      {
        mStoreBuffer[i] = 0.0f;
      }
      break;
    }

    case ReduceOperator::kC:
    {
      // zero the matrix
      #pragma omp parallel for simd
      for (size_t i = 0; i < mBufferSize; i++)
      {
        mStoreBuffer[i] = 0.0f;
        mStoreBuffer2[i] = 0.0f;
      }
      break;
    }

    case ReduceOperator::kIAvg:
    {
      // zero the matrix
      #pragma omp parallel for simd
      for (size_t i = 0; i < mBufferSize; i++)
      {
        mStoreBuffer[i] = 0.0f;
      }
      break;
    }

    case ReduceOperator::kIAvgC:
    {
      // zero the matrix
      #pragma omp parallel for simd
      for (size_t i = 0; i < mBufferSize; i++)
      {
        mStoreBuffer[i] = 0.0f;
      }
      break;
    }

    case ReduceOperator::kQTerm:
    {
      // zero the matrix
      #pragma omp parallel for simd
      for (size_t i = 0; i < mBufferSize; i++)
      {
        mStoreBuffer[i] = 0.0f;
      }
      break;
    }

    case ReduceOperator::kQTermC:
    {
      // zero the matrix
      #pragma omp parallel for simd
      for (size_t i = 0; i < mBufferSize; i++)
      {
        mStoreBuffer[i] = 0.0f;
      }
      break;
    }

    case ReduceOperator::kRms:
    {
      // zero the matrix
      #pragma omp parallel for simd
      for (size_t i = 0; i < mBufferSize; i++)
      {
        mStoreBuffer[i] = 0.0f;
      }
      break;
    }

    case ReduceOperator::kMax:
    {
      // set the values to the highest negative float value
      #pragma omp parallel for simd
      for (size_t i = 0; i < mBufferSize; i++)
      {
        mStoreBuffer[i] = -1.0f * std::numeric_limits<float>::max();
      }
      break;
    }

    case ReduceOperator::kMin:
    {
      // set the values to the highest float value
      #pragma omp parallel for simd
      for (size_t i = 0; i < mBufferSize; i++)
      {
        mStoreBuffer[i] = std::numeric_limits<float>::max();
      }
      break;
    }
  }// switch
}// end of allocateMemory
//----------------------------------------------------------------------------------------------------------------------

/**
 * Free memory.
 */
void BaseOutputStream::freeMemory()
{
  if (mStoreBuffer)
  {
    _mm_free(mStoreBuffer);
    mStoreBuffer = nullptr;
  }
  if (mStoreBuffer2 && !Parameters::getInstance().getNoCompressionOverlapFlag())
  {
    _mm_free(mStoreBuffer2);
    mStoreBuffer2 = nullptr;
  }
}// end of freeMemory
//----------------------------------------------------------------------------------------------------------------------

/**
 * Check or set local minimal and maximal value and their indices.
 */
void BaseOutputStream::checkOrSetMinMaxValue(ReducedValue &minValue, ReducedValue &maxValue, float value, hsize_t index)
{
  if (minValue.value > value)
  {
    minValue.value = value;
    minValue.index = index;
  }
  if (maxValue.value < value)
  {
    maxValue.value = value;
    maxValue.index = index;
  }
}// end of checkOrSetMinMaxValue
//----------------------------------------------------------------------------------------------------------------------

/**
 * Check or set global (#pragma omp critical) minimal and maximal value and their indices.
 */
void BaseOutputStream::checkOrSetMinMaxValueGlobal(ReducedValue &minValue, ReducedValue &maxValue, ReducedValue minValueLocal, ReducedValue maxValueLocal)
{
  #pragma omp critical
  {
    if (minValue.value > minValueLocal.value)
    {
      minValue.value = minValueLocal.value;
      minValue.index = minValueLocal.index;
    }
  }
  #pragma omp critical
  {
    if (maxValue.value < maxValueLocal.value)
    {
      maxValue.value = maxValueLocal.value;
      maxValue.index = maxValueLocal.index;
    }
  }
}// end of checkOrSetMinMaxValueGlobal
//----------------------------------------------------------------------------------------------------------------------

/**
 * Load minimal and maximal values from dataset attributes.
 */
void BaseOutputStream::loadMinMaxValues(Hdf5File &file, hid_t group, std::string datasetName, ReducedValue &minValue, ReducedValue &maxValue)
{
  // Reload min and max values
  if (mReduceOp == ReduceOperator::kNone || mReduceOp == ReduceOperator::kC)
  {
    //try {
      minValue.value = file.readFloatAttribute(group, datasetName, "min");
      maxValue.value = file.readFloatAttribute(group, datasetName, "max");
      minValue.index = hsize_t(file.readLongLongAttribute(group, datasetName, "min_index"));
      maxValue.index = hsize_t(file.readLongLongAttribute(group, datasetName, "max_index"));
    //} catch (std::exception &)
    //{}
  }
}// end of loadMinMaxValues
//----------------------------------------------------------------------------------------------------------------------

/**
 * Store minimal and maximal values as dataset attributes.
 */
void BaseOutputStream::storeMinMaxValues(Hdf5File &file, hid_t group, std::string datasetName, ReducedValue minValue, ReducedValue maxValue)
{
  if (mReduceOp == ReduceOperator::kNone || mReduceOp == ReduceOperator::kC)
  {
    file.writeFloatAttribute(group, datasetName, "min", minValue.value);
    file.writeFloatAttribute(group, datasetName, "max", maxValue.value);
    file.writeLongLongAttribute(group, datasetName, "min_index", ssize_t(minValue.index));
    file.writeLongLongAttribute(group, datasetName, "max_index", ssize_t(maxValue.index));
  }
}// end of storeMinMaxValues
//----------------------------------------------------------------------------------------------------------------------

/**
 * Load checkpoint compression coefficients and average intensity.
 */
void BaseOutputStream::loadCheckpointCompressionCoefficients()
{
  if (mReduceOp == ReduceOperator::kC)
  {
    Hdf5File& checkpointFile = Parameters::getInstance().getCheckpointFile();
    checkpointFile.readCompleteDataset(checkpointFile.getRootGroup(),
                                       "Temp_" + mRootObjectName + "_1",
                                       DimensionSizes(mBufferSize, 1, 1),
                                       mStoreBuffer);
    checkpointFile.readCompleteDataset(checkpointFile.getRootGroup(),
                                       "Temp_" + mRootObjectName + "_2",
                                       DimensionSizes(mBufferSize, 1, 1),
                                       mStoreBuffer2);
  }
  if (mReduceOp == ReduceOperator::kIAvgC)
  {
    Hdf5File& checkpointFile = Parameters::getInstance().getCheckpointFile();
    checkpointFile.readCompleteDataset(checkpointFile.getRootGroup(),
                                       "Temp_" + mRootObjectName,
                                       DimensionSizes(mBufferSize, 1, 1),
                                       mStoreBuffer);
  }
}// end of loadCheckpointCompressionCoefficients
//----------------------------------------------------------------------------------------------------------------------

/**
 * Store checkpoint compression coefficients and average intensity.
 */
void BaseOutputStream::storeCheckpointCompressionCoefficients()
{
  // Store temp compression coefficients
  if (mReduceOp == ReduceOperator::kC)
  {
    Hdf5File& checkpointFile = Parameters::getInstance().getCheckpointFile();
    hid_t dataset1 = checkpointFile.createDataset(checkpointFile.getRootGroup(),
                                       "Temp_" + mRootObjectName + "_1",
                                       DimensionSizes(mBufferSize, 1, 1),
                                       DimensionSizes(mBufferSize, 1, 1),
                                       Hdf5File::MatrixDataType::kFloat,
                                       Parameters::getInstance().getCompressionLevel());

    checkpointFile.writeHyperSlab(dataset1, DimensionSizes(0, 0, 0), DimensionSizes(mBufferSize, 1, 1), mStoreBuffer);
    checkpointFile.closeDataset(dataset1);

    hid_t dataset2 = checkpointFile.createDataset(checkpointFile.getRootGroup(),
                                       "Temp_" + mRootObjectName + "_2",
                                       DimensionSizes(mBufferSize, 1, 1),
                                       DimensionSizes(mBufferSize, 1, 1),
                                       Hdf5File::MatrixDataType::kFloat,
                                       Parameters::getInstance().getCompressionLevel());

    checkpointFile.writeHyperSlab(dataset2, DimensionSizes(0, 0, 0), DimensionSizes(mBufferSize, 1, 1), mStoreBuffer2);
    checkpointFile.closeDataset(dataset2);

    // Write data and domain type
    checkpointFile.writeMatrixDataType  (checkpointFile.getRootGroup(), "Temp_" + mRootObjectName + "_1", Hdf5File::MatrixDataType::kFloat);
    checkpointFile.writeMatrixDomainType(checkpointFile.getRootGroup(), "Temp_" + mRootObjectName + "_1", Hdf5File::MatrixDomainType::kReal);
    checkpointFile.writeMatrixDataType  (checkpointFile.getRootGroup(), "Temp_" + mRootObjectName + "_2", Hdf5File::MatrixDataType::kFloat);
    checkpointFile.writeMatrixDomainType(checkpointFile.getRootGroup(), "Temp_" + mRootObjectName + "_2", Hdf5File::MatrixDomainType::kReal);
  }
  // Store temp compression average intensity
  if (mReduceOp == ReduceOperator::kIAvgC)
  {
    Hdf5File& checkpointFile = Parameters::getInstance().getCheckpointFile();
    hid_t dataset1 = checkpointFile.createDataset(checkpointFile.getRootGroup(),
                                                  "Temp_" + mRootObjectName,
                                                  DimensionSizes(mBufferSize, 1, 1),
                                                  DimensionSizes(mBufferSize, 1, 1),
                                                  Hdf5File::MatrixDataType::kFloat,
                                                  Parameters::getInstance().getCompressionLevel());

    checkpointFile.writeHyperSlab(dataset1, DimensionSizes(0, 0, 0), DimensionSizes(mBufferSize, 1, 1), mStoreBuffer);
    checkpointFile.closeDataset(dataset1);
    // Write data and domain type
    checkpointFile.writeMatrixDataType  (checkpointFile.getRootGroup(), "Temp_" + mRootObjectName, Hdf5File::MatrixDataType::kFloat);
    checkpointFile.writeMatrixDomainType(checkpointFile.getRootGroup(), "Temp_" + mRootObjectName, Hdf5File::MatrixDomainType::kReal);
  }
}// end of storeCheckpointCompressionCoefficients
//----------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Private methods --------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//
