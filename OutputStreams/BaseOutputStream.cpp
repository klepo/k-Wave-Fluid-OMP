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


//--------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------- Constants -----------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Public methods ---------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Constructor - there is no sensor mask by default!
 */
BaseOutputStream::BaseOutputStream(Hdf5File&            file,
                                   MatrixName&          rootObjectName,
                                   const RealMatrix&    sourceMatrix,
                                   const ReduceOperator reduceOp,
                                   float*               bufferToReuse)
  : mFile(file),
    mRootObjectName(rootObjectName),
    mSourceMatrix(sourceMatrix),
    mReduceOp(reduceOp),
    mBufferReuse(bufferToReuse != nullptr),
    mBufferSize(0),
    mStoreBuffer(bufferToReuse)
{
  // Set compression variables
  if (mReduceOp == ReduceOperator::kC)
  {
    // Set compression helper
    mCompressHelper = &CompressHelper::getInstance();
  }
};
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



//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------ Protected methods -------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Allocate memory using a proper memory alignment.
 */
void BaseOutputStream::allocateMemory()
{
  mStoreBuffer = (float*) _mm_malloc(mBufferSize * sizeof(float), kDataAlignment);

  if (!mStoreBuffer)
  {
    throw std::bad_alloc();
  }

  if (mReduceOp == ReduceOperator::kC)
  {
    mStoreBuffer2 = (float*) _mm_malloc(mBufferSize * sizeof(float), kDataAlignment);

    if (!mStoreBuffer2)
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
  if (mStoreBuffer2)
  {
    _mm_free(mStoreBuffer2);
    mStoreBuffer2 = nullptr;
  }
  if (minValue)
  {
    _mm_free(minValue);
    minValue = nullptr;
  }
  if (maxValue)
  {
    _mm_free(maxValue);
    maxValue = nullptr;
  }
  if (minValueIndex)
  {
    _mm_free(minValueIndex);
    minValueIndex = nullptr;
  }
  if (maxValueIndex)
  {
    _mm_free(maxValueIndex);
    maxValueIndex = nullptr;
  }
}// end of freeMemory
//----------------------------------------------------------------------------------------------------------------------

void BaseOutputStream::checkOrSetMinMaxValue(float &minV, float &maxV, float value, hsize_t &minVIndex, hsize_t &maxVIndex, hsize_t index)
{
  if (minV > value)
  { // TODO: think about this
    #pragma omp critical
    {
      if (minV > value)
      {
        minV = value;
        minVIndex = index;
      }
    }
  }

  if (maxV < value)
  { // TODO: think about this
    #pragma omp critical
    {
      if (maxV < value)
      {
        maxV = value;
        maxVIndex = index;
      }
    }
  }
}

void BaseOutputStream::allocateMinMaxMemory(hsize_t items)
{
  this->items = items;
  maxValue = (float*) _mm_malloc(items * sizeof(float), kDataAlignment);
  minValue = (float*) _mm_malloc(items * sizeof(float), kDataAlignment);
  maxValueIndex = (hsize_t*) _mm_malloc(items * sizeof(hsize_t), kDataAlignment);
  minValueIndex = (hsize_t*) _mm_malloc(items * sizeof(hsize_t), kDataAlignment);

  for (size_t i = 0; i < items; i++)
  {
    maxValue[i] = std::numeric_limits<float>::min();
    minValue[i] = std::numeric_limits<float>::max();
    maxValueIndex[i] = 0;
    minValueIndex[i] = 0;
  }
}

void BaseOutputStream::loadMinMaxValues(Hdf5File &file, hid_t group, std::string datasetName, size_t index, bool checkpoint)
{
  std::string suffix = checkpoint ? "_" + std::to_string(index) : "";

  // Reload min and max values
  if (mReduceOp == ReduceOperator::kNone || mReduceOp == ReduceOperator::kC)
  {
    //try {
      minValue[index] = file.readFloatAttribute(group, datasetName, "min" + suffix);
      maxValue[index] = file.readFloatAttribute(group, datasetName, "max" + suffix);
      minValueIndex[index] = hsize_t(file.readLongLongAttribute(group, datasetName, "min_index" + suffix));
      maxValueIndex[index] = hsize_t(file.readLongLongAttribute(group, datasetName, "max_index" + suffix));
    //} catch (std::exception &) {
    //}
  }
}

void BaseOutputStream::storeMinMaxValues(Hdf5File &file, hid_t group, std::string datasetName, size_t index, bool checkpoint)
{
  std::string suffix = checkpoint ? "_" + std::to_string(index) : "";
  if (mReduceOp == ReduceOperator::kNone || mReduceOp == ReduceOperator::kC)
  {
    file.writeFloatAttribute(group, datasetName, "min" + suffix, minValue[index]);
    file.writeFloatAttribute(group, datasetName, "max" + suffix, maxValue[index]);
    file.writeLongLongAttribute(group, datasetName, "min_index" + suffix, ssize_t(minValueIndex[index]));
    file.writeLongLongAttribute(group, datasetName, "max_index" + suffix, ssize_t(maxValueIndex[index]));
  }
}

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
}

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
}

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Private methods --------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//
