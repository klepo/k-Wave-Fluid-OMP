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
}// end of freeMemory
//----------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Private methods --------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//
