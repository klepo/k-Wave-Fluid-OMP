/**
 * @file      WholeDomainOutputStream.cpp
 *
 * @author    Jiri Jaros \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            jarosjir@fit.vutbr.cz
 *
 * @brief     The implementation file of the class saving RealMatrix data into the output
 *            HDF5 file, e.g. p_max_all.
 *
 * @version   kspaceFirstOrder 2.17
 *
 * @date      26 August    2017, 17:03 (created) \n
 *            08 February  2023, 12:00 (revised)
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

#include <OutputStreams/WholeDomainOutputStream.h>
#include <Parameters/Parameters.h>

//--------------------------------------------------------------------------------------------------------------------//
//--------------------------------------------------- Constants ------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Public methods ---------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Constructor - links the HDF5 dataset and SourceMatrix.
 * @param [in] file          - HDF5 file to write the output to
 * @param [in] datasetName   - The name of the HDF5 group. This group contains datasets for particular cuboids
 * @param [in] sourceMatrix  - Source matrix to be sampled
 * @param [in] reductionOp   - Reduction operator
 * @param [in] bufferToReuse - If there is a memory space to be reused, provide a pointer
 */
WholeDomainOutputStream::WholeDomainOutputStream(Hdf5File& file,
  MatrixName& datasetName,
  const RealMatrix& sourceMatrix,
  const ReduceOperator reductionOp,
  float* bufferToReuse,
  OutputStreamContainer* outputStreamContainer,
  bool doNotSaveFlag)
  : BaseOutputStream(file, datasetName, sourceMatrix, reductionOp, bufferToReuse, outputStreamContainer, doNotSaveFlag),
    mDataset(H5I_BADID), mSampledTimeStep(0)
{

} // end of WholeDomainOutputStream
//----------------------------------------------------------------------------------------------------------------------

/**
 * Destructor.
 */
WholeDomainOutputStream::~WholeDomainOutputStream()
{
  close();
  // free memory only if it was allocated
  if (!mBufferReuse)
    freeMemory();
} // end of Destructor
//----------------------------------------------------------------------------------------------------------------------

/**
 * Create a HDF5 stream for the whole domain and allocate data for it.
 */
void WholeDomainOutputStream::create()
{
  DimensionSizes chunkSize(mSourceMatrix.getDimensionSizes().nx, mSourceMatrix.getDimensionSizes().ny, 1);

  // Create a dataset under the root group
  mDataset = mFile.createDataset(mFile.getRootGroup(),
    mRootObjectName,
    mSourceMatrix.getDimensionSizes(),
    chunkSize,
    Hdf5File::MatrixDataType::kFloat,
    Parameters::getInstance().getCompressionLevel());

  // Write dataset parameters
  mFile.writeMatrixDomainType(mFile.getRootGroup(), mRootObjectName, Hdf5File::MatrixDomainType::kReal);
  mFile.writeMatrixDataType(mFile.getRootGroup(), mRootObjectName, Hdf5File::MatrixDataType::kFloat);

  // Set buffer size
  mBufferSize = mSourceMatrix.size();

  // Allocate memory if needed
  if (!mBufferReuse)
    allocateMemory();
} // end of create
//----------------------------------------------------------------------------------------------------------------------

/**
 * Reopen the output stream after restart and reload data.
 */
void WholeDomainOutputStream::reopen()
{
  const Parameters& params = Parameters::getInstance();

  // Set buffer size
  mBufferSize = mSourceMatrix.size();

  // Allocate memory if needed
  if (!mBufferReuse)
    allocateMemory();

  // Open the dataset under the root group
  mDataset = mFile.openDataset(mFile.getRootGroup(), mRootObjectName);

  mSampledTimeStep = 0;
  if (mReduceOp == ReduceOperator::kNone)
  { // seek in the dataset
    mSampledTimeStep = (params.getTimeIndex() < params.getSamplingStartTimeIndex())
                         ? 0
                         : (params.getTimeIndex() - params.getSamplingStartTimeIndex());
  }
  else
  { // reload data
    if (params.getTimeIndex() > params.getSamplingStartTimeIndex())
    {
      mFile.readCompleteDataset(mFile.getRootGroup(), mRootObjectName, mSourceMatrix.getDimensionSizes(), mStoreBuffer);
    }
  }
} // end of reopen
//----------------------------------------------------------------------------------------------------------------------

/**
 * Sample all grid points, line them up in the buffer an flush to the disk unless a reduction operator is applied.
 */
void WholeDomainOutputStream::sample()
{
  const float* sourceData = mSourceMatrix.getData();

  switch (mReduceOp)
  {
  case ReduceOperator::kNone:
  {
    // We sample it as a single cuboid of full dimensions.
    /* We use here direct HDF5 offload using MEMSPACE - seems to be faster for bigger datasets*/
    const DimensionSizes datasetPosition(0, 0, 0, mSampledTimeStep); // 4D position in the dataset

    DimensionSizes cuboidSize(mSourceMatrix.getDimensionSizes()); // Size of the cuboid
    cuboidSize.nt = 1;

    // iterate over all cuboid to be sampled
    mFile.writeCuboidToHyperSlab(mDataset,
      datasetPosition,
      DimensionSizes(0, 0, 0, 0), // position in the SourceMatrix
      cuboidSize,
      mSourceMatrix.getDimensionSizes(),
      mSourceMatrix.getData());

    mSampledTimeStep++; // Move forward in time

    break;
  } // case kNone

  case ReduceOperator::kRms:
  {
#pragma omp parallel for simd
    for (size_t i = 0; i < mBufferSize; i++)
    {
      mStoreBuffer[i] += (sourceData[i] * sourceData[i]);
    }
    break;
  } // case kRms

  case ReduceOperator::kMax:
  {
#pragma omp parallel for simd
    for (size_t i = 0; i < mBufferSize; i++)
    {
      mStoreBuffer[i] = std::max(mStoreBuffer[i], sourceData[i]);
    }
    break;
  } // case roMAX

  case ReduceOperator::kMin:
  {
#pragma omp parallel for simd
    for (size_t i = 0; i < mBufferSize; i++)
    {
      mStoreBuffer[i] = std::min(mStoreBuffer[i], sourceData[i]);
    }
    break;
  } // case kMin

  default:
  {
    break;
  }
  } // switch
} // end of sample
//----------------------------------------------------------------------------------------------------------------------

/**
 * Post sampling step, can work with other filled stream buffers
 */
void WholeDomainOutputStream::postSample()
{

} // end of postSample
//----------------------------------------------------------------------------------------------------------------------

/**
 * Apply post-processing on the buffer and flush it to the file.
 */
void WholeDomainOutputStream::postProcess()
{
  // run inherited method
  BaseOutputStream::postProcess();
  // When no reduction operator is applied, the data is flushed after every time step
  if (mReduceOp != ReduceOperator::kNone)
    flushBufferToFile();
} // end of postProcessing
//----------------------------------------------------------------------------------------------------------------------

/**
 * Checkpoint the stream
 */
void WholeDomainOutputStream::checkpoint()
{
  // raw data has already been flushed, others has to be flushed here.
  if (mReduceOp != ReduceOperator::kNone)
    flushBufferToFile();
} // end of checkpoint
//----------------------------------------------------------------------------------------------------------------------

/**
 * Close stream (apply post-processing if necessary, flush data and close).
 */
void WholeDomainOutputStream::close()
{
  // the dataset is still opened
  if (mDataset != H5I_BADID)
  {
    mFile.closeDataset(mDataset);
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
void WholeDomainOutputStream::flushBufferToFile()
{
  DimensionSizes size = mSourceMatrix.getDimensionSizes();
  DimensionSizes position(0, 0, 0);

  // Not used for roNONE now!
  if (mReduceOp == ReduceOperator::kNone)
  {
    position.nt = mSampledTimeStep;
    size.nt     = mSampledTimeStep;
  }

  mFile.writeHyperSlab(mDataset, position, size, mStoreBuffer);
  mSampledTimeStep++;
} // end of flushBufferToFile
//----------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Private methods --------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//
