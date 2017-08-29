/**
 * @file        IndexOutputStream.cpp
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The implementation file of the class saving data based on index senor mask into
 *              the output HDF5 file.
 *
 * @version     kspaceFirstOrder3D 2.16
 *
 * @date        26 August    2017, 17:03 (created) \n
 *              29 August    2017, 10:18 (revised)
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

#include <OutputStreams/IndexOutputStream.h>
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
IndexOutputStream::IndexOutputStream(Hdf5File&            file,
                                     MatrixName&          datasetName,
                                     const RealMatrix&    sourceMatrix,
                                     const IndexMatrix&   sensorMask,
                                     const ReduceOperator reduceOp,
                                     float *              bufferToReuse)
  : BaseOutputStream(file, datasetName, sourceMatrix, reduceOp, bufferToReuse),
    sensorMask(sensorMask),
    mDataset(H5I_BADID),
    mSampledTimeStep(0)
{

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
  size_t nSampledElementsPerStep = sensorMask.size();

  const Parameters& params = Parameters::getInstance();

  // Derive dataset dimension sizes
  DimensionSizes datasetSize(nSampledElementsPerStep,
                             (mReduceOp == ReduceOperator::kNone) ?
                               params.getNt() - params.getSamplingStartTimeIndex() : 1,
                             1);

  // Set HDF5 chunk size
  DimensionSizes chunkSize(nSampledElementsPerStep, 1, 1);
  // for data bigger than 32 MB
  if (nSampledElementsPerStep > (kChunkSize4MB * 8))
  {
    chunkSize.nx = kChunkSize4MB; // set chunk size to MB
  }

  // Create a dataset under the root group
  mDataset = mFile.createDataset(mFile.getRootGroup(),
                                 mRootObjectName,
                                 datasetSize,
                                 chunkSize,
                                 Hdf5File::MatrixDataType::kFloat,
                                 params.getCompressionLevel());

  // Write dataset parameters
  mFile.writeMatrixDomainType(mFile.getRootGroup(), mRootObjectName, Hdf5File::MatrixDomainType::kReal);
  mFile.writeMatrixDataType  (mFile.getRootGroup(), mRootObjectName, Hdf5File::MatrixDataType::kFloat);

  // Sampled time step
  mSampledTimeStep = 0;

  // Set buffer size
  mBufferSize = nSampledElementsPerStep;

  // Allocate memory if needed
  if (!mBufferReuse) allocateMemory();
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

  // Allocate memory if needed
  if (!mBufferReuse) allocateMemory();

  // Reopen the dataset
  mDataset = mFile.openDataset(mFile.getRootGroup(), mRootObjectName);

  if (mReduceOp == ReduceOperator::kNone)
  { // raw time series - just seek to the right place in the dataset
    mSampledTimeStep = (params.getTimeIndex() < params.getSamplingStartTimeIndex()) ?
                          0 : (params.getTimeIndex() - params.getSamplingStartTimeIndex());
  }
  else
  { // aggregated quantities - reload data
    mSampledTimeStep = 0;
    // read only if it is necessary (it is anything to read).
    if (params.getTimeIndex() > params.getSamplingStartTimeIndex())
    {
      // Since there is only a single timestep in the dataset, I can read the whole dataset
      mFile.readCompleteDataset(mFile.getRootGroup(),
                                mRootObjectName,
                                DimensionSizes(mBufferSize, 1, 1),
                                mStoreBuffer);
    }
  }
}// end of reopen
//----------------------------------------------------------------------------------------------------------------------

/**
 * Sample grid points, line them up in the buffer an flush to the disk unless a reduction operator is applied.
 */
void IndexOutputStream::sample()
{
  const float*  sourceData = mSourceMatrix.getData();
  const size_t* sensorData = sensorMask.getData();

  switch (mReduceOp)
  {
    case ReduceOperator::kNone :
    {
      #pragma omp parallel for
      for (size_t i = 0; i < mBufferSize; i++)
      {
        mStoreBuffer[i] = sourceData[sensorData[i]];
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

    case ReduceOperator::kRms  :
    {
      #pragma omp parallel for
      for (size_t i = 0; i < mBufferSize; i++)
      {
        mStoreBuffer[i] += (sourceData[sensorData[i]] * sourceData[sensorData[i]]);
      }
      break;
    }// case kRms

    case ReduceOperator::kMax  :
    {
      #pragma omp parallel for
      for (size_t i = 0; i < mBufferSize; i++)
      {
        mStoreBuffer[i] = std::max(mStoreBuffer[i], sourceData[sensorData[i]]);
      }
      break;
    }// case kMax

    case ReduceOperator::kMin  :
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
 * Apply post-processing on the buffer and flush it to the file.
 */
void IndexOutputStream::postProcess()
{
  // run inherited method
  BaseOutputStream::postProcess();
  // When no reduction operator is applied, the data is flushed after every time step
  if (mReduceOp != ReduceOperator::kNone)
  {
    flushBufferToFile();
  }
}// end of postProcessing
//----------------------------------------------------------------------------------------------------------------------

/**
 * Checkpoint the stream and close.
 */
void IndexOutputStream::checkpoint()
{
  // raw data has already been flushed, others has to be flushed here
  if (mReduceOp != ReduceOperator::kNone)
  {
    flushBufferToFile();
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
    mFile.closeDataset(mDataset);
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
void IndexOutputStream::flushBufferToFile()
{
  mFile.writeHyperSlab(mDataset,
                       DimensionSizes(0, mSampledTimeStep, 0),
                       DimensionSizes(mBufferSize, 1, 1),
                       mStoreBuffer);
  mSampledTimeStep++;
}// end of flushToFile
//----------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Private methods --------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

