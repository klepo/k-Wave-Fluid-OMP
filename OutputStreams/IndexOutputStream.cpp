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

#include <OutputStreams/IndexOutputStream.h>

#include <Parameters/Parameters.h>



//----------------------------------------------------------------------------//
//                 TIndexOutputHDF5Stream implementation                      //
//                              public methods                                //
//----------------------------------------------------------------------------//

/**
 * Constructor - links the HDF5 dataset, source (sampled matrix), Sensor mask
 * and the reduction operator together. The constructor DOES NOT allocate memory
 * because the size of the sensor mask is not known at the time the instance of
 * the class is being created.
 *
 * @param [in] HDF5_File       - Handle to the HDF5 (output) file
 * @param [in] HDF5_ObjectName - The dataset's name (index based sensor data
 *                                is store in a single dataset)
 * @param [in] SourceMatrix    - The source matrix (only real matrices are
 *                               supported)
 * @param [in] SensorMask      - Index based sensor mask
 * @param [in] ReductionOp     - Reduction operator
 * @param [in] BufferToReuse   - An external buffer can be used to line up
 *                               the grid points
 */
TIndexOutputHDF5Stream::TIndexOutputHDF5Stream(Hdf5File &             HDF5_File,
                                               const char *             HDF5_ObjectName,
                                               const RealMatrix &      SourceMatrix,
                                               const IndexMatrix &     SensorMask,
                                               const TReductionOperator ReductionOp,
                                               float *                  BufferToReuse)
        : TBaseOutputHDF5Stream(HDF5_File, HDF5_ObjectName, SourceMatrix, ReductionOp, BufferToReuse),
          SensorMask(SensorMask),
          HDF5_DatasetId(H5I_BADID),
          SampledTimeStep(0)
{

}// end of TIndexOutputHDF5Stream
//------------------------------------------------------------------------------


/**
 * Destructor.
 * If the file is still opened, it applies the post processing and flush the data.
 * Then, the object memory is freed and the object destroyed.
 */
TIndexOutputHDF5Stream::~TIndexOutputHDF5Stream()
{
  Close();
  // free memory only if it was allocated
  if (!BufferReuse) FreeMemory();
}// end of Destructor
//------------------------------------------------------------------------------



/**
 * Create a HDF5 stream, create a dataset, and allocate data for it.
 */
void TIndexOutputHDF5Stream::Create()
{
  size_t NumberOfSampledElementsPerStep = SensorMask.size();

  Parameters& Params = Parameters::getInstance();

  // Derive dataset dimension sizes
  DimensionSizes DatasetSize(NumberOfSampledElementsPerStep,
                              (ReductionOp == roNONE) ?  Params.getNt() - Params.getSamplingStartTimeIndex() : 1,
                              1);

  // Set HDF5 chunk size
  DimensionSizes ChunkSize(NumberOfSampledElementsPerStep, 1, 1);
  // for chunks bigger than 32 MB
  if (NumberOfSampledElementsPerStep > (ChunkSize_4MB * 8))
  {
    ChunkSize.nx = ChunkSize_4MB; // set chunk size to MB
  }

  // Create a dataset under the root group
  HDF5_DatasetId = HDF5_File.createDataset(HDF5_File.getRootGroup(),
                                           HDF5_RootObjectName,
                                           DatasetSize,
                                           ChunkSize,
                                           Hdf5File::MatrixDataType::kFloat,
                                           Params.getCompressionLevel());

  // Write dataset parameters
  HDF5_File.writeMatrixDomainType(HDF5_File.getRootGroup(),
                                  HDF5_RootObjectName,
                                  Hdf5File::MatrixDomainType::kReal);
  HDF5_File.writeMatrixDataType  (HDF5_File.getRootGroup(),
                                  HDF5_RootObjectName,
                                  Hdf5File::MatrixDataType::kFloat);


  // Sampled time step
  SampledTimeStep = 0;

  // Set buffer size
  BufferSize = NumberOfSampledElementsPerStep;

  // Allocate memory if needed
  if (!BufferReuse) AllocateMemory();
}// end of Create
//------------------------------------------------------------------------------

/**
 * Reopen the output stream after restart.
 */
void TIndexOutputHDF5Stream::Reopen()
{
  // Get parameters
  Parameters& Params = Parameters::getInstance();

  // Set buffer size
  BufferSize = SensorMask.size();

  // Allocate memory if needed
  if (!BufferReuse) AllocateMemory();

  // Reopen the dataset
  HDF5_DatasetId = HDF5_File.openDataset(HDF5_File.getRootGroup(),
                                         HDF5_RootObjectName);


  if (ReductionOp == roNONE)
  { // raw time series - just seek to the right place in the dataset
    SampledTimeStep = (Params.getTimeIndex() < Params.getSamplingStartTimeIndex()) ?
                              0 : (Params.getTimeIndex() - Params.getSamplingStartTimeIndex());

  }
  else
  { // aggregated quantities - reload data
    SampledTimeStep = 0;
    // read only if it is necessary (it is anything to read).
    if (Params.getTimeIndex() > Params.getSamplingStartTimeIndex())
    {
      // Since there is only a single timestep in the dataset, I can read the whole dataset
      HDF5_File.readCompleteDataset(HDF5_File.getRootGroup(),
                                    HDF5_RootObjectName,
                                    DimensionSizes(BufferSize, 1, 1),
                                    StoreBuffer);
    }
  }
}// end of Reopen
//------------------------------------------------------------------------------


/**
 * Sample grid points, line them up in the buffer an flush to the disk unless a
 * reduction operator is applied.
 */
void TIndexOutputHDF5Stream::Sample()
{
  const float  * SourceData = SourceMatrix.getData();
  const size_t * SensorData = SensorMask.getData();

  switch (ReductionOp)
  {
    case roNONE :
    {
      #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
      for (size_t i = 0; i < BufferSize; i++)
      {
        StoreBuffer[i] = SourceData[SensorData[i]];
      }
      // only raw time series are flushed down to the disk every time step
      FlushBufferToFile();
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
    }// case roNONE

    case roRMS  :
    {
      #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
      for (size_t i = 0; i < BufferSize; i++)
      {
        StoreBuffer[i] += (SourceData[SensorData[i]] * SourceData[SensorData[i]]);
      }
      break;
    }// case roRMS

    case roMAX  :
    {
      #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
      for (size_t i = 0; i < BufferSize; i++)
      {
        if (StoreBuffer[i] < SourceData[SensorData[i]])
          StoreBuffer[i] = SourceData[SensorData[i]];
      }
      break;
    }// case roMAX

    case roMIN  :
    {
      #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
      for (size_t i = 0; i < BufferSize; i++)
      {
        if (StoreBuffer[i] > SourceData[SensorData[i]])
          StoreBuffer[i] = SourceData[SensorData[i]];
      }
      break;
    } //case roMIN
  }// switch
}// end of Sample
//------------------------------------------------------------------------------

/**
 * Apply post-processing on the buffer and flush it to the file.
 */
void TIndexOutputHDF5Stream::PostProcess()
{
  // run inherited method
  TBaseOutputHDF5Stream::PostProcess();
  // When no reduction operator is applied, the data is flushed after every time step
  if (ReductionOp != roNONE) FlushBufferToFile();
}// end of PostProcessing
//------------------------------------------------------------------------------


/**
 * Checkpoint the stream and close.
 */
void TIndexOutputHDF5Stream::Checkpoint()
{
  // raw data has already been flushed, others has to be fushed here
  if (ReductionOp != roNONE) FlushBufferToFile();

}// end of Checkpoint
//------------------------------------------------------------------------------

/**
 * Close stream (apply post-processing if necessary, flush data and close).
 */
void TIndexOutputHDF5Stream::Close()
{
  // the dataset is still opened
  if (HDF5_DatasetId != H5I_BADID)
  {
    HDF5_File.closeDataset(HDF5_DatasetId);
  }

  HDF5_DatasetId = H5I_BADID;
}// end of Close
//------------------------------------------------------------------------------

//----------------------------------------------------------------------------//
//                 TIndexOutputHDF5Stream implementation                      //
//                            protected methods                               //
//----------------------------------------------------------------------------//


/**
 * Flush the buffer down to the file at the actual position.
 */
void TIndexOutputHDF5Stream::FlushBufferToFile()
{
  HDF5_File.writeHyperSlab(HDF5_DatasetId,
                           DimensionSizes(0,SampledTimeStep,0),
                           DimensionSizes(BufferSize,1,1),
                           StoreBuffer);
  SampledTimeStep++;
}// end of FlushToFile
//------------------------------------------------------------------------------


