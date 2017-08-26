/**
 * @file        WholeDomainOutputStream.cpp
 * @author      Jiri Jaros \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The implementation file of the class saving RealMatrix data into the output
 *              HDF5 file, e.g. p_max_all.
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
#include <OutputStreams/WholeDomainOutputStream.h>
#include <Parameters/Parameters.h>


//----------------------------------------------------------------------------//
//             TWholeDomainOutputHDF5Stream implementation                    //
//                              public methods                                //
//----------------------------------------------------------------------------//

/**
 * Constructor - links the HDF5 dataset and SourceMatrix.
 * @param [in] HDF5_File        - HDF5 file to write the output to
 * @param [in] HDF5_DatasetName - The name of the HDF5 group. This group contains datasets for particular cuboids
 * @param [in] SourceMatrix     - Source matrix to be sampled
 * @param [in] ReductionOp      - Reduction operator
 * @param [in] BufferToReuse    - If there is a memory space to be reused, provide a pointer
 */
TWholeDomainOutputHDF5Stream::TWholeDomainOutputHDF5Stream(Hdf5File &             HDF5_File,
                                                           const char *             HDF5_DatasetName,
                                                           const RealMatrix &      SourceMatrix,
                                                           const TReductionOperator ReductionOp,
                                                           float *                  BufferToReuse)
        : TBaseOutputHDF5Stream(HDF5_File, HDF5_DatasetName, SourceMatrix, ReductionOp, BufferToReuse),
          HDF5_DatasetId(H5I_BADID),
          SampledTimeStep(0)
{

}// end of TWholeDomainOutputHDF5Stream
//------------------------------------------------------------------------------

/**
 * Destructor.
 * if the file is still opened, it applies the post processing and flush the data.
 * Then, the object memory is freed and the object destroyed.
 */
TWholeDomainOutputHDF5Stream::~TWholeDomainOutputHDF5Stream()
{
  Close();
  // free memory only if it was allocated
  if (!BufferReuse) FreeMemory();
}// end of Destructor
//------------------------------------------------------------------------------

/**
 * Create a HDF5 stream for the whole domain and allocate data for it.
 */
void TWholeDomainOutputHDF5Stream::Create()
{
  DimensionSizes ChunkSize(SourceMatrix.getDimensionSizes().nx, SourceMatrix.getDimensionSizes().ny, 1);

  // Create a dataset under the root group
  HDF5_DatasetId = HDF5_File.createDataset(HDF5_File.getRootGroup(),
                                           HDF5_RootObjectName,
                                           SourceMatrix.getDimensionSizes(),
                                           ChunkSize,
                                           Hdf5File::MatrixDataType::kFloat,
                                           Parameters::getInstance().getCompressionLevel());

  // Write dataset parameters
  HDF5_File.writeMatrixDomainType(HDF5_File.getRootGroup(),
                                  HDF5_RootObjectName,
                                  Hdf5File::MatrixDomainType::kReal);
  HDF5_File.writeMatrixDataType  (HDF5_File.getRootGroup(),
                                  HDF5_RootObjectName,
                                  Hdf5File::MatrixDataType::kFloat);

  // Set buffer size
  BufferSize = SourceMatrix.size();

  // Allocate memory if needed
  if (!BufferReuse) AllocateMemory();
}//end of Create
//------------------------------------------------------------------------------


/**
 * Reopen the output stream after restart and reload data.
 */
void TWholeDomainOutputHDF5Stream::Reopen()
{
  Parameters& Params = Parameters::getInstance();

  // Set buffer size
  BufferSize = SourceMatrix.size();

  // Allocate memory if needed
  if (!BufferReuse) AllocateMemory();

  // Open the dataset under the root group
  HDF5_DatasetId = HDF5_File.openDataset(HDF5_File.getRootGroup(),
                                         HDF5_RootObjectName);

  SampledTimeStep = 0;
  if (ReductionOp == roNONE)
  { // seek in the dataset
    SampledTimeStep = (Params.getTimeIndex() < Params.getSamplingStartTimeIndex()) ?
                        0 : (Params.getTimeIndex() - Params.getSamplingStartTimeIndex());
  }
  else
  { // reload data
    if (Params.getTimeIndex() > Params.getSamplingStartTimeIndex())
    {
      HDF5_File.readCompleteDataset(HDF5_File.getRootGroup(),
                                    HDF5_RootObjectName,
                                    SourceMatrix.getDimensionSizes(),
                                    StoreBuffer);
    }
  }
}// end of Reopen
//------------------------------------------------------------------------------


/**
 * Sample all grid points, line them up in the buffer an flush to the disk unless
 * a reduction operator is applied.
 */
void TWholeDomainOutputHDF5Stream::Sample()
{
  const float * SourceData = SourceMatrix.getData();

  switch (ReductionOp)
  {
    case roNONE :
    {

       /* We use here direct HDF5 offload using MEMSPACE - seems to be faster for bigger datasets*/
      const DimensionSizes DatasetPosition(0,0,0,SampledTimeStep); //4D position in the dataset

      DimensionSizes CuboidSize(SourceMatrix.getDimensionSizes());// Size of the cuboid
      CuboidSize.nt = 1;

      // iterate over all cuboid to be sampled
      HDF5_File.writeCuboidToHyperSlab(HDF5_DatasetId,
                                       DatasetPosition,
                                       DimensionSizes(0,0,0,0), // position in the SourceMatrix
                                       CuboidSize,
                                       SourceMatrix.getDimensionSizes(),
                                       SourceMatrix.getData());

      SampledTimeStep++;   // Move forward in time

      // This version is obsolete (needs one more data movement)
      /*
       #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
      for (size_t i = 0; i < BufferSize; i++)
      {
        StoreBuffer[i] = SourceData[i];
      }
      // only raw time series are flushed down to the disk every time step
      FlushBufferToFile();
      */

      break;
    }// case roNONE

    case roRMS  :
    {
      #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
      for (size_t i = 0; i < BufferSize; i++)
      {
        StoreBuffer[i] += (SourceData[i] * SourceData[i]);
      }
      break;
    }// case roRMS

    case roMAX  :
    {
      #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
      for (size_t i = 0; i < BufferSize; i++)
      {
        if (StoreBuffer[i] < SourceData[i])  StoreBuffer[i] = SourceData[i];
      }
      break;
    }//case roMAX

    case roMIN  :
    {
      #pragma omp parallel for if (BufferSize > MinGridpointsToSampleInParallel)
      for (size_t i = 0; i < BufferSize; i++)
      {
        if (StoreBuffer[i] > SourceData[i]) StoreBuffer[i] = SourceData[i];
      }
      break;
    } //case roMIN
  }// switch
}// end of Sample
//------------------------------------------------------------------------------


/**
 * Apply post-processing on the buffer and flush it to the file.
 */
void TWholeDomainOutputHDF5Stream::PostProcess()
{
  // run inherited method
  TBaseOutputHDF5Stream::PostProcess();
  // When no reduction operator is applied, the data is flushed after every time step
  if (ReductionOp != roNONE) FlushBufferToFile();
}// end of PostProcessing
//------------------------------------------------------------------------------

/**
 * Checkpoint the stream
 */
void TWholeDomainOutputHDF5Stream::Checkpoint()
{
  // raw data has already been flushed, others has to be flushed here.
  if (ReductionOp != roNONE) FlushBufferToFile();
}// end of Checkpoint
//------------------------------------------------------------------------------

/**
 * Close stream (apply post-processing if necessary, flush data and close).
 */
void TWholeDomainOutputHDF5Stream::Close()
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
//                TWholeDomainOutputHDF5Stream implementation                 //
//                            protected methods                               //
//----------------------------------------------------------------------------//


/**
 * Flush the buffer down to the file at the actual position.
 */
void TWholeDomainOutputHDF5Stream::FlushBufferToFile()
{
  DimensionSizes Size = SourceMatrix.getDimensionSizes();
  DimensionSizes Position(0,0,0);

  // Not used for roNONE now!
  if (ReductionOp == roNONE)
  {
    Position.nt = SampledTimeStep;
    Size.nt = SampledTimeStep;
  }

  HDF5_File.writeHyperSlab(HDF5_DatasetId,
                           Position,
                           Size,
                           StoreBuffer);
  SampledTimeStep++;
}// end of FlushToFile
//------------------------------------------------------------------------------

