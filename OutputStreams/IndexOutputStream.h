/**
 * @file        IndexOutputStream.h
 * @author      Jiri Jaros \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The header file of the class saving data based on the index senor mask into the output HDF5 file.
 *
 * @version     kspaceFirstOrder3D 2.16
 *
 * @date        26 August    2017, 16:55 (created) \n
 *              26 August    2017, 16:55 (revised)
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

#ifndef INDEX_OUTPUT_STREAM_H
#define INDEX_OUTPUT_STREAM_H

#include <OutputStreams/BaseOutputStream.h>


/**
 * @class TIndexOutputHDF5Stream.
 * @brief   Output stream for quantities sampled by an index sensor mask.
 * @details Output stream for quantities sampled by an index sensor mask.
 *        This class writes data to a single dataset in a root group of the HDF5
 *        file (time-series as well as aggregations).
 *
 */
class TIndexOutputHDF5Stream : public TBaseOutputHDF5Stream
{
  public:

    /// Constructor - links the HDF5 dataset, SourceMatrix, and SensorMask together.
    TIndexOutputHDF5Stream(Hdf5File &             HDF5_File,
                           const char *             HDF5_ObjectName,
                           const RealMatrix &      SourceMatrix,
                           const IndexMatrix &     SensorMask,
                           const TReductionOperator ReductionOp,
                           float *                  BufferToReuse = NULL);


    /// Destructor.
    virtual ~TIndexOutputHDF5Stream();

    /// Create a HDF5 stream and allocate data for it.
    virtual void Create();

    /// Reopen the output stream after restart and reload data.
    virtual void Reopen();

    /// Sample data into buffer, apply reduction or flush to disk - based on a sensor mask.
    virtual void Sample();

    /// Apply post-processing on the buffer and flush it to the file.
    virtual void PostProcess();

    /// Checkpoint the stream.
    virtual void Checkpoint();

    /// Close stream (apply post-processing if necessary, flush data and close).
    virtual void Close();

  protected:

    /// Flush the buffer to the file.
    virtual void FlushBufferToFile();

    /// Sensor mask to sample data.
    const IndexMatrix & SensorMask;
    /// Handle to a HDF5 dataset.
    hid_t  HDF5_DatasetId;

    /// Time step to store (N/A for aggregated).
    size_t SampledTimeStep;
}; // end of TIndexOutputHDF5Stream
//------------------------------------------------------------------------------


#endif	/* INDEX_OUTPUT_STREAM_H */

