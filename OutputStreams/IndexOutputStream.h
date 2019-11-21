/**
 * @file      IndexOutputStream.h
 *
 * @author    Jiri Jaros \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            jarosjir@fit.vutbr.cz
 *
 * @brief     The header file of the class saving data based on the index senor mask into the output HDF5 file.
 *
 * @version   kspaceFirstOrder 2.17
 *
 * @date      26 August    2017, 16:55 (created) \n
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

#ifndef INDEX_OUTPUT_STREAM_H
#define INDEX_OUTPUT_STREAM_H

#include <OutputStreams/BaseOutputStream.h>

/**
 * @class   IndexOutputStream.
 * @brief   Output stream for quantities sampled by an index sensor mask.
 *
 * Output stream for quantities sampled by an index sensor mask. This class writes data to a single dataset in a
 * root group of the HDF5 file (time-series as well as aggregations).
 */
class IndexOutputStream : public BaseOutputStream
{
  public:
    /// Default constructor not allowed.
    IndexOutputStream() = delete;

    /**
     * @brief Constructor links the HDF5 dataset, SourceMatrix, and SensorMask together.
     *
     * Constructor - links the HDF5 dataset, source (sampled matrix), Sensor mask and the reduction operator together.
     * The constructor DOES NOT allocate memory because the size of the sensor mask is not known at the time the
     * instance of the class is being created.
     *
     * @param [in] file          - Handle to the HDF5 (output) file.
     * @param [in] datasetName   - The dataset's name (index based sensor data is store in a single dataset).
     * @param [in] sourceMatrix  - The source matrix (only real matrices are supported).
     * @param [in] sensorMask    - Index based sensor mask.
     * @param [in] reduceOp      - Reduction operator.
     * @param [in] bufferToReuse - An external buffer can be used to line up the grid points.
     */
    IndexOutputStream(Hdf5File&               file,
                      MatrixName&             datasetName,
                      const RealMatrix&       sourceMatrix,
                      const IndexMatrix&      sensorMask,
                      const ReduceOperator    reduceOp,
                      float*                  bufferToReuse = nullptr,
                      OutputStreamContainer*  outputStreamContainer = nullptr,
                      bool                    doNotSaveFlag = false);

    /// Copy constructor not allowed.
    IndexOutputStream(const IndexOutputStream&) = delete;

    /**
     * @brief Destructor.
     *
     * If the file is still opened, it applies the post processing and flush the data.
     * Then, the object memory is freed and the object destroyed.
     */
    virtual ~IndexOutputStream();

    /// Create a HDF5 stream and allocate data for it.
    virtual void create();

    /// Reopen the output stream after restart and reload data.
    virtual void reopen();

    /// Sample data into buffer, apply reduction or flush to disk - based on a sensor mask.
    virtual void sample();

    /// Post sampling step, can work with other filled stream buffers
    virtual void postSample();

    /// Apply post-processing on the buffer and flush it to the file.
    virtual void postProcess();

    /// Apply post-processing 2 on the buffer and flush it to the file.
    virtual void postProcess2();

    /// Checkpoint the stream.
    virtual void checkpoint();

    /// Close stream (apply post-processing if necessary, flush data and close).
    virtual void close();

  protected:
    /// Flush the buffer to the file.
    virtual void flushBufferToFile(float *bufferToFlush = nullptr);

    /// Sensor mask to sample data.
    const IndexMatrix& sensorMask;
    /// Handle to a HDF5 dataset.
    hid_t  mDataset;

    /// Time step to store (N/A for aggregated).
    size_t mSampledTimeStep;

    /// Maximal value
    ReducedValue    mMaxValue;
    /// Minimal value
    ReducedValue    mMinValue;

}; // end of IndexOutputStream
//----------------------------------------------------------------------------------------------------------------------

#endif	/* INDEX_OUTPUT_STREAM_H */
