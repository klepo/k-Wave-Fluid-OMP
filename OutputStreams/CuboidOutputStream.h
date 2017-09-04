/**
 * @file      CuboidOutputStream.h
 *
 * @author    Jiri Jaros \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            jarosjir@fit.vutbr.cz
 *
 * @brief     The header file of classes responsible for storing output quantities based on the
 *            cuboid sensor mask into the output HDF5 file.
 *
 * @version   kspaceFirstOrder3D 2.16
 *
 * @date      26 August    2017, 16:55 (created) \n
 *            04 September 2017, 11:10 (revised)
 *
 * @copyright Copyright (C) 2017 Jiri Jaros and Bradley Treeby.
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

#ifndef CUBOID_OUTPUT_STREAM_H
#define CUBOID_OUTPUT_STREAM_H

#include <vector>

#include <OutputStreams/BaseOutputStream.h>


/**
 * @class CuboidOutputStream
 * @brief Output stream for quantities sampled by a cuboid corner sensor mask.
 *
 * Output stream for quantities sampled by a cuboid corner sensor mask. This class writes data into separated datasets
 * (one per cuboid) under a given dataset in the HDF5 file (time-series as well as aggregations).
 *
 */
class CuboidOutputStream : public BaseOutputStream
{
  public:
    /// Default constructor not allowed
    CuboidOutputStream() = delete;

    /**
     * @brief  Constructor links the HDF5 dataset, SourceMatrix, and SensorMask together.
     *
     * @param [in] file          - HDF5 file to write the output to.
     * @param [in] groupName     - The name of the HDF5 group. This group contains datasets for particular cuboids.
     * @param [in] sourceMatrix  - Source matrix to be sampled.
     * @param [in] sensorMask    - Sensor mask with the cuboid coordinates.
     * @param [in] reduceOp      - Reduction operator.
     * @param [in] bufferToReuse - If there is a memory space to be reused, provide a pointer.
     */
    CuboidOutputStream(Hdf5File&            file,
                       MatrixName&          groupName,
                       const RealMatrix&    sourceMatrix,
                       const IndexMatrix&   sensorMask,
                       const ReduceOperator ReduceOp,
                       float*               bufferToReuse = nullptr);

    /// Copy constructor is not allowed.
    CuboidOutputStream(const CuboidOutputStream&) = delete;

    /**
     * @brief Destructor.
     *
     * If the file is still opened, it applies the post processing and flush the data.
     * Then, the object memory is freed and the object destroyed.
     */
    virtual ~CuboidOutputStream();

    /// operator= is not allowed.
    CuboidOutputStream& operator=(const CuboidOutputStream&) = delete;

    /// Create a HDF5 stream and allocate data for it.
    virtual void create();

    /// Reopen the output stream after restart and reload data.
    virtual void reopen();

    /// Sample data into buffer and apply reduction, or flush to disk - based on a sensor mask.
    virtual void sample();

    /// Apply post-processing on the buffer and flush it to the file.
    virtual void postProcess();

    /// Checkpoint the stream and close.
    virtual void checkpoint();

    /// Close stream (apply post-processing if necessary, flush data and close).
    virtual void close();

  protected:
    /**
     * @struct CuboidInfo
     * @brief  This structure information about a HDF5 dataset (one cuboid).
     */
    struct CuboidInfo
    {
      /// Id of the dataset storing the given cuboid.
      hid_t  cuboidId;
      /// Having a single buffer for all cuboids, where this one starts.
      size_t startingPossitionInBuffer;
    };

    /**
     * @brief Create a new dataset for a given cuboid specified by index (order).
     * @param [in] cuboidIdx - Index of the cuboid in the sensor mask.
     * @return Handle to the HDF5 dataset.
     */
    virtual hid_t createCuboidDataset(const size_t cuboidIdx);

    /**
     * @brief  Sample aggregated values.
     * @tparam reduceOp - Reduction operator
     */
    template<BaseOutputStream::ReduceOperator reduceOp>
    void sampleAggregated();

    /// Flush the buffer to the file.
    virtual void flushBufferToFile();

    /// Sensor mask to sample data.
    const IndexMatrix&      mSensorMask;

    /// Handle to a HDF5 dataset.
    hid_t                   mGroup;

    /// vector keeping handles and positions of all cuboids.
    std::vector<CuboidInfo> mCuboidsInfo;

    /// Timestep to store (N/A for aggregated).
    size_t                  mSampledTimeStep;

};// end of CuboidOutputStream
//----------------------------------------------------------------------------------------------------------------------

#endif	/* CUBOID_OUTPUT_STREAM_H */

