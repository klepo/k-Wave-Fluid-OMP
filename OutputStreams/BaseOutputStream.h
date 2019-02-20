/**
 * @file      BaseOutputStream.h
 *
 * @author    Jiri Jaros \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            jarosjir@fit.vutbr.cz
 *
 * @brief     The header file of the class saving RealMatrix data into the output HDF5 file.
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

#ifndef BASE_OUTPUT_STREAM_H
#define BASE_OUTPUT_STREAM_H


#include <MatrixClasses/RealMatrix.h>
#include <MatrixClasses/IndexMatrix.h>
#include <Hdf5/Hdf5File.h>

/**
 * @class   BaseOutputStream
 * @brief   Abstract base class for output data streams (sampled data).
 *
 * Data are sampled based on the the sensor mask and the reduction operator. The sampled data is stored in the output
 * HDF5 file.
 */
class BaseOutputStream
{
  public:

    /**
     * @enum  ReduceOperator
     * @brief How to aggregate data.
     */
    enum class ReduceOperator
    {
      /// Store actual data (time series).
      kNone,
       /// Calculate root mean square.
      kRms,
      /// Store maximum.
      kMax,
      /// Store minimum.
      kMin
    };

    /// Default constructor not allowed.
    BaseOutputStream() = delete;

    /**
     * @brief Constructor
     *
     * There is no sensor mask by default to support both sensor mask and whole domain sampling.
     * The constructor links the HDF5 dataset, source (sampled matrix) and the reduction operator together. The
     * constructor DOES NOT allocate memory because the size of the sensor mask is not known at the time the instance
     * of the class is being created.
     *
     * @param [in] file           - Handle to the output HDF5 file.
     * @param [in] rootObjectName - The root object that stores the sample data (dataset or group).
     * @param [in] sourceMatrix   - The source matrix (only real matrices  are supported).
     * @param [in] reduceOp       - Reduction operator.
     * @param [in] bufferToReuse  - An external buffer can be used to line up the grid points.
     */
    BaseOutputStream(Hdf5File&            file,
                     MatrixName&          rootObjectName,
                     const RealMatrix&    sourceMatrix,
                     const ReduceOperator reduceOp,
                     float*               bufferToReuse = nullptr);

    /// Copy constructor not allowed.
    BaseOutputStream(const BaseOutputStream& src);
    /**
     * @brief Destructor.
     *
     * If the file is still opened, it applies the post processing and flush the data.
     * Then, the object memory is freed and the object destroyed.
     */
    virtual ~BaseOutputStream() {};

    /// Operator = not allowed (we don't want any data movements).
    BaseOutputStream& operator = (const BaseOutputStream& src);

    /// Create a HDF5 stream and allocate data for it.
    virtual void create() = 0;

    /// Reopen the output stream after restart.
    virtual void reopen() = 0;

    /// Sample data into buffer, apply reduction or flush to disk - based on a sensor mask.
    virtual void sample() = 0;

    /// Apply post-processing on the buffer and flush it to the file.
    virtual void postProcess();

    /// Checkpoint the stream.
    virtual void checkpoint() = 0;

    /// Close stream (apply post-processing if necessary, flush data and close).
    virtual void close() = 0;

  protected:

    /**
     * @brief    Allocate memory using proper memory alignment.
     * @throw    std::bad_alloc - If there's not enough memory.
     * @warning  This can routine is not used in the base class (should be used in derived ones).
     */
    virtual void allocateMemory();
    /**
     * @brief   Free memory.
     * @warning This can routine is not used in the base class (should be used in derived ones).
     */
    virtual void freeMemory();

    /// Handle to HDF5 output file.
    Hdf5File&            mFile;
    /// HDF5 group/dataset in the output file where to store data in.
    std::string          mRootObjectName;
    /// Source matrix to be sampled.
    const RealMatrix&    mSourceMatrix;
    /// Reduction operator.
    const ReduceOperator mReduceOp;

    /// if true, the container reuses another matrix as scratch place, e.g. kTemp1RealND, kTemp1RealND, kTemp1RealND.
    bool   mBufferReuse;
    /// Buffer size.
    size_t mBufferSize;
    /// Temporary buffer for store - only if Buffer Reuse = false!
    float* mStoreBuffer;

    /// chunk size of 4MB in number of float elements.
    static constexpr size_t kChunkSize4MB = 1048576;
};// end of BaseOutputStream
//----------------------------------------------------------------------------------------------------------------------

#endif	/* BASE_OUTPUT_STREAM_H */
