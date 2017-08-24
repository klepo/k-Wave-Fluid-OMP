/**
 * @file        OutputHDF5Stream.h
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The header file of classes responsible for storing output
 *              quantities into the output HDF5 file.
 *
 * @version     kspaceFirstOrder3D 2.16
 *
 * @date        11 July      2012, 10:30 (created) \n
 *              24 August    2017, 12:21 (revised)
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


#ifndef OUTPUTHDF5STREAM_H
#define	OUTPUTHDF5STREAM_H

#include <string>
#include <vector>
#include <stdexcept>

#include <MatrixClasses/RealMatrix.h>
#include <MatrixClasses/IndexMatrix.h>

#include <Hdf5/Hdf5File.h>

using namespace std;

/**
 * @class TBaseOutputHDF5Stream
 * @brief   Abstract base class for output data streams (sampled data).
 * @details Abstract base class for output data streams (sampled data).
 *
 */
class TBaseOutputHDF5Stream
{
  public:

    /**
     * @enum TReductionOperator
     * @brief How to aggregate data.
     * @details How to aggregate data \n
     *           roNONE - store actual data (time series)
     *           roRMS  - calculate root mean square \n
     *           roMAX  - store maximum
     *           roMIN  - store minimum
     */
    enum TReductionOperator
    {
      roNONE, roRMS, roMAX, roMIN
    };

    /**
     * @brief   Constructor - there is no sensor mask by default!
     * @details Constructor - there is no sensor mask by default!
     * it links the HDF5 dataset, source (sampled matrix) and the reduction
     * operator together. The constructor DOES NOT allocate memory because the
     * size of the sensor mask is not known at the time the instance of
     * the class is being created.
     *
     * @param [in] HDF5_File           - Handle to the HDF5 (output) file
     * @param [in] HDF5_RootObjectName - The root object that stores the sample
     *                                   data (dataset or group)
     * @param [in] SourceMatrix        - The source matrix (only real matrices
     *                                   are supported)
     * @param [in] ReductionOp         - Reduction operator
     * @param [in] BufferToReuse       - An external buffer can be used to line
     *                                   up the grid points
     */
    TBaseOutputHDF5Stream(Hdf5File &             HDF5_File,
                          const char *             HDF5_RootObjectName,
                          const TRealMatrix &      SourceMatrix,
                          const TReductionOperator ReductionOp,
                          float *                  BufferToReuse = NULL)
            : HDF5_File          (HDF5_File),
              HDF5_RootObjectName(NULL),
              SourceMatrix       (SourceMatrix),
              ReductionOp        (ReductionOp),
              BufferReuse        (BufferToReuse != NULL),
              BufferSize         (0),
              StoreBuffer        (BufferToReuse)
    {
      // copy the dataset name (just for sure)
      this->HDF5_RootObjectName = new char[strlen(HDF5_RootObjectName)];
      strcpy(this->HDF5_RootObjectName, HDF5_RootObjectName);
    };

    /**
     * @brief Destructor.
     * @details Destructor.
     */
    virtual ~TBaseOutputHDF5Stream()
    {
      delete [] HDF5_RootObjectName;
    };

    /// Create a HDF5 stream and allocate data for it.
    virtual void Create() = 0;

    /// Reopen the output stream after restart.
    virtual void Reopen() = 0;

    /// Sample data into buffer, apply reduction or flush to disk - based on a sensor mask.
    virtual void Sample() = 0;

    /// Apply post-processing on the buffer and flush it to the file.
    virtual void PostProcess();

    /// Checkpoint the stream.
    virtual void Checkpoint() = 0;

    /// Close stream (apply post-processing if necessary, flush data and close).
    virtual void Close() = 0;

  protected:
    /// Default constructor not allowed.
    TBaseOutputHDF5Stream();
    /// Copy constructor not allowed.
    TBaseOutputHDF5Stream(const TBaseOutputHDF5Stream & src);
    /// Operator = not allowed (we don't want any data movements).
    TBaseOutputHDF5Stream & operator = (const TBaseOutputHDF5Stream & src);

    /// A generic function to allocate memory - not used in the base class.
    virtual void AllocateMemory();
    /// A generic function to free memory - not used in the base class.
    virtual void FreeMemory();

    /// HDF5 file handle.
    Hdf5File &             HDF5_File;
    /// Dataset name.
    char *                   HDF5_RootObjectName;
    /// Source matrix to be sampled.
    const TRealMatrix&       SourceMatrix;
    /// Reduction operator.
    const TReductionOperator ReductionOp;

    /// if true, the container reuses e.g. Temp_1_RS3D, Temp_2_RS3D, Temp_3_RS3D.
    bool    BufferReuse;
    /// Buffer size.
    size_t  BufferSize;
    /// Temporary buffer for store - only if Buffer Reuse = false!
    float * StoreBuffer;

    /// chunk size of 4MB in number of float elements.
    static const size_t ChunkSize_4MB = 1048576;

    /// The minimum number of elements to start sampling in parallel (4MB).
    static const size_t MinGridpointsToSampleInParallel = 1048576;
};// end of TOutputHDF5Stream
//------------------------------------------------------------------------------


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
                           const TRealMatrix &      SourceMatrix,
                           const TIndexMatrix &     SensorMask,
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
    const TIndexMatrix & SensorMask;
    /// Handle to a HDF5 dataset.
    hid_t  HDF5_DatasetId;

    /// Time step to store (N/A for aggregated).
    size_t SampledTimeStep;
}; // end of TIndexOutputHDF5Stream
//------------------------------------------------------------------------------


/**
 * @class TCuboidOutputHDF5Stream
 * @brief Output stream for quantities sampled by a cuboid corner sensor mask.
 * @details Output stream for quantities sampled by a cuboid corner sensor mask.
 *          This class writes data into separated datasets (one per cuboid) under
 *          a given dataset in the HDF5 file (time-series as well as aggregations).
 *
 */
class TCuboidOutputHDF5Stream : public TBaseOutputHDF5Stream
{
  public:
    /// Constructor - links the HDF5 File, SourceMatrix, and SensorMask together.
    TCuboidOutputHDF5Stream(Hdf5File &             HDF5_File,
                            const char *             HDF5_GroupName,
                            const TRealMatrix &      SourceMatrix,
                            const TIndexMatrix &     SensorMask,
                            const TReductionOperator ReductionOp,
                            float *                  BufferToReuse = NULL);

    /// Destructor.
    virtual ~TCuboidOutputHDF5Stream();

    /// Create a HDF5 stream and allocate data for it.
    virtual void Create();

    /// Reopen the output stream after restart and reload data.
    virtual void Reopen();

    /// Sample data into buffer and apply reduction, or flush to disk - based on a sensor mask.
    virtual void Sample();

    /// Apply post-processing on the buffer and flush it to the file.
    virtual void PostProcess();

    /// Checkpoint the stream and close.
    virtual void Checkpoint();

    /// Close stream (apply post-processing if necessary, flush data and close).
    virtual void Close();

  protected:
    /**
     * @struct TCuboidInfo
     * @brief This structure information about a HDF5 dataset (one cuboid).
     * Namely, its HDF5_ID, Starting position in a lineup buffer.
     */
    struct TCuboidInfo
    {
      /// ID of the dataset storing the given cuboid.
      hid_t  HDF5_CuboidId;
      /// Having a single buffer for all cuboids, where this one starts.
      size_t StartingPossitionInBuffer;
    };

    /// Create a new dataset for a given cuboid specified by index (order).
    virtual hid_t CreateCuboidDataset(const size_t Index);

    /// Flush the buffer to the file.
    virtual void FlushBufferToFile();

    /// Sensor mask to sample data.
    const TIndexMatrix &     SensorMask;

    /// Handle to a HDF5 dataset.
    hid_t                    HDF5_GroupId;

    /// vector keeping handles and positions of all cuboids
    std::vector<TCuboidInfo> CuboidsInfo;

    /// Timestep to store (N/A for aggregated).
    size_t                   SampledTimeStep;

};// end of TCubodiOutputHDF5Stream
//------------------------------------------------------------------------------




/**
 * @class TWholeDomainOutputHDF5Stream
 * @brief Output stream for quantities sampled in the whole domain.
 * @details Output stream for quantities sampled in the whole domain.
 *          The data is stored in a single dataset (aggregated quantities only).
 */
class TWholeDomainOutputHDF5Stream : public TBaseOutputHDF5Stream
{

  public:
    /// Constructor - links the HDF5 File, SourceMatrix, and SensorMask together.
    TWholeDomainOutputHDF5Stream(Hdf5File &             HDF5_File,
                                 const char *             HDF5_DatasetName,
                                 const TRealMatrix &      SourceMatrix,
                                 const TReductionOperator ReductionOp,
                                 float *                  BufferToReuse = NULL);

    /// Destructor.
    virtual ~TWholeDomainOutputHDF5Stream();

    /// Create a HDF5 stream and allocate data for it.
    virtual void Create();

    /// Reopen the output stream after restart and reload data.
    virtual void Reopen();

    /// Sample data into buffer and apply reduction, or flush to disk (no sensor mask here).
    virtual void Sample();

    /// Apply post-processing on the buffer and flush it to the file.
    virtual void PostProcess();

    ///Checkpoint the stream and close.
    virtual void Checkpoint();

    /// Close stream (apply post-processing if necessary, flush data and close).
    virtual void Close();

  protected:
    /// Flush the buffer to the file.
    virtual void FlushBufferToFile();

    /// Handle to a HDF5 dataset.
    hid_t  HDF5_DatasetId;

    /// Time step to store (N/A for aggregated).
    size_t SampledTimeStep;
};// end of TWholeDomainOutputHDF5Stream
//------------------------------------------------------------------------------

#endif	/* OUTPUTREALSTREAM_H */
