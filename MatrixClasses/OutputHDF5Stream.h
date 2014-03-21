/**
 * @file        OutputHDF5Stream.h
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au   \n
 * 
 * @brief       The header file of classes responsible for storing output 
 *              quantities into the output HDF5 file
 * 
 * @version     kspaceFirstOrder3D 2.14
 * @date        11 July  2012, 10:30 (created) \n
 *              18 March 2014, 15:00 (revised)
 * 
 * @section License
 * This file is part of the C++ extension of the k-Wave Toolbox (http://www.k-wave.org).\n
 * Copyright (C) 2012 Jiri Jaros and Bradley Treeby
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

#include <MatrixClasses/RealMatrix.h>
#include <MatrixClasses/LongMatrix.h>

#include <HDF5/HDF5_File.h>

using namespace std;

/**
 * @class TBaseOutputHDF5Stream
 * @brief Base class for output data streams (sampled data)
 *       
 *        
 */
class TBaseOutputHDF5Stream
{
 public:
   
    /**
     * @enum TOutputHDF5StreamReductionOperator 
     * @brief How to aggregate data \n
     *        roNONE - store actual data
     *        roRMS  - calculate root mean square \n
     *        roMAX  - store maximum
     *        roMIN  - store minimum     
     */
    enum TReductionOperator 
    {
      roNONE, roRMS, roMAX, roMIN, 
    };
  
    // Constructor (by default there is no sensor mask)
    TBaseOutputHDF5Stream(THDF5_File &             HDF5_File,
                          const char *             HDF5ObjectName,
                          const TRealMatrix &      SourceMatrix, 
                          const TReductionOperator ReductionOp,
                          float *                  BufferToReuse = NULL) 
            : HDF5_File     (HDF5_File), 
              HDF5ObjectName(HDF5ObjectName),
              SourceMatrix  (SourceMatrix),
              ReductionOp   (ReductionOp),
              BufferReuse   (BufferToReuse == NULL),
              StoringBuffer (BufferToReuse)
    {};
                                                        
    /// Destructor
    virtual ~TBaseOutputHDF5Stream() = 0;
    
    
    /// Create a HDF5 stream and allocate data for it - a virtual method
    virtual void Create(const size_t NumberOfSampledElementsPerStep) = 0;
    
    /// Reopen the output stream
    /// @TODO - Will be implemented with checkpoint-restart
    virtual void Reopen(const size_t NumberOfSampledElementsPerStep) = 0;
                       
    /// Sample data into buffer, apply reduction or flush to disk - based on a sensor mask
    virtual void Sample() = 0;
                
    /// Checkpoint the stram and close
    virtual void Checkpoint() = 0;
    
    /// Close stream (apply post-processing if necessary, flush data and close)
    virtual void Close() = 0;
    
 protected:
    /// Default constructor not allowed
    TBaseOutputHDF5Stream();
    /// Copy constructor not allowed
    TBaseOutputHDF5Stream(const TBaseOutputHDF5Stream & src);
    /// Operator = not allowed (we don't want any data movements)
    TBaseOutputHDF5Stream & operator = (const TBaseOutputHDF5Stream & src);
  
 
    /// HDF5 file handle
    const THDF5_File &       HDF5_File;            
    /// Dataset name
    const char *             HDF5ObjectName;
    /// Source matrix to be sampled
    const TRealMatrix&       SourceMatrix;    
    /// Reduction operator    
    const TReductionOperator ReductionOp;
        
    /// if true, the container reuses e.g. Temp_1_RS3D, Temp_2_RS3D, Temp_3_RS3D
    bool BufferReuse;
    /// Temporary buffer for store - only if Buffer Reuse = false!
    float * StoringBuffer;    
        
    /// chunk size of 4MB in number of float elements
    static const size_t ChunkSize_4MB = 1048576; 
      
};// end of TOutputHDF5Stream
//------------------------------------------------------------------------------








/**
 * @class TIndexOutputHDF5Stream
 * @brief Output stream for quantities sampled by an index sensor mask.
 *        This class writes data to a single dataset in a root group of the HDF5 
 *        file (time-series as well as aggregations)
 * 
 */
class TIndexOutputHDF5Stream
{
  
};// end of TIndexOutputHDF5Stream
//------------------------------------------------------------------------------


/**
 * @class TCubodidOutputHDF5Stream
 * @brief Output stream for quantities sampled by a cuboid sensor mask.
 *        This class writes data into separated datasets (one per cuboid) under 
 *        a given dataset in the HDF5 file (time-series as well as aggregations)
 * 
 */
class TCubodidOutputHDF5Stream
{
  
};// end of TCubodiOutputHDF5Stream
//------------------------------------------------------------------------------

/**
 * @class TWholeDomainOutputHDF5Stream 
 * @brief Output stream for quantities sampled in the whole domain.
 *        The data is stored in a single dataset (aggregated quantities only)
 */
class TWholeDomainOutputHDF5Stream
{
  
};// end of TWholeDomainOutputHDF5Stream
//------------------------------------------------------------------------------
















class TOutputHDF5Stream
{
  public:
    
    /**
     * @enum TOutputHDF5StreamReductionOperator 
     * @brief How to aggregate data
     */
    enum TReductionOperator 
    {
      roRMS, roMAX, roMIN
    };
    
    /// Constructor
    TOutputHDF5Stream() : HDF5_File(NULL), HDF5_Dataset_id(H5I_BADID), 
                          Position(0,0,0), BufferSize(0,0,0), TotalDatasetSize(0,0,0),  
                          BufferReuse(false), StoringBuffer (NULL), 
                          DatasetName(NULL),  WholeDomainStream(false) 
    {
    };

    /// Create stream for a sensor mask based data
    virtual void CreateStream(THDF5_File & HDF5_File, 
                              const char * DatasetName,
                              const size_t NumberOfSensorElements, 
                              const size_t NumberOfTimeSteps,
                              float *      BufferToReuse = NULL);

    /// Create stream operating on the whole domain (no sensor_mask), only for aggregated quantities
    virtual void CreateWholeDomainAggregatedStream(THDF5_File & HDF5_File,
                                                   const char * DatasetName);

    /// Reopen the output stream
    /// @TODO - Will be implemented with checkpoint-restart
    virtual void ReopenStream(THDF5_File & HDF5_File,
                              const char * DatasetName,
                              const size_t NumberOfSensorElements,
                              const long   ActTimeStep = 0)
    {};
    
    /// Reopen the AllDomainAggregatedStream
    /// @TODO - Will be implemented with checkpoint-restart
    virtual void ReopenAllDomainAggregatedStream(THDF5_File & HDF5_File, 
                                                 const char * DatasetName,
                                                 const TDimensionSizes &RealMatrixDimensionSizes)
    {};
    
            
    /// Close stream
    virtual void CloseStream();
        
    /// Write data directly to the stream
    virtual void WriteDataDirectly( const TRealMatrix& SourceMatrix, 
                                    const TLongMatrix& SensorMask);
    
    /// Add data into reduce buffer (rms, max, min)    
    virtual void ReduceDataIntoBuffer(const TRealMatrix& SourceMatrix,
                                      const TLongMatrix& SensorMask,
                                      const TReductionOperator  ReductionOp);

    /// Add data into reduce buffer (rms, max, min), all domain    
    virtual void ReduceDataIntoBuffer(const TRealMatrix& SourceMatrix,                                                                             
                                      const TReductionOperator ReductionOp);
    
    
    /// RMS post-processing
    virtual void RMSPostprocessing(const float ScalingCoeff);
        
    /// Write reduced buffer (after post-processing)
    virtual void WriteRudcedBuffer();

    /// Reload ReduceBuffer
    /// @TODO - Will be implemented with checkpoint-restart
    virtual void ReloadRudcedBuffer()
    {};
           
    /// Destructor
    virtual ~TOutputHDF5Stream();
    
  protected:
    /**
     *  Copy constructor is not allowed for public
     * @param src
     */
    TOutputHDF5Stream(const TOutputHDF5Stream& src);

    /// operator = is not allowed for public
    TOutputHDF5Stream & operator =(const TOutputHDF5Stream& src);

    /// Allocate memory
    virtual void AllocateMemory();
    /// Free allocated memory
    virtual void FreeMemory();

    /// Write buffer into file
    virtual void WriteBuffer(const float * Buffer);

    /// HDF5 file handle
    THDF5_File* HDF5_File;
    /// HDF5 dataset handle
    hid_t HDF5_Dataset_id;

    /// Position in the dataset
    TDimensionSizes Position;
    /// Buffer size
    TDimensionSizes BufferSize;
    /// Total HDF5 dataset size
    TDimensionSizes TotalDatasetSize;    
    
    /// if true, the container reuses one of Temp_1_RS3D, Temp_2_RS3D, Temp_3_RS3D
    bool BufferReuse;
    /// Temporary buffer for store - only if Buffer Reuse = false!
    float * StoringBuffer;    
    /// Dataset name
    char  * DatasetName;
    /// All domain stream? - 
    bool WholeDomainStream;
    

    /// chunk size of 4MB in number of float elements
    static const size_t ChunkSize_4MB = 1048576; 
}; // end of TOutputHDF5Stream
//------------------------------------------------------------------------------

#endif	/* OUTPUTREALSTREAM_H */

