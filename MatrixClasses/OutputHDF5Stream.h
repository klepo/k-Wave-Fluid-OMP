/**
 * @file        OutputHDF5Stream.h
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au   \n
 * 
 * @brief       The header file of the class saving RealMatrix data into 
 *              the output HDF5 file
 * 
 * @version     kspaceFirstOrder3D 2.14
 * @date        11 July 2012, 10:30 (created) \n
 *              25 February 2014, 16:40 (revised)
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
 * @class TOutputHDF5Stream
 * @brief Output stream for sensor data. 
 *        This class write data to HDF5 file (time-series as well as aggregations)
 */
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

    /// Create stream for sensor mask based date
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

