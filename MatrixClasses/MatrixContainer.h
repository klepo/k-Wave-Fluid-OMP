/**
 * @file        MatrixContainer.h
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au
 *
 * @brief       The header file containing the matrix container
 *
 * @version     kspaceFirstOrder3D 2.14
 * @date        14 September 2012, 14:33 (created) \n
 *              08 July      2014, 13:57 (revised)
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

#ifndef MATRIXCONTAINER_H
#define	MATRIXCONTAINER_H

#include <string.h>
#include <map>

#include <MatrixClasses/BaseMatrix.h>
#include <MatrixClasses/BaseFloatMatrix.h>
#include <MatrixClasses/RealMatrix.h>
#include <MatrixClasses/ComplexMatrix.h>
#include <MatrixClasses/FFTWComplexMatrix.h>
#include <MatrixClasses/UXYZ_SGXYZMatrix.h>
#include <MatrixClasses/LongMatrix.h>
#include <MatrixClasses/OutputHDF5Stream.h>

#include <Utils/MatrixNames.h>
#include <Utils/DimensionSizes.h>

/**
 * @enum TMatrixID
 * @brief Matrix identifers of all matrices in the k-space code
 */
enum TMatrixID
{
    kappa, c2, p,

    ux_sgx,uy_sgy, uz_sgz,
    duxdx, duydy, duzdz,
    dxudxn    , dyudyn    , dzudzn,
    dxudxn_sgx, dyudyn_sgy, dzudzn_sgz,

    rhox, rhoy, rhoz, rho0,
    dt_rho0_sgx, dt_rho0_sgy, dt_rho0_sgz,


    p0_source_input, sensor_mask_index, sensor_mask_corners,
    ddx_k_shift_pos, ddy_k_shift_pos, ddz_k_shift_pos,
    ddx_k_shift_neg, ddy_k_shift_neg, ddz_k_shift_neg,

    pml_x_sgx, pml_y_sgy, pml_z_sgz,
    pml_x    , pml_y    , pml_z,

    absorb_tau, absorb_eta, absorb_nabla1, absorb_nabla2, BonA,

    ux_source_input, uy_source_input, uz_source_input,
    p_source_input,

    u_source_index, p_source_index, transducer_source_input,
    delay_mask,

    //---------------- output matrices -------------//
    p_sensor_raw,  p_sensor_rms, p_sensor_max, p_sensor_min,
    p_sensor_max_all, p_sensor_min_all, p_sensor_i_1_raw,
    ux_sensor_raw, uy_sensor_raw, uz_sensor_raw,
    ux_sensor_i_1_agr_2, uy_sensor_i_1_agr_2, uz_sensor_i_1_agr_2,// aggregated values form two points

    ux_sensor_rms, uy_sensor_rms, uz_sensor_rms,
    ux_sensor_max, uy_sensor_max, uz_sensor_max,
    ux_sensor_min, uy_sensor_min, uz_sensor_min,
    ux_sensor_max_all, uy_sensor_max_all, uz_sensor_max_all,
    ux_sensor_min_all, uy_sensor_min_all, uz_sensor_min_all,

    Ix_sensor_avg, Iy_sensor_avg, Iz_sensor_avg,
    Ix_sensor_max, Iy_sensor_max, Iz_sensor_max,


    //--------------Temporary matrices -------------//
    Temp_1_RS3D, Temp_2_RS3D, Temp_3_RS3D,
    FFT_X_temp, FFT_Y_temp, FFT_Z_temp
}; // enum TMatrixID
//------------------------------------------------------------------------------


/**
 * @struct TMatrixRecord
 * @brief  A structure storing details about the matrix. The matrix container
 * stores this structures.
 */
struct TMatrixRecord
{
    /**
     * @enum TMatrixDataType
     * @brief All possible types of the matrix
     */
    enum TMatrixDataType {mdtReal, mdtComplex, mdtIndex, mdtFFTW, mdtUxyz};

    /// Pointer to the matrix object
    TBaseMatrix   * MatrixPtr;
    /// Matrix data type
    TMatrixDataType MatrixDataType;
    /// Matrix dimension sizes
    TDimensionSizes DimensionSizes;
    /// Is the matrix content loaded from the HDF5 file
    bool            LoadData;
    /// Is the matrix necessary to be preserver when checkpoint is enabled
    bool            Checkpoint;
    /// HDF5 matrix name
    string          HDF5MatrixName;

    /// Default constructor
    TMatrixRecord() : MatrixPtr(NULL), MatrixDataType(mdtReal),
                      DimensionSizes(), LoadData(false), Checkpoint(false),
                      HDF5MatrixName("")
                      {};

    /// Copy constructor
    TMatrixRecord(const TMatrixRecord& src);

    /// operator =
    TMatrixRecord& operator = (const TMatrixRecord& src);

    /// Set all values of the record
    void SetAllValues(TBaseMatrix *          MatrixPtr,
                      const TMatrixDataType  MatrixDataType,
                      const TDimensionSizes  DimensionSizes,
                      const bool             LoadData,
                      const bool             Checkpoint,
                      const string           HDF5MatrixName);

    virtual ~TMatrixRecord() {};


};// end of TMatrixRecord
//------------------------------------------------------------------------------



/**
 * @class TMatrixContainer
 * @brief Class implementing the matrix container
 */
class TMatrixContainer
{
  public:

    /// Constructor
    TMatrixContainer() {}
    /// Destructor
    virtual ~TMatrixContainer();

    /**
     * Get number of matrices in the container
     * @return number of matrices in the container
     */
    size_t size ()  {return MatrixContainer.size(); };
    /**
     * Is the container empty?
     * @return true if the container is empty
     */
    bool   empty()  {return MatrixContainer.empty();};


    /// Create instances of all objects in the container
    void CreateAllObjects();
    /// Load all matrices from the HDF5 file
    void LoadDataFromInputHDF5File(THDF5_File & HDF5_File);
    /// Load all matrices from the HDF5 file
    void LoadDataFromCheckpointHDF5File(THDF5_File & HDF5_File);
    /// Store selected matrices into the checkpoint file
    void StoreDataIntoCheckpointHDF5File(THDF5_File & HDF5_File);

    /// Free all matrices - destroy them
    void FreeAllMatrices();


    /// Set all matrices recored - populate the container
    void AddMatricesIntoContainer();

    /**
     * Get matrix record
     * @param [in] MatrixID - Matrix identifier
     * @return the matrix record
     */
    TMatrixRecord& GetMatrixRecord(const TMatrixID MatrixID)
    {
     return MatrixContainer[MatrixID];
    };

    /**
     * operator []
     * @param [in]  MatrixID - Matrix identifier
     * @return the matrix record
     */
    TMatrixRecord& operator [] (const TMatrixID MatrixID)
    {
      return MatrixContainer[MatrixID];
    };


    /**
     * Get BaseMatrix from the container
     * @param [in] MatrixID - Matrix identifier
     * @return Base Matrix
     */
    TBaseMatrix& GetBaseMatrix(const TMatrixID MatrixID)
    {
      return static_cast<TBaseMatrix &> (*(MatrixContainer[MatrixID].MatrixPtr));
    };

    /**
     * Get BaseFloatMatrix from the container
     * @param [in] MatrixID - Matrix identifier
     * @return BaseFloatMatrix
     */
    TBaseFloatMatrix& GetBaseFloatMatrix(const TMatrixID MatrixID)
    {
      return static_cast<TBaseFloatMatrix &> (*(MatrixContainer[MatrixID].MatrixPtr));
    };

    /**
     * Get RealMatrix from the container
     * @param [in] MatrixID - Matrix identifier
     * @return RealMatrix
     */
    TRealMatrix& GetRealMatrix(const TMatrixID MatrixID)
    {
      return static_cast<TRealMatrix &> (*(MatrixContainer[MatrixID].MatrixPtr));
    };

    /**
     * Get Uxyz_sgzMatrix from the container
     * @param [in] MatrixID - Matrix identifier
     * @return  Uxyz_sgzMatrix
     */
    Tuxyz_sgxyzMatrix& GetUxyz_sgxyzMatrix(const TMatrixID MatrixID)
    {
      return static_cast<Tuxyz_sgxyzMatrix &>(*(MatrixContainer[MatrixID].MatrixPtr));
    };

    /**
     * Get ComplexMatrix from the container
     * @param [in] MatrixID - Matrix identifier
     * @return ComplexMatrix
     */
    TComplexMatrix& GetComplexMatrix(const TMatrixID MatrixID)
    {
      return static_cast<TComplexMatrix &>(*(MatrixContainer[MatrixID].MatrixPtr));
    };

    /**
     * GetFFTWComplexMatrix from the container
     * @param [in] MatrixID - Matrix identifier
     * @return FFTWComplexMatrix
     */
    TFFTWComplexMatrix& GetFFTWComplexMatrix(const TMatrixID MatrixID)
    {
      return static_cast<TFFTWComplexMatrix &>  (*(MatrixContainer[MatrixID].MatrixPtr));
    };

    /**
     * Get LongMatrix matrix from the container
     * @param [in] MatrixID - Matrix identifier
     * @return LongMatrix
     */
    TLongMatrix& GetLongMatrix(const TMatrixID MatrixID)
    {
      return static_cast<TLongMatrix &>(*(MatrixContainer[MatrixID].MatrixPtr));
    };


  protected:


  private:

    /// Datatype for map associating the matrix ID enum and matrix record
    typedef map<TMatrixID, TMatrixRecord> TMatrixRecordContainer;

    /// Map holding the container
    TMatrixRecordContainer MatrixContainer;

    /// Copy constructor is not allowed for public
    TMatrixContainer(const TMatrixContainer& src);

    /// Operator = is not allowed for public
    TMatrixContainer & operator = (const TMatrixContainer& src);

    /// Print error and throw an exception
    void PrintErrorAndThrowException(const char * FMT,
                                     const string HDF5MatrixName,
                                     const char * File,
                                     const int Line);

};// end of TMatrixContainer
//------------------------------------------------------------------------------


/**
 * @class TOutputStreamContainer
 * @brief A container for output streams
 */
class TOutputStreamContainer
{
  public:
    /// Constructor
    TOutputStreamContainer() {};
    /// Destructor
    virtual ~TOutputStreamContainer();

    /// Get size of the container
    size_t size() const
    {
      return OutputStreamContainer.size();
    };

    /// Is the container empty?
    bool empty() const
    {
      return OutputStreamContainer.empty();
    };

    /**
     * @brief Operator []
     * @details Operator []
     * @param MatrixID
     * @return
     */
    TBaseOutputHDF5Stream & operator [] (const TMatrixID MatrixID)
    {
      return (* (OutputStreamContainer[MatrixID]));
    };

    /// Create all streams
    void AddStreamsIntoContainer(TMatrixContainer & MatrixContainer);

    /// Create all streams
    void CreateStreams();
    /// Sample all streams
    void SampleStreams();

    /// Close all streams
    void CloseStreams();

    /// Free all streams - destroy them
    void FreeAllStreams();

  protected:
    // Create a new output stream
    TBaseOutputHDF5Stream * CreateNewOutputStream(TMatrixContainer &       MatrixContainer,
                                                  const TMatrixID          SampledMatrixID,
                                                  const char *             HDF5_DatasetName,
                                                  const TBaseOutputHDF5Stream::TReductionOperator ReductionOp,
                                                  float *                  BufferToReuse = NULL);


    /// Copy constructor not allowed for public
    TOutputStreamContainer(const TOutputStreamContainer &);
    /// Operator = not allowed for public
    TOutputStreamContainer & operator = (TOutputStreamContainer &);

  private:
    /// Output stream map
    typedef map < TMatrixID, TBaseOutputHDF5Stream * > TOutputStreamMap;
    /// Map with output streams
    TOutputStreamMap OutputStreamContainer;

}; // end of TOutputStreamContainer
//------------------------------------------------------------------------------

#endif	/* MATRIXCONTAINER_H */

