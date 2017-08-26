/**
 * @file        MatrixContainer.h
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The header file containing the matrix container.
 *
 * @version     kspaceFirstOrder3D 2.16
 *
 * @date        14 September 2012, 14:33 (created) \n
 *              26 August    2017, 22:48 (revised)
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

#ifndef MATRIXCONTAINER_H
#define	MATRIXCONTAINER_H

#include <string.h>
#include <map>

#include <MatrixClasses/BaseMatrix.h>
#include <MatrixClasses/BaseFloatMatrix.h>
#include <MatrixClasses/RealMatrix.h>
#include <MatrixClasses/ComplexMatrix.h>
#include <MatrixClasses/FftwComplexMatrix.h>
#include <MatrixClasses/VelocityMatrix.h>
#include <MatrixClasses/IndexMatrix.h>

#include <OutputStreams/BaseOutputStream.h>
#include <OutputStreams/IndexOutputStream.h>
#include <OutputStreams/CuboidOutputStream.h>
#include <OutputStreams/WholeDomainOutputStream.h>

#include <Utils/MatrixNames.h>
#include <Utils/DimensionSizes.h>

using namespace std;

/**
 * @enum TMatrixID
 * @brief Matrix identifers of all matrices in the k-space code
 */
enum TMatrixID
{
    kappa, c2, p,

    ux_sgx,uy_sgy, uz_sgz,
    ux_shifted, uy_shifted, uz_shifted,
    duxdx, duydy, duzdz,
    dxudxn    , dyudyn    , dzudzn,
    dxudxn_sgx, dyudyn_sgy, dzudzn_sgz,

    rhox, rhoy, rhoz, rho0,
    dt_rho0_sgx, dt_rho0_sgy, dt_rho0_sgz,

    p0_source_input, sensor_mask_index, sensor_mask_corners,
    ddx_k_shift_pos, ddy_k_shift_pos, ddz_k_shift_pos,
    ddx_k_shift_neg, ddy_k_shift_neg, ddz_k_shift_neg,
    x_shift_neg_r, y_shift_neg_r, z_shift_neg_r,
    pml_x_sgx, pml_y_sgy, pml_z_sgz,
    pml_x    , pml_y    , pml_z,

    absorb_tau, absorb_eta, absorb_nabla1, absorb_nabla2, BonA,

    ux_source_input, uy_source_input, uz_source_input,
    p_source_input,

    u_source_index, p_source_index, transducer_source_input,
    delay_mask,

    //---------------- output matrices -------------//
    p_sensor_raw,  p_sensor_rms, p_sensor_max, p_sensor_min,
    p_sensor_max_all, p_sensor_min_all,
    ux_sensor_raw, uy_sensor_raw, uz_sensor_raw,

    ux_shifted_sensor_raw, uy_shifted_sensor_raw, uz_shifted_sensor_raw, //non_staggered
    ux_sensor_rms, uy_sensor_rms, uz_sensor_rms,
    ux_sensor_max, uy_sensor_max, uz_sensor_max,
    ux_sensor_min, uy_sensor_min, uz_sensor_min,
    ux_sensor_max_all, uy_sensor_max_all, uz_sensor_max_all,
    ux_sensor_min_all, uy_sensor_min_all, uz_sensor_min_all,


    //--------------Temporary matrices -------------//
    Temp_1_RS3D, Temp_2_RS3D, Temp_3_RS3D,
    FFT_X_temp, FFT_Y_temp, FFT_Z_temp, FFT_shift_temp
}; // enum TMatrixID
//------------------------------------------------------------------------------


/**
 * @struct TMatrixRecord
 * @brief  A structure storing details about the matrix. The matrix container
 * stores this structures.
 * @details A structure storing details about the matrix. The matrix container
 * stores the list of these records with the data.
 */
struct TMatrixRecord
{
  /**
   * @enum TMatrixDataType
   * @brief All possible types of the matrix.
   */
  enum TMatrixDataType {mdtReal, mdtComplex, mdtIndex, mdtFFTW, mdtUxyz};

  /// Pointer to the matrix object.
  BaseMatrix   * MatrixPtr;
  /// Matrix data type.
  TMatrixDataType MatrixDataType;
  /// Matrix dimension sizes.
  DimensionSizes dimensionSizes;
  /// Is the matrix content loaded from the HDF5 file.
  bool            LoadData;
  /// Is the matrix necessary to be preserver when checkpoint is enabled.
  bool            Checkpoint;
  /// HDF5 matrix name.
  string          HDF5MatrixName;

  /// Default constructor.
  TMatrixRecord() : MatrixPtr(NULL), MatrixDataType(mdtReal),
          dimensionSizes(), LoadData(false), Checkpoint(false),
          HDF5MatrixName("")
  {};

  /// Copy constructor.
  TMatrixRecord(const TMatrixRecord& src);

  /// operator =
  TMatrixRecord& operator = (const TMatrixRecord& src);

  /// Set all values of the record.
  void SetAllValues(BaseMatrix *          MatrixPtr,
                    const TMatrixDataType  MatrixDataType,
                    const DimensionSizes  dimensionSizes,
                    const bool             LoadData,
                    const bool             Checkpoint,
                    const string           HDF5MatrixName);

  // Destructor.
  virtual ~TMatrixRecord() {};
};// end of TMatrixRecord
//------------------------------------------------------------------------------



/**
 * @class TMatrixContainer
 * @brief Class implementing the matrix container.
 * @details This container is responsible to maintain all the matrices in the
 *          code except the output streams. The matrices are allocated, freed, loaded
 *          stored and checkpointed from here.
 */
class TMatrixContainer
{
  public:

    /// Constructor.
    TMatrixContainer() {}
    /// Destructor.
    virtual ~TMatrixContainer();

    /**
     * @brief Get number of matrices in the container.
     * @details Get number of matrices in the container.
     * @return number of matrices in the container.
     */
    size_t size() const
    {
      return MatrixContainer.size();
    };

    /**
     * @brief Is the container empty?
     * @details Is the container empty?
     * @return true if the container is empty
     */
    bool empty() const
    {
      return MatrixContainer.empty();
    };

    /// Create instances of all objects in the container.
    void CreateAllObjects();

    /// Load all matrices from the HDF5 file.
    void LoadDataFromInputHDF5File(Hdf5File & HDF5_File);
    /// Load all matrices from the HDF5 file.
    void LoadDataFromCheckpointHDF5File(Hdf5File & HDF5_File);
    /// Store selected matrices into the checkpoint file.
    void StoreDataIntoCheckpointHDF5File(Hdf5File & HDF5_File);

    /// Free all matrices - destroy them.
    void FreeAllMatrices();

    /// Set all matrices recored - populate the container.
    void AddMatricesIntoContainer();

    /**
     * @brief Get matrix record (data and information).
     * @details Get matrix record (data and information).
     * @param [in] MatrixID - Matrix identifier
     * @return the matrix record.
     */
    TMatrixRecord& GetMatrixRecord(const TMatrixID MatrixID)
    {
      return MatrixContainer[MatrixID];
    };

    /**
     * @brief operator [].
     * @details operator [].
     * @param [in]  MatrixID - Matrix identifier
     * @return the matrix record.
     */
    TMatrixRecord& operator [] (const TMatrixID MatrixID)
    {
      return MatrixContainer[MatrixID];
    };

    /**
     * @brief Get the matrix with a specific type from the container.
     * @details This template routine returns the reference to the matrix recasted to
     * the specific class.
     * @param [in] MatrixID - Matrix identifier
     * @return Base Matrix
     */
    template <typename T>
    inline T& GetMatrix(const TMatrixID MatrixID)
    {
      return static_cast<T &> (*(MatrixContainer[MatrixID].MatrixPtr));
    };

  protected:


  private:

    /// Datatype for map associating the matrix ID enum and matrix record.
    typedef map<TMatrixID, TMatrixRecord> TMatrixRecordContainer;

    /// Map holding the container.
    TMatrixRecordContainer MatrixContainer;

    /// Copy constructor is not allowed for public.
    TMatrixContainer(const TMatrixContainer& src);

    /// Operator = is not allowed for public.
    TMatrixContainer & operator = (const TMatrixContainer& src);

    /// Print error and throw an exception.
    void PrintErrorAndThrowException(const char * FMT,
                                     const string HDF5MatrixName,
                                     const char * File,
                                     const int Line);

};// end of TMatrixContainer
//------------------------------------------------------------------------------


/**
 * @class TOutputStreamContainer
 * @brief A container for output streams.
 * @details The output stream container maintains matrices used for sampling data.
 * These may or may not require some scratch place or reuse temp matrices.
 */
class TOutputStreamContainer
{
  public:
    /// Constructor.
    TOutputStreamContainer() {};
    /// Destructor.
    virtual ~TOutputStreamContainer();

    /**
     * @brief Get size of the container.
     * @details Get size of the container.
     */
    size_t size() const
    {
      return OutputStreamContainer.size();
    };

    /**
     * @brief  Is the container empty?
     * @details  Is the container empty?
     */
    bool empty() const
    {
      return OutputStreamContainer.empty();
    };

    /**
     * @brief Operator [].
     * @details Operator [].
     * @param [in] MatrixID
     * @return Ouptut stream
     */
    BaseOutputStream & operator [] (const TMatrixID MatrixID)
    {
      return (* (OutputStreamContainer[MatrixID]));
    };

    /// Create all streams in container (no file manipulation).
    void AddStreamsIntoContainer(TMatrixContainer & MatrixContainer);

    /// Create all streams - opens the datasets.
    void CreateStreams();
    /// Reopen streams after checkpoint file (datasets).
    void ReopenStreams();

    /// Sample all streams.
    void SampleStreams();
    /// Post-process all streams and flush them to the file.
    void PostProcessStreams();
    /// Checkpoint streams.
    void CheckpointStreams();

    /// Close all streams.
    void CloseStreams();

    /// Free all streams - destroy them.
    void FreeAllStreams();

  protected:
    /// Create a new output stream.
    BaseOutputStream * CreateNewOutputStream(TMatrixContainer & MatrixContainer,
                                                  const TMatrixID    SampledMatrixID,
                                                  const char *       HDF5_DatasetName,
                                                  const BaseOutputStream::ReduceOperator ReductionOp,
                                                  float *            BufferToReuse = NULL);

    /// Copy constructor not allowed for public.
    TOutputStreamContainer(const TOutputStreamContainer &);
    /// Operator = not allowed for public.
    TOutputStreamContainer & operator = (TOutputStreamContainer &);

  private:
    /// Output stream map.
    typedef map < TMatrixID, BaseOutputStream * > TOutputStreamMap;
    /// Map with output streams.
    TOutputStreamMap OutputStreamContainer;

}; // end of TOutputStreamContainer
//------------------------------------------------------------------------------

#endif	/* MATRIXCONTAINER_H */
