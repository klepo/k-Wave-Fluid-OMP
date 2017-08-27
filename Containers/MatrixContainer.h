/**
 * @file        MatrixContainer.h
 * @author      Jiri Jaros \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The header file containing the matrix container.
 *
 * @version     kspaceFirstOrder3D 2.16
 *
 * @date        14 September 2012, 14:33 (created) \n
 *              27 August    2017, 09:29 (revised)
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

#ifndef MATRIX_CONTAINER_H
#define MATRIX_CONTAINER_H

#include <string.h>
#include <map>

#include <MatrixClasses/BaseMatrix.h>
#include <Containers/MatrixRecord.h>

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

#endif	/* MATRIX_CONTAINER_H */
