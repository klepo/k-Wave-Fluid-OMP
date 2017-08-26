/**
 * @file        UXYZ_SGXYZMatrix.h
 * @author      Jiri Jaros
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The header file containing the particle velocity matrix.
 *
 * @version     kspaceFirstOrder3D 2.16
 *
 * @date        28 July      2011, 11:37 (created) \n
 *              25 August    2017, 22:02 (revised)
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


#ifndef UXYZ_SGXYZ_REAL_MATRIX_H
#define UXYZ_SGXYZ_REAL_MATRIX_H


#include <MatrixClasses/RealMatrix.h>
#include <MatrixClasses/FFTWComplexMatrix.h>
#include <MatrixClasses/IndexMatrix.h>

/**
 * @class   VelocityMatrix.
 * @brief   The velocity matrix
 * @details The velocity matrix. This class implements a couple of kernels that modify the particle velocity.
 */
class VelocityMatrix : public RealMatrix
{
  public:

   // Default constructor not allowed for public.
    VelocityMatrix() = delete;
    /**
     * @brief   Constructor.
     * @details Constructor allocating memory.
     * @param [in] dimensionSizes - Dimension sizes
     */
    VelocityMatrix(const DimensionSizes& dimensionSizes) : RealMatrix(dimensionSizes) {};
    /// Copy constructor is not allowed.
    VelocityMatrix(const VelocityMatrix& src) = delete;
    /// Destructor
    virtual ~VelocityMatrix() {};
    /// operator = not allowed for public.
    VelocityMatrix& operator= (const VelocityMatrix& src) = delete;


   /**
    * @brief Compute velocity for the initial pressure problem, heterogeneous medium, uniform grid.
    *
    * <b> Matlab code: </b>
    *
    * \verbatim
        ux_sgx = dt ./ rho0_sgx .* ifft(ux_sgx).
        uy_sgy = dt ./ rho0_sgy .* ifft(uy_sgy).
        uz_sgz = dt ./ rho0_sgz .* ifft(uz_sgz).
     \endverbatim
     *
     * @param [in] dtRho0Sgxyz - Density matrix in x, y or z direction
     * @param [in] fftTemp     - temporary FFT matrix.
     */
   void computeInitialVelocity(const RealMatrix&  dtRho0Sgxyz,
                               FftwComplexMatrix& fftTemp);
   /**
    * @brief Compute velocity for the initial pressure problem, homogeneous medium, uniform grid.
    *
    * <b> Matlab code: </b>
    *
    * \verbatim
        ux_sgx = dt ./ rho0_sgx .* ifft(ux_sgx).
        uy_sgy = dt ./ rho0_sgy .* ifft(uy_sgy).
        uz_sgz = dt ./ rho0_sgz .* ifft(uz_sgz).
    \endverbatim
    *
    * @param [in] dtRho0Sgxyz - Scalar density in x, y or z direction
    * @param [in] fftTemp     - temporary FFT matrix.
    */
   void computeInitialVelocityHomogeneousUniform(const float        dtRho0Sgxyz,
                                                 FftwComplexMatrix& fftTemp);

   /**
    * @brief Compute acoustic velocity for initial pressure problem, homogenous medium, non-uniform grid, x direction.
    *
    * <b> Matlab code: </b>
    *
    * \verbatim
        ux_sgx = dt ./ rho0_sgx .* dxudxn_sgx .* ifft(ux_sgx).
    \endverbatim
    *
    * @param [in] dtRho0Sgx - Scalar density in x direction
    * @param [in] dxudxnSgx - Non uniform grid shift in x direction.
    * @param [in] fftTemp   - temporary FFT matrix.
    */
   void computeInitialVelocityXHomogeneousNonuniform(const float        dtRho0Sgx,
                                                     const RealMatrix&  dxudxnSgx,
                                                     FftwComplexMatrix& fftTemp);
   /**
    * @brief Compute acoustic velocity for initial pressure problem, homogenous medium, non-uniform grid, y direction.
    *
    * <b> Matlab code: </b>
    *
    * \verbatim
        uy_sgy = dt ./ rho0_sgy .* dyudxn_sgy .* ifft(uy_sgy).
    \endverbatim
    *
    * @param [in] dtRho0Sgy - Scalar density in y direction
    * @param [in] dyudynSgy - Non uniform grid shift in y direction.
    * @param [in] fftTemp   - temporary FFT matrix.
    */
   void computeInitialVelocityYHomogeneousNonuniform(const float        dtRho0Sgy,
                                                     const RealMatrix&  dyudynSgy,
                                                     FftwComplexMatrix& fftTemp);
   /**
    * @brief Compute acoustic velocity for initial pressure problem, homogenous medium, non-uniform grid, z direction.
    *
    * <b> Matlab code: </b>
    *
    * \verbatim
        uz_sgz = dt ./ rho0_sgz .* dzudzn_sgz .* ifft(uz_sgz).
    \endverbatim
    *
    * @param [in] dtRho0Sgz - Scalar density in z direction
    * @param [in] dzudznSgz - Non uniform grid shift in z direction.
    * @param [in] fftTemp   - temporary FFT matrix.
    */
   void computeInitialVelocityZHomogeneousNonuniform(const float        dtRho0Sgz,
                                                     const RealMatrix&  dzudznSgz,
                                                     FftwComplexMatrix& fftTemp);

    //------------------------------------------------- X dimension --------------------------------------------------//
    /**
     * @brief Compute acoustic velocity for heterogeneous medium and a uniform grid, x direction.
     *
     * @param [in] ifftX     - ifftn( bsxfun(\@times, ddx_k_shift_pos, kappa .* p_k))
     * @param [in] dtRho0Sgx - Acoustic density on staggered grid in x direction.
     * @param [in] pmlX      - Perfectly matched layer in x direction.
     */
    void computeVelocityX(const RealMatrix& ifftX,
                          const RealMatrix& dtRho0Sgx,
                          const RealMatrix& pmlX);
    /**
     * @brief Compute acoustic velocity for homogeneous medium and a uniform grid, x direction.
     *
     * @param [in] ifftX  - ifftn( bsxfun(\@times, ddx_k_shift_pos, kappa .* p_k))
     * @param [in] dtRho0 - precomputed dt .* rho0 scalar value
     * @param [in] pmlX    - Perfectly matched layer in x direction.
     */
    void computeVelocityXHomogeneousUniform(const RealMatrix& ifftX,
                                            const float       dtRho0,
                                            const RealMatrix& pmlX);
    /**
     * @brief Compute acoustic velocity for homogenous medium and non-uniform grid, x direction.
     *
     * @param [in]     ifftX     - ifftn( bsxfun(\@times, ddx_k_shift_pos, kappa .* p_k))
     * @param [in]     dtRho0    - precomputed dt .* rho0 scalar value
     * @param [in]     dxudxnSgx - Non uniform grid shift in x direction.
     * @param [in]     pmlX      - Perfectly matched layer in x direction.
     */
    void computeVelocityXHomogeneousNonuniform(const RealMatrix& ifftX,
                                               const float       dtRho0,
                                               const RealMatrix& dxudxnSgx,
                                               const RealMatrix& pmlX);

    //------------------------------------------------- Y dimension --------------------------------------------------//
    /**
     * @brief Compute acoustic velocity for heterogeneous medium and a uniform grid, y direction.
     *
     * @param [in] ifftY     - ifftn( bsxfun(\@times, ddy_k_shift_pos, kappa .* p_k))
     * @param [in] dtRho0Sgy - Acoustic density on staggered grid in y direction.
     * @param [in] pmlY      - Perfectly matched layer in y direction.
     */
    void computeVelocityY(const RealMatrix& ifftY,
                          const RealMatrix& dtRho0Sgy,
                          const RealMatrix& pmlY);
    /**
     * @brief Compute acoustic velocity for homogeneous medium and a uniform grid, y direction.
     *
     * @param [in] ifftY  - ifftn( bsxfun(\@times, ddy_k_shift_pos, kappa .* p_k))
     * @param [in] dtRho0 - precomputed dt .* rho0 scalar value
     * @param [in] pmlY    - Perfectly matched layer in y direction.
     */
    void computeVelocityYHomogeneousUniform(const RealMatrix& ifftY,
                                            const float       dtRho0,
                                            const RealMatrix& pmlY);
    /**
     * @brief Compute acoustic velocity for homogenous medium and non-uniform grid, y direction.
     *
     * @param [in]     ifftY     - ifftn( bsxfun(\@times, ddy_k_shift_pos, kappa .* p_k))
     * @param [in]     dtRho0    - precomputed dt .* rho0 scalar value
     * @param [in]     dyudynSgy - Non uniform grid shift in y direction.
     * @param [in]     pmlY      - Perfectly matched layer in y direction.
     */
    void computeVelocityYHomogeneousNonuniform(const RealMatrix& ifftY,
                                               const float       dtRho0,
                                               const RealMatrix& dyudynSgy,
                                               const RealMatrix& pmlY);

    //------------------------------------------------- Z dimension --------------------------------------------------//
    /**
     * @brief Compute acoustic velocity for heterogeneous medium and a uniform grid, z direction.
     * @param [in] ifftZ     - ifftn( bsxfun(\@times, ddz_k_shift_pos, kappa .* p_k))
     * @param [in] dtRho0Sgz - Acoustic density on staggered grid in z direction.
     * @param [in] pmlZ      - Perfectly matched layer in z direction.
     */
    void computeVelocityZ(const RealMatrix& ifftZ,
                          const RealMatrix& dtRho0Sgz,
                          const RealMatrix& pmlZ);
    /**
     * @brief Compute acoustic velocity for homogeneous medium and a uniform grid, z direction
     *
     * @param [in] ifftZ  - ifftn( bsxfun(\@times, ddz_k_shift_pos, kappa .* p_k))
     * @param [in] dtRho0 - precomputed dt .* rho0 scalar value
     * @param [in] pmlZ    - Perfectly matched layer in z direction.
     */
    void computeVelocityZHomogeneousUniform(const RealMatrix& ifftZ,
                                            const float       dtRho0,
                                            const RealMatrix& pmlZ);
    /**
     * @brief Compute acoustic velocity for homogenous medium and non-uniform grid, z direction.
     *
     * @param [in]     ifftZ     - ifftn( bsxfun(\@times, ddz_k_shift_pos, kappa .* p_k))
     * @param [in]     dtRho0    - precomputed dt .* rho0 scalar value
     * @param [in]     dzudznSgz - Non uniform grid shift in z direction.
     * @param [in]     pmlZ      - Perfectly matched layer in z direction.
     */
    void computeVelocityZHomogeneousNonuniform(const RealMatrix& ifftZ,
                                               const float       dtRho0,
                                               const RealMatrix& dzudznSgz,
                                               const RealMatrix& pmlZ);

    //--------------------------------------------------- Sources ----------------------------------------------------//
    /**
     * @brief Add transducer data source to velocity x component.
     *
     * @param [in] velocitySourceIndex   - Where to add the signal (source geometry).
     * @param [in] transducerSourceInput - Transducer signal.
     * @param [in] delayMask             - Delay mask to push the signal in the domain (incremented per invocation).
     * @param [in] timeIndex             - Actual time step.
     */
     void addTransducerSource(const IndexMatrix& velocitySourceIndex,
                              const RealMatrix&  transducerSourceInput,
                              const IndexMatrix& delayMask,
                              const size_t       timeIndex);

     /**
      * @brief Add in velocity source terms.
      * @param [in] velocitySourceInput - Source input to add.
      * @param [in] velocitySourceIndex - Source geometry index matrix.
      * @param [in] timeIndex           - Actual time step.
      * @param [in] velocitySourceMode  - Velocity source mode (0 = Dirichlet boundary, 1 = add in).
      * @param [in] velocitySourceMany  - Velocity source mode (0 = One series, 1 = multiple series).
      */
     void addVelocitySource(const RealMatrix & velocitySourceInput,
                            const IndexMatrix& velocitySourceIndex,
                            const size_t       timeIndex,
                            const size_t       velocitySourceMode,
                            const size_t       velocitySourceMany);

  protected:

  private:

}; // end of VelocityMatrix
//----------------------------------------------------------------------------------------------------------------------

#endif	/* UXYZ_SGXYZ_REAL_MATRIX_H */

