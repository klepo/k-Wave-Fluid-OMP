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
 *              22 August    2017, 13:17 (revised)
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


#ifndef UXYZ_SGXYZREALMATRIX_H
#define	UXYZ_SGXYZREALMATRIX_H


#include <MatrixClasses/RealMatrix.h>
#include <MatrixClasses/FFTWComplexMatrix.h>
#include <MatrixClasses/IndexMatrix.h>

/**
 * @class Tuxyz_sgxyzMatrix.
 * @brief The velocity matrix
 * @details The velocity matrix. This class implements a couple of kernels that
 *          modify the particle velocity.
 */
class Tuxyz_sgxyzMatrix : public TRealMatrix
{
  public:

    /**
     * @brief Constructor.
     * @details Constructor allocating memory.
     * @param [in] DimensionSizes
     */
    Tuxyz_sgxyzMatrix(struct DimensionSizes DimensionSizes) :
            TRealMatrix(DimensionSizes)
    {};


   /// Compute dt ./ rho0_sgx .* ifft (FFT).
   void Compute_dt_rho_sg_mul_ifft_div_2(const TRealMatrix & dt_rho_0_sgx,
                                         TFFTWComplexMatrix& FFT);
   /// Compute dt ./ rho0_sgx .* ifft (FFT), if rho0_sgx is scalar, uniform grid.
   void Compute_dt_rho_sg_mul_ifft_div_2(const float dt_rho_0_sgx,
                                         TFFTWComplexMatrix& FFT);

   /// Compute dt ./ rho0_sgx .* ifft (FFT), if rho0_sgx is scalar, non uniform grid, x component.
   void Compute_dt_rho_sg_mul_ifft_div_2_scalar_nonuniform_x(const float         dt_rho_0_sgx,
                                                             const TRealMatrix & dxudxn_sgx,
                                                             TFFTWComplexMatrix& FFT);
   /// Compute dt./rho0_sgx .* ifft (FFT), if rho0_sgx is scalar, non uniform grid, y component.
   void Compute_dt_rho_sg_mul_ifft_div_2_scalar_nonuniform_y(const float          dt_rho_0_sgy,
                                                             const TRealMatrix & dyudyn_sgy,
                                                             TFFTWComplexMatrix& FFT);
   /// Compute dt./rho0_sgx .* ifft (FFT), if rho0_sgx is scalar, non uniform grid, z component.
   void Compute_dt_rho_sg_mul_ifft_div_2_scalar_nonuniform_z(const float         dt_rho_0_sgz,
                                                             const TRealMatrix & dzudzn_sgz,
                                                             TFFTWComplexMatrix& FFT);

   //------------------------- X dimension -----------------------------------//
   /// Compute a new value of ux_sgx, default case
   void Compute_ux_sgx_normalize(const TRealMatrix& FFT_p,
                                 const TRealMatrix& dt_rho0,
                                 const TRealMatrix& pml);
   /// Compute a new value of ux_sgx, scalar, uniform case
   void Compute_ux_sgx_normalize_scalar_uniform(const TRealMatrix& FFT_p,
                                                const float        dt_rho0,
                                                const TRealMatrix& pml);
   /// Compute a new value of ux_sgx, scalar, non-uniform case
   void Compute_ux_sgx_normalize_scalar_nonuniform(const TRealMatrix& FFT_p,
                                                   const float        dt_rho0,
                                                   const TRealMatrix& dxudxn_sgx,
                                                   const TRealMatrix& pml);

  //------------------------- Y dimension ------------------------------------//
   /// Compute a new value of uy_sgy, default case
   void Compute_uy_sgy_normalize(const TRealMatrix& FFT_p,
                                 const TRealMatrix& dt_rho0,
                                 const TRealMatrix& pml);
   /// Compute a new value of uy_sgy, scalar, uniform case
   void Compute_uy_sgy_normalize_scalar_uniform(const TRealMatrix& FFT_p,
                                                const float        dt_rho0,
                                                const TRealMatrix& pml);
   /// Compute a new value of uy_sgy, scalar, non-uniform case
   void Compute_uy_sgy_normalize_scalar_nonuniform(const TRealMatrix& FFT_p,
                                                   const float        dt_rho0,
                                                   const TRealMatrix& dyudyn_sgy,
                                                   const TRealMatrix& pml);

   //------------------------- Z dimension -----------------------------------//
   /// Compute a new value for uz_sgz, default case.
   void Compute_uz_sgz_normalize(const TRealMatrix& FFT_p,
                                 const TRealMatrix& dt_rho0,
                                 const TRealMatrix& pml);
   /// Compute a new value for uz_sgz, scalar, uniform case.
   void Compute_uz_sgz_normalize_scalar_uniform(const TRealMatrix& FFT_p,
                                                const float        dt_rho0,
                                                const TRealMatrix& pml);
   /// Compute a new value for uz_sgz, scalar, non-uniform case.
   void Compute_uz_sgz_normalize_scalar_nonuniform(const TRealMatrix& FFT_p,
                                                   const float        dt_rho0,
                                                   const TRealMatrix& dzudzn_sgz,
                                                   const TRealMatrix& pml);


   /// Add transducer data  source to X component.
   void AddTransducerSource(const TIndexMatrix& u_source_index,
                                  TIndexMatrix& delay_mask,
                            const TRealMatrix & transducer_signal);

   /// Add in velocity source terms.
   void Add_u_source(const TRealMatrix & u_source_input,
                     const TIndexMatrix& u_source_index,
                     const size_t        t_index,
                     const size_t        u_source_mode,
                     const size_t        u_source_many);

   /// Destructor
   virtual ~Tuxyz_sgxyzMatrix() {};

protected:
    // Default constructor not allowed for public.
    Tuxyz_sgxyzMatrix() : TRealMatrix() {};

    /**
     * Copy constructor not allowed for public.
     * @param src
     */
    Tuxyz_sgxyzMatrix(const Tuxyz_sgxyzMatrix& src);

    /// operator = not allowed for public.
    Tuxyz_sgxyzMatrix& operator = (const Tuxyz_sgxyzMatrix& src);
private:

}; // end of Tuxyz_sgxyzMatrix

#endif	/* UX_SGREALMATRIX_H */

