/**
 * @file        UXYZ_SGXYZMatrix.h
 * @author      Jiri Jaros
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au   \n
 * 
 * @brief       The header file containing the particle velocity matrix
 * 
 * @version     kspaceFirstOrder3D 2.13
 * @date        28 July 2011, 11:37             (created) \n
 *              17 September 2012, 15:25        (revised)
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


#ifndef UXYZ_SGXYZREALMATRIX_H
#define	UXYZ_SGXYZREALMATRIX_H


#include <MatrixClasses/RealMatrix.h>
#include <MatrixClasses/FFTWComplexMatrix.h>

/**
 * @class Tuxyz_sgxyzMatrix
 * @brief The velocity matrix 
 */
class Tuxyz_sgxyzMatrix : public TRealMatrix {
public:
    
    /**
     * @brief Constructor
     * @param [in] DimensionSizes
     */
    Tuxyz_sgxyzMatrix(struct TDimensionSizes DimensionSizes) :
                       TRealMatrix(DimensionSizes) {};    
    
    
 
   /// compute this formula dt./rho0_sgx .* ifft (FFT)
   void Compute_dt_rho_sg_mul_ifft_div_2(TRealMatrix& dt_rho_0_sgx, TFFTWComplexMatrix& FFT);      
   /// compute this formula dt./rho0_sgx .* ifft (FFT),  if rho0_sgx is scalar, uniform grid
   void Compute_dt_rho_sg_mul_ifft_div_2(float dt_rho_0_sgx, TFFTWComplexMatrix& FFT);
   /// compute this formula dt./rho0_sgx .* ifft (FFT),  if rho0_sgx is scalar, non uniform grid, x component
   void Compute_dt_rho_sg_mul_ifft_div_2_scalar_nonuniform_x(float dt_rho_0_sgx, TRealMatrix & dxudxn_sgx, TFFTWComplexMatrix& FFT);
   /// compute this formula dt./rho0_sgx .* ifft (FFT),  if rho0_sgx is scalar, non uniform grid, y component
   void Compute_dt_rho_sg_mul_ifft_div_2_scalar_nonuniform_y(float dt_rho_0_sgy, TRealMatrix & dyudyn_sgy, TFFTWComplexMatrix& FFT);
   /// compute this formula dt./rho0_sgx .* ifft (FFT),  if rho0_sgx is scalar, non uniform grid, z component
   void Compute_dt_rho_sg_mul_ifft_div_2_scalar_nonuniform_z(float dt_rho_0_sgz, TRealMatrix & dzudzn_sgz, TFFTWComplexMatrix& FFT);

   /// Compute new value of ux_sgx, default case
   void Compute_ux_sgx_normalize(TRealMatrix& FFT_p, TRealMatrix& dt_rho0, TRealMatrix& pml);
   /// Compute new value of ux_sgx, scalar, uniform case
   void Compute_ux_sgx_normalize_scalar_uniform(TRealMatrix& FFT_p, float dt_rho0, TRealMatrix& pml);
   /// Compute new value of ux_sgx, scalar, non-uniform case
   void Compute_ux_sgx_normalize_scalar_nonuniform(TRealMatrix& FFT_p, float dt_rho0, TRealMatrix & dxudxn_sgx, TRealMatrix& pml);
      
   /// Compute new value of uy_sgy, default case
   void Compute_uy_sgy_normalize(TRealMatrix& FFT_p, TRealMatrix& dt_rho0, TRealMatrix& pml);
   /// Compute new value of uy_sgy, scalar, uniform case
   void Compute_uy_sgy_normalize_scalar_uniform(TRealMatrix& FFT_p, float dt_rho0, TRealMatrix& pml);   
   /// Compute new value of uy_sgy, scalar, non-uniform case
   void Compute_uy_sgy_normalize_scalar_nonuniform(TRealMatrix& FFT_p, float dt_rho0,TRealMatrix & dyudyn_sgy, TRealMatrix& pml);   
   
   /// Compute new value for uz_sgz, default case
   void Compute_uz_sgz_normalize(TRealMatrix& FFT_p, TRealMatrix& dt_rho0, TRealMatrix& pml);
   /// Compute new value for uz_sgz, scalar, uniform case
   void Compute_uz_sgz_normalize_scalar_uniform(TRealMatrix& FFT_p, float& dt_rho0, TRealMatrix& pml);
   /// Compute new value for uz_sgz, scalar, non-uniform case
   void Compute_uz_sgz_normalize_scalar_nonuniform(TRealMatrix& FFT_p, float& dt_rho0,TRealMatrix & dzudzn_sgz, TRealMatrix& pml);
   
   
   /// Add transducer data  source to X component
   void AddTransducerSource(TLongMatrix& u_source_index, TLongMatrix& delay_mask, TRealMatrix& transducer_signal);
   
   /// Add in velocity source terms
   void Add_u_source(TRealMatrix &u_source_input, TLongMatrix & u_source_index, int t_index, long u_source_mode, long u_source_many);
   
   /// Destructor
   virtual ~Tuxyz_sgxyzMatrix() {};
    
protected:
    // Default constructor not allowed for public
    Tuxyz_sgxyzMatrix() : TRealMatrix() {};
    /**
     * Copy constructor not allowed for public
     * @param src
     */
    Tuxyz_sgxyzMatrix(const Tuxyz_sgxyzMatrix& src);
    
    /// operator = not allowed for public
    Tuxyz_sgxyzMatrix& operator = (const Tuxyz_sgxyzMatrix& src);
private:
    
}; // end of Tuxyz_sgxyzMatrix

#endif	/* UX_SGREALMATRIX_H */

