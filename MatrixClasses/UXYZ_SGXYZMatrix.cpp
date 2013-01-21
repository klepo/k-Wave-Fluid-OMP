/**
 * @file        UXYZ_SGXYZMatrix.cpp
 * @author      Jiri Jaros
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au   \n
 * 
 * @brief       The implementation file containing the particle velocity matrix
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



#include <MatrixClasses/UXYZ_SGXYZMatrix.h>
#include <MatrixClasses/FFTWComplexMatrix.h>


using namespace std;
//----------------------------------------------------------------------------//
//                              Constants                                     //
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
//                              Implementation                                //
//                              public methods                                //
//----------------------------------------------------------------------------//



/**
 * Compute dt./rho0_sgx .* the content of the matrix 
 * @param [in] dt_rho0_sg - matrix with the component of dt .* rho0_sg{x,y,z}
 * @param [in] FFT        - FFT matrix
 */
void Tuxyz_sgxyzMatrix::Compute_dt_rho_sg_mul_ifft_div_2(TRealMatrix& dt_rho0_sg, TFFTWComplexMatrix& FFT){
    
    FFT.Compute_iFFT_3D_C2R(*this);
        
    
    const float Divider = 1.0f/(2.0f *pTotalElementCount);
    //dt_rho0_sgx .* real...
    
#ifndef __NO_OMP__            
    #pragma omp parallel for schedule (static)     
#endif
    for (size_t i = 0; i < pTotalElementCount; i++){
        pMatrixData[i] = dt_rho0_sg[i] * (pMatrixData[i] * Divider);       
    }
    

}// end of Compute_dt_div_rho_sgx_mul
//------------------------------------------------------------------------------


   
/**
 * Compute dt./rho0_sgx .* ifft (FFT), if rho0_sgx is scalar, uniform case
 * @param [in] dt_rho_0_sgx     - scalar value
 * @param [in] FFT              - FFT matrix
 */
void Tuxyz_sgxyzMatrix::Compute_dt_rho_sg_mul_ifft_div_2(float dt_rho_0_sgx, TFFTWComplexMatrix& FFT){
    FFT.Compute_iFFT_3D_C2R(*this);
        
    
    const float Divider = 1.0f/(2.0f *pTotalElementCount) * dt_rho_0_sgx;
    //dt_rho0_sgx .* real...
    
#ifndef __NO_OMP__            
    #pragma omp parallel for schedule (static)     
#endif
    for (size_t i = 0; i < pTotalElementCount; i++){
        pMatrixData[i] = pMatrixData[i] * Divider;       
    }
      
    
}// end of Compute_dt_rho_sg_mul_ifft_div_2
//------------------------------------------------------------------------------
   

/**
 * Compute dt./rho0_sgx .* ifft (FFT), when rho0_sgx is scalar, nonuniform
 * @param [in] dt_rho_0_sgx     - scalar value
 * @param [in] dxudxn_sgx       - non-uniform mapping
 * @param [in] FFT              - FFT matrix
 */
void Tuxyz_sgxyzMatrix::Compute_dt_rho_sg_mul_ifft_div_2_scalar_nonuniform_x
        (float dt_rho_0_sgx, TRealMatrix & dxudxn_sgx, TFFTWComplexMatrix& FFT){
    
    
    
    FFT.Compute_iFFT_3D_C2R(*this);        
    
    const float Divider = 1.0f/(2.0f *pTotalElementCount) * dt_rho_0_sgx;
    //dt_rho0_sgx .* real...
    
        
    #ifndef __NO_OMP__            
        #pragma omp for schedule (static) 
    #endif    
    for (size_t z = 0; z < pDimensionSizes.Z; z++){

        register size_t i = z* pDimensionSizes.Y * pDimensionSizes.X;                        
        for (size_t y = 0; y< pDimensionSizes.Y; y++){                                               
            for (size_t x = 0; x < pDimensionSizes.X; x++){              
                pMatrixData[i] = pMatrixData[i] * Divider * dxudxn_sgx[x];
                i++;                    
            } // x
        } // y
      } // z

    
    
    
}// end of Compute_dt_rho_sg_mul_ifft_div_2_scalar_nonuniform_x
//------------------------------------------------------------------------------


/**
 * Compute dt./rho0_sgy .* ifft (FFT), if rho0_sgx is scalar, nonuniform
 * @param [in] dt_rho_0_sgy     - scalar value
 * @param [in] dyudyn_sgy       - non-uniform mapping
 * @param [in] FFT              - FFT matrix
 */
void Tuxyz_sgxyzMatrix::Compute_dt_rho_sg_mul_ifft_div_2_scalar_nonuniform_y
        (float dt_rho_0_sgy, TRealMatrix & dyudyn_sgy, TFFTWComplexMatrix& FFT){
    
    
    
    FFT.Compute_iFFT_3D_C2R(*this);        
    
    const float Divider = 1.0f/(2.0f *pTotalElementCount) * dt_rho_0_sgy;
    //dt_rho0_sgx .* real...
    
        
    #ifndef __NO_OMP__            
        #pragma omp for schedule (static) 
    #endif    
    for (size_t z = 0; z < pDimensionSizes.Z; z++){

        register size_t i = z* pDimensionSizes.Y * pDimensionSizes.X;                        
        for (size_t y = 0; y< pDimensionSizes.Y; y++){                                               
            const float dyudyn_sgy_data = dyudyn_sgy[y] * Divider;
            for (size_t x = 0; x < pDimensionSizes.X; x++){              
                pMatrixData[i] = pMatrixData[i] * dyudyn_sgy_data;
                i++;                    
            } // x
        } // y
      } // z

    
    
    
}// end of Compute_dt_rho_sg_mul_ifft_div_2_scalar_nonuniform_y
//------------------------------------------------------------------------------

/**
 * Compute dt./rho0_sgz .* ifft (FFT), if rho0_sgx is scalar, uniform
 * @param [in] dt_rho_0_sgz - scalar value
 * @param [in] dzudzn_sgz
 * @param [in] FFT          - FFT matrix
 */
void Tuxyz_sgxyzMatrix::Compute_dt_rho_sg_mul_ifft_div_2_scalar_nonuniform_z
    (float dt_rho_0_sgz, TRealMatrix & dzudzn_sgz, TFFTWComplexMatrix& FFT){
        
    
    
    
    FFT.Compute_iFFT_3D_C2R(*this);        
    
    const float Divider = 1.0f/(2.0f *pTotalElementCount) * dt_rho_0_sgz;
    //dt_rho0_sgx .* real...
    
        
    #ifndef __NO_OMP__            
        #pragma omp for schedule (static) 
    #endif    
    for (size_t z = 0; z < pDimensionSizes.Z; z++){

        register size_t i = z* pDimensionSizes.Y * pDimensionSizes.X;                        
        const float dzudzn_sgz_data = dzudzn_sgz[z] * Divider;
        
        for (size_t y = 0; y< pDimensionSizes.Y; y++){                                                           
            for (size_t x = 0; x < pDimensionSizes.X; x++){              
                pMatrixData[i] = pMatrixData[i] * dzudzn_sgz_data;
                i++;                    
            } // x
        } // y
      } // z

    
        
 }// end of Compute_dt_rho_sg_mul_ifft_div_2_scalar_nonuniform_z
//------------------------------------------------------------------------------
   


/**
 * Compute new value for ux_sgx.
 *
 * @param [in] FFT_p   - fft of pressure 
 * @param [in] dt_rho0 - dt_rho0_sgx
 * @param [in] pml     - pml_x
 */
void Tuxyz_sgxyzMatrix::Compute_ux_sgx_normalize(TRealMatrix& FFT_p, TRealMatrix& dt_rho0, TRealMatrix& pml){

    const float Divider = 1.0f / pTotalElementCount;              
    
    #ifndef __NO_OMP__            
        #pragma omp for schedule (static) 
    #endif
    for (size_t z = 0; z < pDimensionSizes.Z; z++){        
        
        register size_t i = z* p2DDataSliceSize;        
        for (size_t y = 0; y< pDimensionSizes.Y; y++){                    
            for (size_t x = 0; x < pDimensionSizes.X; x++){
                                                
                register float pMatrixElement = pMatrixData[i];
                
                //FFT_p.ElementMultiplyMatrices(dt_rho0);
                const float FFT_p_el = Divider * FFT_p[i] * dt_rho0[i];
                
                //BSXElementRealMultiply_1D_X(abc);                
                pMatrixElement *= pml[x];
                 
                 //ElementSubMatrices(FFT_p);
                pMatrixElement -= FFT_p_el;
                
                //BSXElementRealMultiply_1D_X(abc);
                pMatrixData[i] = pMatrixElement * pml[x];
                        
                i++;
            } // x
        } // y
    } // z
    
}// end of Compute_ux_sgx_normalize_Optimized
//------------------------------------------------------------------------------

/**
 * Compute_ux_sgx_normalize if rho0 is a scalar.
 * @param  [in] FFT_p   - matrix
 * @param  [in] dt_rho0 - scalar
 * @param  [in] pml     - matrix
 */
void Tuxyz_sgxyzMatrix::Compute_ux_sgx_normalize_scalar_uniform(TRealMatrix& FFT_p, float dt_rho0, TRealMatrix& pml){
    
    const float Divider = dt_rho0 / pTotalElementCount;              
    
    #ifndef __NO_OMP__            
        #pragma omp for schedule (static) 
    #endif
    for (size_t z = 0; z < pDimensionSizes.Z; z++){        
        
        register size_t i = z* p2DDataSliceSize;        
        for (size_t y = 0; y< pDimensionSizes.Y; y++){                    
            for (size_t x = 0; x < pDimensionSizes.X; x++){
                                                
                register float pMatrixElement = pMatrixData[i];
                
                //FFT_p.ElementMultiplyMatrices(dt_rho0);
                const float FFT_p_el = Divider * FFT_p[i];
                
                //BSXElementRealMultiply_1D_X(abc);                
                pMatrixElement *= pml[x];
                 
                 //ElementSubMatrices(FFT_p);
                pMatrixElement -= FFT_p_el;
                
                //BSXElementRealMultiply_1D_X(abc);
                pMatrixData[i] = pMatrixElement * pml[x];
                        
                i++;
            } // x
        } // y
    } // z
    
    
}// end of Compute_ux_sgx_normalize
//------------------------------------------------------------------------------
  

/**
 * Compute_ux_sgx_normalize if rho0 is a scalar, non uniform.
 * @param [in] FFT_p      - matrix
 * @param [in] dt_rho0    - scalar
 * @param [in] dxudxn_sgx - scalar
 * @param [in] pml        - matrix
 */
void Tuxyz_sgxyzMatrix::Compute_ux_sgx_normalize_scalar_nonuniform(TRealMatrix& FFT_p, float dt_rho0, TRealMatrix & dxudxn_sgx, TRealMatrix& pml){
    
    const float Divider = dt_rho0 / pTotalElementCount;              
    
    #ifndef __NO_OMP__            
        #pragma omp for schedule (static) 
    #endif
    for (size_t z = 0; z < pDimensionSizes.Z; z++){        
        
        register size_t i = z* p2DDataSliceSize;        
        for (size_t y = 0; y< pDimensionSizes.Y; y++){                    
            for (size_t x = 0; x < pDimensionSizes.X; x++){
                                                
                register float pMatrixElement = pMatrixData[i];
                
                //FFT_p.ElementMultiplyMatrices(dt_rho0);
                const float FFT_p_el = (Divider * dxudxn_sgx[x]) * FFT_p[i];
                
                //BSXElementRealMultiply_1D_X(abc);                
                pMatrixElement *= pml[x];
                 
                 //ElementSubMatrices(FFT_p);
                pMatrixElement -= FFT_p_el;
                
                //BSXElementRealMultiply_1D_X(abc);
                pMatrixData[i] = pMatrixElement * pml[x];
                        
                i++;
            } // x
        } // y
    } // z
    
    
}// end of Compute_ux_sgx_normalize_scalar_nonuniform
//------------------------------------------------------------------------------


/**
 *  Compute new value for uy_sgy 
 * 
 * @param [in] FFT_p   - fft of pressure 
 * @param [in] dt_rho0 - dt_rh0_sgy
 * @param [in] pml     - pml_y
 */
void Tuxyz_sgxyzMatrix::Compute_uy_sgy_normalize(TRealMatrix& FFT_p, TRealMatrix& dt_rho0, TRealMatrix& pml){
  
 
    const float Divider = 1.0f / pTotalElementCount;    
    #ifndef __NO_OMP__            
        #pragma omp for schedule (static) 
    #endif
    for (size_t z = 0; z < pDimensionSizes.Z; z++){
                
        size_t i = z* p2DDataSliceSize;
        for (size_t y = 0; y< pDimensionSizes.Y; y++){                    
            const float pml_y = pml[y];
            
            for (size_t x = 0; x < pDimensionSizes.X; x++){
                
                register float pMatrixElement = pMatrixData[i];
                
                //FFT_p.ElementMultiplyMatrices(dt_rho0);
                const float FFT_p_el = Divider * FFT_p[i] * dt_rho0[i] ;
                
                //BSXElementRealMultiply_1D_X(abc);                
                pMatrixElement *= pml_y;
                 
                 //ElementSubMatrices(FFT_p);
                pMatrixElement -= FFT_p_el;
                
                //BSXElementRealMultiply_1D_X(abc);
                pMatrixData[i] = pMatrixElement * pml_y;
                
                i++;
                        
            } // x
        } // y
    } // z
         
}// end of Compute_uy_sgy_normalize_Optimized
//------------------------------------------------------------------------------


/**
 * Compute_uy_sgy_normalize if rho0 is a scalar
 * @param [in] FFT_p    - matrix
 * @param [in] dt_rho0  - scalar
 * @param [in] pml      - matrix
 */
void Tuxyz_sgxyzMatrix::Compute_uy_sgy_normalize_scalar_uniform(TRealMatrix& FFT_p, float dt_rho0, TRealMatrix& pml){
    
    const float Divider = dt_rho0 / pTotalElementCount;    
    #ifndef __NO_OMP__            
        #pragma omp for schedule (static) 
    #endif
    for (size_t z = 0; z < pDimensionSizes.Z; z++){
                
        size_t i = z* p2DDataSliceSize;
        for (size_t y = 0; y< pDimensionSizes.Y; y++){                    
            const float pml_y = pml[y];
            
            for (size_t x = 0; x < pDimensionSizes.X; x++){
                
                register float pMatrixElement = pMatrixData[i];
                
                //FFT_p.ElementMultiplyMatrices(dt_rho0);
                const float FFT_p_el = Divider * FFT_p[i] ;
                
                //BSXElementRealMultiply_1D_X(abc);                
                pMatrixElement *= pml_y;
                 
                 //ElementSubMatrices(FFT_p);
                pMatrixElement -= FFT_p_el;
                
                //BSXElementRealMultiply_1D_X(abc);
                pMatrixData[i] = pMatrixElement * pml_y;
                
                i++;
                        
            } // x
        } // y
    } // z
    
    
}// end of Compute_ux_sgx_normalize
//------------------------------------------------------------------------------


/**
 * Compute_uy_sgy_normalize if rho0 is a scalar, non uniform
 * @param [in] FFT_p      - matrix
 * @param [in] dt_rho0    - scalar
 * @param [in] dyudyn_sgy - scalar
 * @param [in] pml        - matrix
 */
void Tuxyz_sgxyzMatrix::Compute_uy_sgy_normalize_scalar_nonuniform(TRealMatrix& FFT_p, float dt_rho0, TRealMatrix & dyudyn_sgy, TRealMatrix& pml){
   
        const float Divider = dt_rho0 / pTotalElementCount;    
    #ifndef __NO_OMP__            
        #pragma omp for schedule (static) 
    #endif
    for (size_t z = 0; z < pDimensionSizes.Z; z++){
                
        size_t i = z* p2DDataSliceSize;
        for (size_t y = 0; y< pDimensionSizes.Y; y++){                    
            const float pml_y = pml[y];
            const float dyudyn_sgy_data = dyudyn_sgy[y];
            for (size_t x = 0; x < pDimensionSizes.X; x++){
                
                register float pMatrixElement = pMatrixData[i];
                
                //FFT_p.ElementMultiplyMatrices(dt_rho0);
                const float FFT_p_el = (Divider * dyudyn_sgy_data) * FFT_p[i] ;
                
                //BSXElementRealMultiply_1D_X(abc);                
                pMatrixElement *= pml_y;
                 
                 //ElementSubMatrices(FFT_p);
                pMatrixElement -= FFT_p_el;
                
                //BSXElementRealMultiply_1D_X(abc);
                pMatrixData[i] = pMatrixElement * pml_y;
                
                i++;
                        
            } // x
        } // y
    } // z
    
    
}// end of Compute_uy_sgy_normalize_scalar_nonuniform
//------------------------------------------------------------------------------


/**
 * Compute new value of uz_sgz
 * 
 * @param [in] FFT_p   - fft of pressure 
 * @param [in] dt_rho0 - dt_rh0_sgz
 * @param [in] pml     - pml_z
 */
void Tuxyz_sgxyzMatrix::Compute_uz_sgz_normalize(TRealMatrix& FFT_p, TRealMatrix& dt_rho0, TRealMatrix& pml){

      const float Divider = 1.0f / pTotalElementCount;
      
    #ifndef __NO_OMP__            
        #pragma omp for schedule (static) 
    #endif
      for (size_t z = 0; z < pDimensionSizes.Z; z++){        
            size_t i = z* p2DDataSliceSize;
            const float pml_z = pml[z];

            for (size_t y = 0; y< pDimensionSizes.Y; y++){                    
                for (size_t x = 0; x < pDimensionSizes.X; x++){                

                    register float pMatrixElement = pMatrixData[i];

                    //FFT_p.ElementMultiplyMatrices(dt_rho0);
                    const float FFT_p_el = Divider * FFT_p[i] * dt_rho0[i];

                    //BSXElementRealMultiply_1D_X(abc);                
                    pMatrixElement *= pml_z;

                     //ElementSubMatrices(FFT_p);
                    pMatrixElement -= FFT_p_el;

                    //BSXElementRealMultiply_1D_X(abc);

                    pMatrixData[i] = pMatrixElement * pml_z;

                    i++;

                } // x
            } // y
        } // z

}// end of Compute_uz_sgz_normalize_fast
//------------------------------------------------------------------------------

/**
 * Compute_uz_sgz_normalize if rho0 is a scalar
 * @param [in] FFT_p   - matrix
 * @param [in] dt_rho0 - scalar
 * @param [in] pml     - matrix
 */
void Tuxyz_sgxyzMatrix::Compute_uz_sgz_normalize_scalar_uniform(TRealMatrix& FFT_p, float& dt_rho0, TRealMatrix& pml){
    
      const float Divider = dt_rho0 / pTotalElementCount;
      
    #ifndef __NO_OMP__            
        #pragma omp for schedule (static) 
    #endif
      for (size_t z = 0; z < pDimensionSizes.Z; z++){        
            size_t i = z* p2DDataSliceSize;
            const float pml_z = pml[z];

            for (size_t y = 0; y< pDimensionSizes.Y; y++){                    
                for (size_t x = 0; x < pDimensionSizes.X; x++){                

                    register float pMatrixElement = pMatrixData[i];

                    //FFT_p.ElementMultiplyMatrices(dt_rho0);
                    const float FFT_p_el = Divider * FFT_p[i];

                    //BSXElementRealMultiply_1D_X(abc);                
                    pMatrixElement *= pml_z;

                     //ElementSubMatrices(FFT_p);
                    pMatrixElement -= FFT_p_el;

                    //BSXElementRealMultiply_1D_X(abc);

                    pMatrixData[i] = pMatrixElement * pml_z;

                    i++;

                } // x
            } // y
        } // z

    
    
}// end of Compute_uz_sgz_normalize_scalar_uniform
//------------------------------------------------------------------------------


/**
 * Compute_uz_sgz_normalize if rho0 is a scalar, non uniform
 * @param [in] FFT_p      - matrix
 * @param [in] dt_rho0    - scalar
 * @param [in] dzudzn_sgz - scalar
 * @param [in] pml        - matrix
 */
void Tuxyz_sgxyzMatrix::Compute_uz_sgz_normalize_scalar_nonuniform(TRealMatrix& FFT_p, float& dt_rho0,TRealMatrix & dzudzn_sgz, TRealMatrix& pml){
     
    const float Divider = dt_rho0 / pTotalElementCount;
      
    #ifndef __NO_OMP__            
        #pragma omp for schedule (static) 
    #endif
      for (size_t z = 0; z < pDimensionSizes.Z; z++){        
            size_t i = z* p2DDataSliceSize;
            const float pml_z = pml[z];
            const float dzudzn_sgz_data = dzudzn_sgz[z];
            
            for (size_t y = 0; y< pDimensionSizes.Y; y++){                    
                for (size_t x = 0; x < pDimensionSizes.X; x++){                

                    register float pMatrixElement = pMatrixData[i];

                    //FFT_p.ElementMultiplyMatrices(dt_rho0);
                    const float FFT_p_el = (Divider * dzudzn_sgz_data)* FFT_p[i];

                    //BSXElementRealMultiply_1D_X(abc);                
                    pMatrixElement *= pml_z;

                     //ElementSubMatrices(FFT_p);
                    pMatrixElement -= FFT_p_el;

                    //BSXElementRealMultiply_1D_X(abc);

                    pMatrixData[i] = pMatrixElement * pml_z;

                    i++;

                } // x
            } // y
        } // z

    
}// end of Compute_uz_sgz_normalize_scalar_nonuniform
//------------------------------------------------------------------------------


/**
 * Add transducer data to X dimension
 * @param [in] u_source_index      - long index matrix
 * @param [in, out] delay_mask     - long index matrix   - modified inside (+1)
 * @param [in] transducer_signal   - transducer signal
 */
void Tuxyz_sgxyzMatrix::AddTransducerSource(TLongMatrix& u_source_index, TLongMatrix& delay_mask, TRealMatrix& transducer_signal){

//#ifndef __NO_OMP__            
//    #pragma omp parallel for schedule (static)     
//#endif        
    for (size_t i = 0; i < u_source_index.GetTotalElementCount(); i++){
        pMatrixData[u_source_index[i]] += transducer_signal[delay_mask[i]];            
        delay_mask[i]++;    
    }
    
    
}// end of AddTransducerSource
//------------------------------------------------------------------------------
   
   


/**
 * Add in velocity source terms
 * 
 * @param [in] u_source_input       - Source input to add
 * @param [in] u_source_index       - long index matrix
 * @param [in] t_index              - actual time step
 * @param [in] u_source_mode        - Mode 0 = dirichlet boundary, 1 = add in
 * @param [in] u_source_many        - 0 = One series, 1 = multiple series
 * 
 */
void Tuxyz_sgxyzMatrix::Add_u_source(TRealMatrix &u_source_input, TLongMatrix & u_source_index, int t_index, long u_source_mode, long u_source_many){
    
    
    size_t index2D = t_index;        
    if (u_source_many != 0) { // is 2D
        index2D = t_index * u_source_index.GetTotalElementCount();
    }
    
    
    
    if (u_source_mode == 0){            
        for (size_t i = 0; i < u_source_index.GetTotalElementCount(); i++){
            
            pMatrixData[u_source_index[i]] = u_source_input[index2D];
            
            if (u_source_many != 0) index2D ++;
        }    
    }// end of dirichlet
    
    if (u_source_mode == 1){            
        for (size_t i = 0; i < u_source_index.GetTotalElementCount(); i++){
            
            pMatrixData[u_source_index[i]] += u_source_input[index2D];           
            
            if (u_source_many != 0) index2D ++;    
        }            
    }// end of add
    
}// end of Add_u_source
//------------------------------------------------------------------------------   







   
//----------------------------------------------------------------------------//
//                              Implementation                                //
//                              private methods                                //
//----------------------------------------------------------------------------//
