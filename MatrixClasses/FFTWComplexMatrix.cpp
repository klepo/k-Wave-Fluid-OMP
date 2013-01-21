/**
 * @file        FFTWComplexMatrix.cpp
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au
 * 
 * @brief       The implementation file containing the class that implements 
 *              3D FFT using the FFTW interface
 * 
 * @version     kspaceFirstOrder3D 2.13
 * @date        9 August 2011, 13:10    (created) \n
 *              17 September 2012, 15:35 (revised)
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

#include <MatrixClasses/FFTWComplexMatrix.h>
#include <MatrixClasses/RealMatrix.h>

#include <Utils/ErrorMessages.h>

using namespace std;



//----------------------------------------------------------------------------//
//                              Constants                                     //
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//                              Public methods                                //
//----------------------------------------------------------------------------//


/**
 * Constructor
 * @param DimensionSizes - Dimension sizes of the reduced complex matrix
 */
TFFTWComplexMatrix::TFFTWComplexMatrix(struct TDimensionSizes DimensionSizes) :
        TComplexMatrix(),
        fftw_plan_R2C_Created(false), fftw_plan_C2R_Created(false)
{
                  
    InitDimensions(DimensionSizes);
    
    AllocateMemory();        
}// end of TFFTWComplexMatrix
//------------------------------------------------------------------------------



/**
 * Destructor
 */
TFFTWComplexMatrix::~TFFTWComplexMatrix(){
        
    //-- Destroy fftw plans --//
    
    FreeMemory();
    
    if (fftw_plan_R2C_Created)   fftwf_destroy_plan(fftw_plan_R2C);    
    if (fftw_plan_C2R_Created)   fftwf_destroy_plan(fftw_plan_C2R);
    
    
}// end of ~TFFTWComplexMatrix()
//------------------------------------------------------------------------------



/**
 * Create FFTW plan for Real-to-Complex.
 * @param [in,out] InMatrix  - RealMatrix of which to create the plan
 * @warning unless FFTW_ESTIMATE flag is specified, the content of the InMatrix will be destroyed!
 */    
void TFFTWComplexMatrix::CreateFFTPlan3D_R2C(TRealMatrix& InMatrix){
            
    
    fftw_plan_R2C = fftwf_plan_dft_r2c_3d(InMatrix.GetDimensionSizes().Z,InMatrix.GetDimensionSizes().Y,InMatrix.GetDimensionSizes().X,
                                                    InMatrix.GetRawData(), (fftwf_complex *) pMatrixData,
                                                    TFFTWComplexMatrix_FFT_FLAG);
    
    fftw_plan_R2C_Created  = true;
    
    
}// end of CreateFFTPlan3D_RealToComplex
//------------------------------------------------------------------------------

/**
 * Create FFTW plan for Complex-to-Real.
 * @param [in, out] OutMatrix - RealMatrix of which to create the plan.  
 * @warning unless FFTW_ESTIMATE flag is specified, the content of the InMatrix will be destroyed!
 */    
void TFFTWComplexMatrix::CreateFFTPlan3D_C2R(TRealMatrix& OutMatrix){
        
    
    fftw_plan_C2R = fftwf_plan_dft_c2r_3d(OutMatrix.GetDimensionSizes().Z,OutMatrix.GetDimensionSizes().Y,OutMatrix.GetDimensionSizes().X,
                                                    (fftwf_complex *)pMatrixData, OutMatrix.GetRawData(), 
                                                     TFFTWComplexMatrix_FFT_FLAG);
    
    fftw_plan_C2R_Created = true;    
    
}//end of CreateFFTPlan3D_ComplexToReal
//------------------------------------------------------------------------------



    
/**
 * Computer forward 3D Real-to-Complex FFT.
 * @param [in] InMatrix - Input Matrix
 */
void TFFTWComplexMatrix::Compute_FFT_3D_R2C(TRealMatrix & InMatrix){
    
    
    // execute real to complex fft
    fftwf_execute_dft_r2c(fftw_plan_R2C, InMatrix.GetRawData(), (fftwf_complex *) pMatrixData);    
    
    
    
}// end of Compute_FFT_3D_r2c
//------------------------------------------------------------------------------

/**
 * Compute inverse 3D Complex to Real FFT.
 * @param [out] OutMatrix
 */
void TFFTWComplexMatrix::Compute_iFFT_3D_C2R(TRealMatrix & OutMatrix){
    
    fftwf_execute_dft_c2r(fftw_plan_C2R,(fftwf_complex *) pMatrixData, OutMatrix.GetRawData());
        
}// end of Compute_iFFT_3D_c2r
//------------------------------------------------------------------------------



//----------------------------------------------------------------------------//
//                           Protected methods                                //
//----------------------------------------------------------------------------//

/**
 * Allocate Memory using fftwf_malloc function to ensure correct alignment
 * 
 */
 void TFFTWComplexMatrix::AllocateMemory(){
     
    /* No memory allocated before this function*/
    
    
    pMatrixData = (float *) fftwf_malloc(pTotalAllocatedElementCount * sizeof (float));
    
    
    if (!pMatrixData) {
        fprintf(stderr,Matrix_ERR_FMT_NotEnoughMemory, "TFFTWComplexMatrix");
        throw bad_alloc();
        
    }
    
    //-- first touch --//
     #ifndef __NO_OMP__       
        #pragma omp parallel for schedule(static)
     #endif  
    for (size_t i=0; i<pTotalAllocatedElementCount; i++){
        pMatrixData[i] = 0.0f;
    }
    
    
 }// end of virtual bool FreeMemory
 //-----------------------------------------------------------------------------
    
 
 /**
  * Free memory using fftwf_free
  */
 void TFFTWComplexMatrix::FreeMemory(){
      
     if (pMatrixData) fftwf_free( pMatrixData);
     pMatrixData = NULL;
     
 }// end of FreeMemory
 //-----------------------------------------------------------------------------
    
 
 
 
 

//----------------------------------------------------------------------------//
//                           Private methods                                  //
//----------------------------------------------------------------------------//