/**
 * @file        FFTWComplexMatrix.cpp
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The implementation file containing the class that implements
 *              3D FFT using the FFTW interface
 *
 * @version     kspaceFirstOrder3D 2.16
 * @date        09 August    2011, 13:10 (created) \n
 *              24 August    2017, 09:25 (revised)
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

#include <MatrixClasses/FFTWComplexMatrix.h>
#include <MatrixClasses/RealMatrix.h>

#include <Parameters/Parameters.h>
#include <Utils/ErrorMessages.h>

#include <stdexcept>
#include <string.h>
#include <cstdio>
using namespace std;



//----------------------------------------------------------------------------//
//                              Constants                                     //
//----------------------------------------------------------------------------//
const string TFFTWComplexMatrix::FFTW_Wisdom_FileName_Extension = "FFTW_Wisdom";

//----------------------------------------------------------------------------//
//                              Public methods                                //
//----------------------------------------------------------------------------//


/**
 * Constructor.
 * @param [in] DimensionSizes - Dimension sizes of the reduced complex matrix
 */
TFFTWComplexMatrix::TFFTWComplexMatrix(const DimensionSizes & DimensionSizes) :
        TComplexMatrix(),
        fftw_plan_3D_R2C(NULL),  fftw_plan_3D_C2R(NULL),
        fftw_plan_1DX_R2C(NULL), fftw_plan_1DY_R2C(NULL), fftw_plan_1DZ_R2C(NULL),
        fftw_plan_1DX_C2R(NULL), fftw_plan_1DY_C2R(NULL), fftw_plan_1DZ_C2R(NULL)
{
  InitDimensions(DimensionSizes);
  AllocateMemory();
}// end of TFFTWComplexMatrix
//------------------------------------------------------------------------------



/**
 * Destructor.
 */
TFFTWComplexMatrix::~TFFTWComplexMatrix()
{
  FreeMemory();

  // free 3D plans
  if (fftw_plan_3D_R2C)  fftwf_destroy_plan(fftw_plan_3D_R2C);
  if (fftw_plan_3D_C2R)  fftwf_destroy_plan(fftw_plan_3D_C2R);

  //free 1D plans.
  if (fftw_plan_1DX_R2C) fftwf_destroy_plan(fftw_plan_1DX_R2C);
  if (fftw_plan_1DY_R2C) fftwf_destroy_plan(fftw_plan_1DY_R2C);
  if (fftw_plan_1DZ_R2C) fftwf_destroy_plan(fftw_plan_1DZ_R2C);

  if (fftw_plan_1DX_C2R) fftwf_destroy_plan(fftw_plan_1DX_C2R);
  if (fftw_plan_1DY_C2R) fftwf_destroy_plan(fftw_plan_1DY_C2R);
  if (fftw_plan_1DZ_C2R) fftwf_destroy_plan(fftw_plan_1DZ_C2R);

  fftw_plan_3D_R2C = NULL;
  fftw_plan_3D_C2R = NULL;

  fftw_plan_1DX_R2C = NULL;
  fftw_plan_1DY_R2C = NULL;
  fftw_plan_1DZ_R2C = NULL;

  fftw_plan_1DX_C2R = NULL;
  fftw_plan_1DY_C2R = NULL;
  fftw_plan_1DZ_C2R = NULL;

}// end of ~TFFTWComplexMatrix()
//------------------------------------------------------------------------------



/**
 * Create an FFTW plan for 3D Real-to-Complex.
 * @param [in,out] InMatrix  - RealMatrix of which to create the plan
 * @warning Unless FFTW_ESTIMATE flag is specified, the content of the InMatrix is destroyed!
 * @throw runtime_error if the plan can't be created.
 */
void TFFTWComplexMatrix::Create_FFT_Plan_3D_R2C(TRealMatrix& InMatrix)
{
  fftw_plan_3D_R2C = fftwf_plan_dft_r2c_3d(InMatrix.GetDimensionSizes().nz,
                                           InMatrix.GetDimensionSizes().ny,
                                           InMatrix.GetDimensionSizes().nx,
                                           InMatrix.GetRawData(),
                                           (fftwf_complex *) pMatrixData,
                                           TFFTWComplexMatrix_FFT_FLAG);

  if (!fftw_plan_3D_R2C)
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage, kErrFmtFftPlanNotCreated, "FFT_3D_R2C");
    throw runtime_error(ErrorMessage);
  }
}// end of Create_FFT_Plan_3D_R2C
//------------------------------------------------------------------------------

/**
 * Create an FFTW plan for 3D Complex-to-Real.
 * @param [in, out] OutMatrix - RealMatrix of which to create the plan.
 * @warning Unless FFTW_ESTIMATE flag is specified, the content of the InMatrix
 * is destroyed!
 * @throw runtime_error if the plan can't be created.
 */
void TFFTWComplexMatrix::Create_FFT_Plan_3D_C2R(TRealMatrix& OutMatrix)
{
  fftw_plan_3D_C2R = fftwf_plan_dft_c2r_3d(OutMatrix.GetDimensionSizes().nz,
                                           OutMatrix.GetDimensionSizes().ny,
                                           OutMatrix.GetDimensionSizes().nx,
                                           (fftwf_complex *) (pMatrixData),
                                           OutMatrix.GetRawData(),
                                           TFFTWComplexMatrix_FFT_FLAG);

  if (!fftw_plan_3D_C2R)
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage, kErrFmtFftPlanNotCreated, "FFT_3D_C2R");
    throw runtime_error(ErrorMessage);
  }
}//end of Create_FFT_Plan_3D_C2R
//------------------------------------------------------------------------------


/**
 * Create an FFTW plan for 1D Real-to-Complex in the X dimension.
 * There are two versions of this routine for GCC+FFTW and ICPC + MKL, otherwise it will not build!
 * The FFTW version processes the whole matrix at one while the MKL slab by slab
 * @param [in,out] InMatrix  - RealMatrix of which to create the plan
 * @warning Unless FFTW_ESTIMATE flag is specified, the content of the InMatrix is destroyed!
 * @throw runtime_error if the plan can't be created.
 */
void TFFTWComplexMatrix::Create_FFT_Plan_1DX_R2C(TRealMatrix& InMatrix)
{
  // the FFTW uses here 32b interface although it is internally 64b, it doesn't mind
  // since the size of 1 domain will never be bigger than 2^31 - however it it not a clear solution :)
  const int X   = static_cast<int> (InMatrix.GetDimensionSizes().nx);
  const int Y   = static_cast<int> (InMatrix.GetDimensionSizes().ny);
  const int Z   = static_cast<int> (InMatrix.GetDimensionSizes().nz);
  const int X_2 = ((X / 2) + 1);

  // 1D FFT definition - over the X axis
  const int  fft_rank = 1;
  fftw_iodim fft_dims[1];

  fft_dims[0].is = 1;
  fft_dims[0].n  = X;
  fft_dims[0].os = 1;

  // GNU Compiler + FFTW
  #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    // How many definition - other dims
    const int  fft_howmany_rank = 2;
    fftw_iodim fft_howmany_dims[2];

    // Z dim
    fft_howmany_dims[0].is = X * Y;
    fft_howmany_dims[0].n  = Z;
    fft_howmany_dims[0].os = X_2 * Y;

    // Y dim
    fft_howmany_dims[1].is = X;
    fft_howmany_dims[1].n  = Y;
    fft_howmany_dims[1].os = X_2;
  #endif

  // Intel Compiler + MKL
  #if (defined(__INTEL_COMPILER))
    const int  fft_howmany_rank = 1;
    fftw_iodim fft_howmany_dims[1];

    // Y dim
    fft_howmany_dims[0].is = X;
    fft_howmany_dims[0].n  = Y;
    fft_howmany_dims[0].os = X_2;
  #endif

  fftw_plan_1DX_R2C = fftwf_plan_guru_dft_r2c(fft_rank,                         // 1D FFT rank
                                              fft_dims,                         // 1D FFT dimensions of X

                                              fft_howmany_rank,                 // how many in Y and Z
                                              fft_howmany_dims,                 // Dims and strides in Y and Z

                                              InMatrix.GetRawData(),            // input data
                                              (fftwf_complex *) pMatrixData,    // output data
                                              TFFTWComplexMatrix_FFT_FLAG);     // flags

  if (!fftw_plan_1DX_R2C)
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage, kErrFmtFftPlanNotCreated, "FFT_1DX_R2C");
    throw runtime_error(ErrorMessage);
  }
}// end of Create_FFT_Plan_1DX_R2C
//------------------------------------------------------------------------------


/**
 * Create an FFTW plan for 1D Real-to-Complex in the Y dimension.
 * There are two versions of this routine for GCC+FFTW and ICPC + MKL, otherwise it will not build!
 * The FFTW version processes the whole matrix at one while the MKL slab by slab
 * @param [in,out] InMatrix  - RealMatrix of which to create the plan
 * @warning Unless FFTW_ESTIMATE flag is specified, the content of the InMatrix
 *          is destroyed!
 * @warning The FFTW matrix must be able to store 2 * (X * (Y/2 + 1) * Z) elements -
 *          possibly more than reduced dims!
 * @throw runtime_error if the plan can't be created.
 */
void TFFTWComplexMatrix::Create_FFT_Plan_1DY_R2C(TRealMatrix& InMatrix)
{
  // the FFTW uses here 32b interface although it is internally 64b, it doesn't mind
  // since the size of 1 domain will never be bigger than 2^31 - however it it not a clear solution :)
  const int X   = static_cast<int> (InMatrix.GetDimensionSizes().nx);
  const int Y   = static_cast<int> (InMatrix.GetDimensionSizes().ny);
  const int Z   = static_cast<int> (InMatrix.GetDimensionSizes().nz);
  const int Y_2 = ((Y / 2) + 1);

  // 1D FFT definition - over the Y axis
  const int  fft_rank = 1;
  fftw_iodim fft_dims[1];

  fft_dims[0].is = X;
  fft_dims[0].n  = Y;
  fft_dims[0].os = X;

  // GNU Compiler + FFTW
  #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    // How many definition - other dims
    const int  fft_howmany_rank = 2;
    fftw_iodim fft_howmany_dims[2];

    // Z dim
    fft_howmany_dims[0].is = X * Y;
    fft_howmany_dims[0].n  = Z;
    fft_howmany_dims[0].os = X * Y_2 ;

    // X dim
    fft_howmany_dims[1].is = 1;
    fft_howmany_dims[1].n  = X;
    fft_howmany_dims[1].os = 1;
  #endif

    // Intel Compiler + MKL
  #if (defined(__INTEL_COMPILER))
    const int  fft_howmany_rank = 1;
    fftw_iodim fft_howmany_dims[1];

    // X dim
    fft_howmany_dims[0].is = 1;
    fft_howmany_dims[0].n  = X;
    fft_howmany_dims[0].os = 1;
  #endif

  fftw_plan_1DY_R2C = fftwf_plan_guru_dft_r2c(fft_rank,                         // 1D FFT rank
                                              fft_dims,                         // 1D FFT dimensions of Y

                                              fft_howmany_rank,                 // how many in X and Z
                                              fft_howmany_dims,                 // Dims and strides in X and Z

                                              InMatrix.GetRawData(),            // input data
                                              (fftwf_complex *) pMatrixData,    // output data
                                              TFFTWComplexMatrix_FFT_FLAG);     // flags

  if (!fftw_plan_1DY_R2C)
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage, kErrFmtFftPlanNotCreated, "FFT_1DY_R2C");
    throw runtime_error(ErrorMessage);
  }
}// end of Create_FFT_Plan_1DY_R2C
//------------------------------------------------------------------------------


/**
 * Create an FFTW plan for 1D Real-to-Complex in the Z dimension.
 * There are two versions of this routine for GCC+FFTW and ICPC + MKL, otherwise it will not build!
 * The FFTW version processes the whole matrix at one while the MKL slab by slab
 * @param [in,out] InMatrix  - RealMatrix of which to create the plan
 * @warning Unless FFTW_ESTIMATE flag is specified, the content of the InMatrix
 *          is destroyed!
 * @warning The FFTW matrix must be able to store 2 * (X * Y * (Z/2+1) ) elements -
 *          possibly more than reduced dims!
 * @throw runtime_error if the plan can't be created.
 */
void TFFTWComplexMatrix::Create_FFT_Plan_1DZ_R2C(TRealMatrix& InMatrix)
{
  // the FFTW uses here 32b interface although it is internally 64b, it doesn't mind
  // since the size of 1 domain will never be bigger than 2^31 - however it it not a clear solution :)
  const int X   = static_cast<int> (InMatrix.GetDimensionSizes().nx);
  const int Y   = static_cast<int> (InMatrix.GetDimensionSizes().ny);
  const int Z   = static_cast<int> (InMatrix.GetDimensionSizes().nz);

  // 1D FFT definition - over the Y axis
  const int  fft_rank = 1;
  fftw_iodim fft_dims[1];

  fft_dims[0].is = X * Y;
  fft_dims[0].n  = Z;
  fft_dims[0].os = X * Y;

  // GNU Compiler + FFTW
  #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    // How many definition - other dims
    const int  fft_howmany_rank = 2;
    fftw_iodim fft_howmany_dims[2];

    // Y dim
    fft_howmany_dims[0].is = X;
    fft_howmany_dims[0].n  = Y;
    fft_howmany_dims[0].os = X;

    // X dim
    fft_howmany_dims[1].is = 1;
    fft_howmany_dims[1].n  = X;
    fft_howmany_dims[1].os = 1;
  #endif

  // Intel Compiler + MKL
  #if (defined(__INTEL_COMPILER))
    const int  fft_howmany_rank = 1;
    fftw_iodim fft_howmany_dims[1];

    // X dim
    fft_howmany_dims[0].is = 1;
    fft_howmany_dims[0].n  = X;
    fft_howmany_dims[0].os = 1;
  #endif

  fftw_plan_1DZ_R2C = fftwf_plan_guru_dft_r2c(fft_rank,                         // 1D FFT rank
                                              fft_dims,                         // 1D FFT dimensions of Y

                                              fft_howmany_rank,                 // how many in X and Z
                                              fft_howmany_dims,                 // Dims and strides in X and Z

                                              InMatrix.GetRawData(),            // input data
                                              (fftwf_complex *) pMatrixData,    // output data
                                              TFFTWComplexMatrix_FFT_FLAG);     // flags

  if (!fftw_plan_1DZ_R2C)
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage, kErrFmtFftPlanNotCreated, "FFT_1DZ_R2C");
    throw runtime_error(ErrorMessage);
  }
}// end of Create_FFT_Plan_1DZ_R2C
//------------------------------------------------------------------------------


/**
 * Create FFTW plan for Complex-to-Real in the X dimension.
 * There are two versions of this routine for GCC+FFTW and ICPC + MKL, otherwise it will not build!
 * The FFTW version processes the whole matrix at one while the MKL slab by slab
 * @param [in, out] OutMatrix - RealMatrix of which to create the plan.
 * @warning Unless FFTW_ESTIMATE flag is specified, the content of the InMatrix
 *          is destroyed!
 * @throw runtime_error if the plan can't be created.
 */
void TFFTWComplexMatrix::Create_FFT_Plan_1DX_C2R(TRealMatrix& OutMatrix)
{
  // the FFTW uses here 32b interface although it is internally 64b, it doesn't mind
  // since the size of 1 domain will never be bigger than 2^31 - however it it not a clear solution :)
  const int X   = static_cast<int> (OutMatrix.GetDimensionSizes().nx);
  const int Y   = static_cast<int> (OutMatrix.GetDimensionSizes().ny);
  const int Z   = static_cast<int> (OutMatrix.GetDimensionSizes().nz);
  const int X_2 = ((X / 2) + 1);

  // 1D FFT definition - over the X axis
  const int  fft_rank = 1;
  fftw_iodim fft_dims[1];

  fft_dims[0].is = 1;
  fft_dims[0].n  = X;
  fft_dims[0].os = 1;

  // GNU Compiler + FFTW
  #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    // How many definition - other dims
    const int  fft_howmany_rank = 2;
    fftw_iodim fft_howmany_dims[2];

    // Z dim
    fft_howmany_dims[0].is = X_2 * Y;
    fft_howmany_dims[0].n  = Z;
    fft_howmany_dims[0].os = X * Y;

    // Y dim
    fft_howmany_dims[1].is = X_2;
    fft_howmany_dims[1].n  = Y;
    fft_howmany_dims[1].os = X;
  #endif

  // Intel Compiler + MKL
  #if (defined(__INTEL_COMPILER))
    const int  fft_howmany_rank = 1;
    fftw_iodim fft_howmany_dims[1];

    // Y dim
    fft_howmany_dims[0].is = X_2;
    fft_howmany_dims[0].n  = Y;
    fft_howmany_dims[0].os = X;
  #endif

  fftw_plan_1DX_C2R = fftwf_plan_guru_dft_c2r(fft_rank,                         // 1D FFT rank
                                              fft_dims,                         // 1D FFT dimensions of X

                                              fft_howmany_rank,                 // how many in Y and Z
                                              fft_howmany_dims,                 // Dims and strides in Y and Z

                                              (fftwf_complex *) pMatrixData,    // input data
                                              OutMatrix.GetRawData(),           // output data
                                              TFFTWComplexMatrix_FFT_FLAG);     // flags

  if (!fftw_plan_1DX_C2R)
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage, kErrFmtFftPlanNotCreated, "FFT_1DX_C2R");
    throw runtime_error(ErrorMessage);
  }
}// end of Create_FFT_Plan_1DX_C2R
//------------------------------------------------------------------------------



/**
 * Create FFTW plan for Complex-to-Real in the Y dimension.
 * There are two versions of this routine for GCC+FFTW and ICPC + MKL, otherwise it will not build!
 * The FFTW version processes the whole matrix at one while the MKL slab by slab
 * @param [in, out] OutMatrix - RealMatrix of which to create the plan.
 * @warning Unless FFTW_ESTIMATE flag is specified, the content of the InMatrix is destroyed!
 * @throw runtime_error if the plan can't be created.
 */
void TFFTWComplexMatrix::Create_FFT_Plan_1DY_C2R(TRealMatrix& OutMatrix)
{
  // the FFTW uses here 32b interface although it is internally 64b, it doesn't mind
  // since the size of 1 domain will never be bigger than 2^31 - however it it not a clear solution :)
  const int X   = static_cast<int> (OutMatrix.GetDimensionSizes().nx);
  const int Y   = static_cast<int> (OutMatrix.GetDimensionSizes().ny);
  const int Z   = static_cast<int> (OutMatrix.GetDimensionSizes().nz);
  const int Y_2 = ((Y / 2) + 1);

  // 1D FFT definition - over the Y axis
  const int  fft_rank = 1;

  fftw_iodim fft_dims[1];
  fft_dims[0].is = X;
  fft_dims[0].n  = Y;
  fft_dims[0].os = X;

  // GNU Compiler + FFTW
  #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    // How many definition - other dims
    const int  fft_howmany_rank = 2;
    fftw_iodim fft_howmany_dims[2];

    // Z dim
    fft_howmany_dims[0].is = X * Y_2;
    fft_howmany_dims[0].n  = Z;
    fft_howmany_dims[0].os = X * Y;

    // X dim
    fft_howmany_dims[1].is = 1;
    fft_howmany_dims[1].n  = X;
    fft_howmany_dims[1].os = 1;
  #endif

// Intel Compiler + MKL
  #if (defined(__INTEL_COMPILER))
    const int  fft_howmany_rank = 1;
    fftw_iodim fft_howmany_dims[1];

    // X dim
    fft_howmany_dims[0].is = 1;
    fft_howmany_dims[0].n  = X;
    fft_howmany_dims[0].os = 1;
  #endif

  fftw_plan_1DY_C2R = fftwf_plan_guru_dft_c2r(fft_rank,                         // 1D FFT rank
                                              fft_dims,                         // 1D FFT dimensions of Y

                                              fft_howmany_rank,                 // how many in X and Z
                                              fft_howmany_dims,                 // Dims and strides in X and Z

                                              (fftwf_complex *) pMatrixData,    // input data
                                              OutMatrix.GetRawData(),           // output data
                                              TFFTWComplexMatrix_FFT_FLAG);     // flags

  if (!fftw_plan_1DY_C2R)
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage, kErrFmtFftPlanNotCreated, "FFT_1DY_C2R");
    throw runtime_error(ErrorMessage);
  }
}// end of Create_FFT_Plan_1DY_C2R
//------------------------------------------------------------------------------


/**
 * Create FFTW plan for Complex-to-Real in the Z dimension.
 * There are two versions of this routine for GCC+FFTW and ICPC + MKL, otherwise it will not build!
 * The FFTW version processes the whole matrix at one while the MKL slab by slab
 * @param [in, out] OutMatrix - RealMatrix of which to create the plan.
 * @warning Unless FFTW_ESTIMATE flag is specified, the content of the InMatrix is destroyed!
 * @throw runtime_error if the plan can't be created.
 */
void TFFTWComplexMatrix::Create_FFT_Plan_1DZ_C2R(TRealMatrix& OutMatrix)
{
  // the FFTW uses here 32b interface although it is internally 64b, it doesn't mind
  // since the size of 1 domain will never be bigger than 2^31 - however it it not a clear solution :)
  const int X   = static_cast<int> (OutMatrix.GetDimensionSizes().nx);
  const int Y   = static_cast<int> (OutMatrix.GetDimensionSizes().ny);
  const int Z   = static_cast<int> (OutMatrix.GetDimensionSizes().nz);

  // 1D FFT definition - over the Y axis
  const int  fft_rank = 1;
  fftw_iodim fft_dims[1];

  fft_dims[0].is = X * Y;
  fft_dims[0].n  = Z;
  fft_dims[0].os = X * Y;

  // GNU Compiler + FFTW
  #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    // How many definition - other dims
    const int  fft_howmany_rank = 2;
    fftw_iodim fft_howmany_dims[2];

    // Y dim
    fft_howmany_dims[0].is = X;
    fft_howmany_dims[0].n  = Y;
    fft_howmany_dims[0].os = X;

    // X dim
    fft_howmany_dims[1].is = 1;
    fft_howmany_dims[1].n  = X;
    fft_howmany_dims[1].os = 1;
  #endif

    // Intel Compiler + MKL
  #if (defined(__INTEL_COMPILER))
    const int fft_howmany_rank = 1;
    fftw_iodim fft_howmany_dims[1];

    // X dim
    fft_howmany_dims[0].is = 1;
    fft_howmany_dims[0].n  = X;
    fft_howmany_dims[0].os = 1;
  #endif

  fftw_plan_1DZ_C2R = fftwf_plan_guru_dft_c2r(fft_rank,                         // 1D FFT rank
                                              fft_dims,                         // 1D FFT dimensions of Y

                                              fft_howmany_rank,                 // how many in X and Z
                                              fft_howmany_dims,                 // Dims and strides in X and Z

                                              (fftwf_complex *) pMatrixData,    // input data
                                              OutMatrix.GetRawData(),           // output data
                                              TFFTWComplexMatrix_FFT_FLAG);     // flags

  if (!fftw_plan_1DZ_C2R)
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage, kErrFmtFftPlanNotCreated, "FFT_1DZ_C2R");
    throw runtime_error(ErrorMessage);
  }
}// end of Create_FFT_Plan_1DZ_C2R
//------------------------------------------------------------------------------


/**
 * Computer forward out-of place 3D Real-to-Complex FFT.
 * @param [in] InMatrix - Input Matrix
 * @throw runtime_error if the plan is not valid.
 */
void TFFTWComplexMatrix::Compute_FFT_3D_R2C(TRealMatrix & InMatrix)
{
  if (fftw_plan_3D_R2C)
  {
    fftwf_execute_dft_r2c(fftw_plan_3D_R2C, InMatrix.GetRawData(), (fftwf_complex *) pMatrixData);
  }
  else //error
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage, kErrFmtFftInvalidPlan, "FFT_3D_R2C");
    throw runtime_error(ErrorMessage);
  }
}// end of Compute_FFT_3D_R2C
//------------------------------------------------------------------------------

/**
 * Compute inverse out-of-place 3D Complex to Real FFT.
 * @param [out] OutMatrix
 * @throw runtime_error if the plan is not valid.
 */
void TFFTWComplexMatrix::Compute_FFT_3D_C2R(TRealMatrix & OutMatrix)
{
  if (fftw_plan_3D_C2R)
  {
    fftwf_execute_dft_c2r(fftw_plan_3D_C2R,(fftwf_complex *) pMatrixData, OutMatrix.GetRawData());
  }
  else // error
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage, kErrFmtFftInvalidPlan, "FFT_3D_C2R");
    throw runtime_error(ErrorMessage);
  }
}// end of Compute_FFT_3D_C2R
//------------------------------------------------------------------------------



/**
 * Compute 1D out-of-place Real-to-Complex FFT in the X dimension.
 * There are two versions of this routine for GCC+FFTW and ICPC + MKL, otherwise it will not build!
 * The FFTW version processes the whole matrix at one while the MKL slab by slab
 * @param [in] InMatrix
 * @throw runtime_error if the plan is not valid.
 */
void TFFTWComplexMatrix::Compute_FFT_1DX_R2C(TRealMatrix& InMatrix)
{
  if (fftw_plan_1DX_R2C)
  {
    // GNU Compiler + FFTW
    #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
      fftwf_execute_dft_r2c(fftw_plan_1DX_R2C,
                            InMatrix.GetRawData(),
                            (fftwf_complex *) pMatrixData);
    #endif

    // Intel Compiler + MKL
    #if (defined(__INTEL_COMPILER))
      const DimensionSizes Dims = TParameters::GetInstance()->GetFullDimensionSizes();
      for (size_t slab_id = 0; slab_id < Dims.Z; slab_id++)
      {
        fftwf_execute_dft_r2c(fftw_plan_1DX_R2C,
                              &InMatrix.GetRawData()[slab_id * Dims.X * Dims.Y],
                              (fftwf_complex *) &pMatrixData[slab_id * 2 * (Dims.X/2 + 1) * Dims.Y]);
      }
    #endif
  }
  else //error
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage, kErrFmtFftInvalidPlan, "FFT_1DX_R2C");
    throw runtime_error(ErrorMessage);
  }
}// end of Compute_FFT_1DX_R2C
//------------------------------------------------------------------------------

/**
 * Compute 1D out-of-place Real-to-Complex FFT in the Y dimension.
 * @param [in] InMatrix
 * @throw runtime_error if the plan is not valid.
 */
void TFFTWComplexMatrix::Compute_FFT_1DY_R2C(TRealMatrix& InMatrix)
{
  if (fftw_plan_1DY_R2C)
  {
    // GNU Compiler + FFTW
    #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
      fftwf_execute_dft_r2c(fftw_plan_1DY_R2C,
                            InMatrix.GetRawData(),
                            (fftwf_complex *) pMatrixData);
    #endif

    // Intel Compiler + MKL
    #if (defined(__INTEL_COMPILER))
      const DimensionSizes Dims = TParameters::GetInstance()->GetFullDimensionSizes();
      for (size_t slab_id = 0; slab_id < Dims.Z; slab_id++)
      {
        fftwf_execute_dft_r2c(fftw_plan_1DY_R2C,
                              &InMatrix.GetRawData()[slab_id * Dims.X * Dims.Y],
                              (fftwf_complex *) &pMatrixData[slab_id * Dims.X * 2 * (Dims.Y/2 + 1)]);
      }
    #endif
  }
  else //error
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage, kErrFmtFftInvalidPlan, "FFT_1DY_R2C");
    throw runtime_error(ErrorMessage);
  }
}// end of Compute_FFT_1DY_R2C
//------------------------------------------------------------------------------

/**
 * Compute 1D out-of-place Real-to-Complex FFT in the Z dimension.
 * @param [in] InMatrix
 * @throw runtime_error if the plan is not valid.
 */
void TFFTWComplexMatrix::Compute_FFT_1DZ_R2C(TRealMatrix& InMatrix)
{
  if (fftw_plan_1DZ_R2C)
  {
    // GNU Compiler + FFTW
    #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
      fftwf_execute_dft_r2c(fftw_plan_1DZ_R2C,
                            InMatrix.GetRawData(),
                            (fftwf_complex *) pMatrixData);
    #endif

    // Intel Compiler + MKL
    #if (defined(__INTEL_COMPILER))
      const DimensionSizes Dims = TParameters::GetInstance()->GetFullDimensionSizes();
      for (size_t slab_id = 0; slab_id < Dims.Y; slab_id++)
      {
        fftwf_execute_dft_r2c(fftw_plan_1DZ_R2C,
                              &InMatrix.GetRawData()[slab_id * Dims.X],
                              (fftwf_complex *) &pMatrixData[slab_id * 2 * Dims.X]);
      }
    #endif
  }
  else //error
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage, kErrFmtFftInvalidPlan, "FFT_1DZ_R2C");
    throw runtime_error(ErrorMessage);
  }
}// end of Compute_FFT_1DZ_R2C
//------------------------------------------------------------------------------


/**
 * Compute 1D out-of-place Complex-to-Real FFT in the X dimension.
 * @param [out] OutMatrix
 * @throw runtime_error if the plan is not valid.
 */
void TFFTWComplexMatrix::Compute_FFT_1DX_C2R(TRealMatrix& OutMatrix)
{
  if (fftw_plan_1DX_C2R)
  {
    // GNU Compiler + FFTW
    #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
      fftwf_execute_dft_c2r(fftw_plan_1DX_C2R,
                            (fftwf_complex *) pMatrixData,
                            OutMatrix.GetRawData());
    #endif

    // Intel Compiler + MKL
    #if (defined(__INTEL_COMPILER))
      const DimensionSizes Dims = TParameters::GetInstance()->GetFullDimensionSizes();
      for (size_t slab_id = 0; slab_id < Dims.Z; slab_id++)
      {
        fftwf_execute_dft_c2r(fftw_plan_1DX_C2R,
                              (fftwf_complex *) &pMatrixData[slab_id * 2 * (Dims.X/2 + 1) * Dims.Y],
                              &OutMatrix.GetRawData()[slab_id * Dims.X * Dims.Y]);
      }
    #endif
  }
  else //error
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage, kErrFmtFftInvalidPlan, "FFT_1DX_C2R");
    throw runtime_error(ErrorMessage);
  }
}// end of Compute_FFT_1DX_C2R
//------------------------------------------------------------------------------

/**
 * Compute 1D out-of-place Complex-to-Real FFT in the Y dimension.
 * @param [out] OutMatrix
 * @throw runtime_error if the plan is not valid.
 */
void TFFTWComplexMatrix::Compute_FFT_1DY_C2R(TRealMatrix& OutMatrix)
{
  if (fftw_plan_1DY_C2R)
  {
    // GNU Compiler + FFTW
    #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
      fftwf_execute_dft_c2r(fftw_plan_1DY_C2R,
                            (fftwf_complex *) pMatrixData,
                            OutMatrix.GetRawData());
    #endif

    // Intel Compiler + MKL
    #if (defined(__INTEL_COMPILER))
      const DimensionSizes Dims = TParameters::GetInstance()->GetFullDimensionSizes();
      for (size_t slab_id = 0; slab_id < Dims.Z; slab_id++)
      {
        fftwf_execute_dft_c2r(fftw_plan_1DY_C2R,
                              (fftwf_complex *) &pMatrixData[slab_id * Dims.X * 2 * (Dims.Y/2 + 1) ],
                              &OutMatrix.GetRawData()[slab_id * Dims.X * Dims.Y]);
      }
    #endif
  }
  else //error
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage, kErrFmtFftInvalidPlan, "FFT_1DY_C2R");
    throw runtime_error(ErrorMessage);
  }
}// end of Compute_FFT_1DY_C2R
//------------------------------------------------------------------------------

/**
 * Compute 1D out-of-place Complex-to-Real FFT in the Z dimension
 * @param [out] OutMatrix
 * @throw runtime_error if the plan is not valid.
 */
void TFFTWComplexMatrix::Compute_FFT_1DZ_C2R(TRealMatrix& OutMatrix)
{
  if (fftw_plan_1DZ_C2R)
  {
    // GNU Compiler + FFTW
    #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
      fftwf_execute_dft_c2r(fftw_plan_1DZ_C2R,
                            (fftwf_complex *) pMatrixData,
                            OutMatrix.GetRawData());
    #endif

    // Intel Compiler + MKL
    #if (defined(__INTEL_COMPILER))
      const DimensionSizes Dims = TParameters::GetInstance()->GetFullDimensionSizes();
      for (size_t slab_id = 0; slab_id < Dims.Y; slab_id++)
      {
        fftwf_execute_dft_c2r(fftw_plan_1DZ_C2R,
                              (fftwf_complex *) &pMatrixData[slab_id * 2 * Dims.X ],
                              &OutMatrix.GetRawData()[slab_id * Dims.X]);
      }
    #endif
  }
  else //error
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage, kErrFmtFftInvalidPlan, "FFT_1DZ_C2R");
    throw runtime_error(ErrorMessage);
  }
}// end of Compute_FFT_1DZ_C2R
//------------------------------------------------------------------------------


/**
 * Export wisdom to the file.
 * @warning This routine is supported only while compiling with the GNU C++ compiler.
 */
void TFFTWComplexMatrix::ExportWisdom()
{
  #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    int success = fftwf_export_wisdom_to_filename(GetWisdomFileName().c_str());

    if (success == 0)
    {
      fprintf(stderr,kErrFmtFftWisdomNotExported);
    }
  #endif
}// end of ExportWisdom
//------------------------------------------------------------------------------



/**
 * Import wisdom from the file.
 * @warning This routine is supported only while compiling with the GNU C++ compiler.
 */
void TFFTWComplexMatrix::ImportWisdom()
{
  #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    int success = fftwf_import_wisdom_from_filename(GetWisdomFileName().c_str());
    if (success == 0)
    {
      fprintf(stderr,ErrFmtFftWisdomNotImported);
    }
  #endif
}// end of Import wisdom
//------------------------------------------------------------------------------

/**
 * Delete stored wisdom (delete the file).
 */
void TFFTWComplexMatrix::DeleteStoredWisdom()
{
  std::remove(GetWisdomFileName().c_str());
}// end of DeleteStoredWisdom
//------------------------------------------------------------------------------


//----------------------------------------------------------------------------//
//                           Protected methods                                //
//----------------------------------------------------------------------------//

/**
 * Allocate Memory using fftwf_malloc function to ensure correct alignment.
 *
 */
void TFFTWComplexMatrix::AllocateMemory()
{
  /* No memory allocated before this function*/
  pMatrixData = (float *) fftwf_malloc(pTotalAllocatedElementCount * sizeof (float));

  if (!pMatrixData)
  {
    fprintf(stderr,kErrFmtNotEnoughMemory, "TFFTWComplexMatrix");
    throw bad_alloc();
  }

  // first touch
  #pragma omp parallel for schedule(static)
  for (size_t i=0; i<pTotalAllocatedElementCount; i++)
  {
    pMatrixData[i] = 0.0f;
  }
}// end of AllocateMemory
//------------------------------------------------------------------------------


 /**
  * Free memory using fftwf_free.
  */
void TFFTWComplexMatrix::FreeMemory()
{
  if (pMatrixData) fftwf_free( pMatrixData);
  pMatrixData = NULL;
 }// end of FreeMemory
 //-----------------------------------------------------------------------------

/**
 * Get Wisdom file name (derive it form the checkpoint filename).
 * @return the filename for wisdom.
 */
string TFFTWComplexMatrix::GetWisdomFileName()
{
  string FileName = TParameters::GetInstance()->GetCheckpointFileName();
  FileName.erase(FileName.find_last_of("."), string::npos);

  FileName.append(".");
  FileName.append(FFTW_Wisdom_FileName_Extension);

  return FileName;
}// end of GetWisdomFileName
//------------------------------------------------------------------------------



//----------------------------------------------------------------------------//
//                           Private methods                                  //
//----------------------------------------------------------------------------//