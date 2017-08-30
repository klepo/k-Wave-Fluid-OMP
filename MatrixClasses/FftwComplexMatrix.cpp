/**
 * @file        FftwComplexMatrix.cpp
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
 *              30 August    2017, 16:09 (revised)
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

#include <stdexcept>

#include <MatrixClasses/FftwComplexMatrix.h>
#include <MatrixClasses/RealMatrix.h>

#include <Parameters/Parameters.h>
#include <Logger/Logger.h>



//--------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------- Constants -----------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

const std::string FftwComplexMatrix::kFftWisdomFileExtension = "FFTW_Wisdom";
//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Public methods ---------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Constructor.
 */
FftwComplexMatrix::FftwComplexMatrix(const DimensionSizes& dimensionSizes)
  : ComplexMatrix(dimensionSizes),
    mR2CFftPlan3D(nullptr),  mC2RFftPlan3D(nullptr),
    mR2CFftPlan1DX(nullptr), mR2CFftPlan1DY(nullptr), mR2CFftPlan1DZ(nullptr),
    mC2RFftPlan1DX(nullptr), mC2RFftPlan1DY(nullptr), mC2RFftPlan1DZ(nullptr)
{
}// end of FftwComplexMatrix
//----------------------------------------------------------------------------------------------------------------------

/**
 * Destructor.
 */
FftwComplexMatrix::~FftwComplexMatrix()
{
  // free 3D plans
  if (mR2CFftPlan3D)  fftwf_destroy_plan(mR2CFftPlan3D);
  if (mC2RFftPlan3D)  fftwf_destroy_plan(mC2RFftPlan3D);

  //free 1D plans.
  if (mR2CFftPlan1DX) fftwf_destroy_plan(mR2CFftPlan1DX);
  if (mR2CFftPlan1DY) fftwf_destroy_plan(mR2CFftPlan1DY);
  if (mR2CFftPlan1DZ) fftwf_destroy_plan(mR2CFftPlan1DZ);

  if (mC2RFftPlan1DX) fftwf_destroy_plan(mC2RFftPlan1DX);
  if (mC2RFftPlan1DY) fftwf_destroy_plan(mC2RFftPlan1DY);
  if (mC2RFftPlan1DZ) fftwf_destroy_plan(mC2RFftPlan1DZ);

  mR2CFftPlan3D = nullptr;
  mC2RFftPlan3D = nullptr;

  mR2CFftPlan1DX = nullptr;
  mR2CFftPlan1DY = nullptr;
  mR2CFftPlan1DZ = nullptr;

  mC2RFftPlan1DX = nullptr;
  mC2RFftPlan1DY = nullptr;
  mC2RFftPlan1DZ = nullptr;

  // Memory will be freed by ~ComplexMatrix
  //freeMemory();

}// end of ~FftwComplexMatrix
//----------------------------------------------------------------------------------------------------------------------

/**
 * Create an FFTW plan for 3D Real-to-Complex.
 */
void FftwComplexMatrix::createR2CFftPlan3D(RealMatrix& inMatrix)
{
  mR2CFftPlan3D = fftwf_plan_dft_r2c_3d(inMatrix.getDimensionSizes().nz,
                                        inMatrix.getDimensionSizes().ny,
                                        inMatrix.getDimensionSizes().nx,
                                        inMatrix.getData(),
                                        reinterpret_cast<fftwf_complex*>(mData),
                                        kFftMeasureFlag);

  if (!mR2CFftPlan3D)
  {
    throw std::runtime_error(kErrFmtCreateR2CFftPlan3D);
  }
}// end of createR2CFftPlan3D
//----------------------------------------------------------------------------------------------------------------------

/**
 * Create an FFTW plan for 3D Complex-to-Real.
 */
void FftwComplexMatrix::createC2RFftPlan3D(RealMatrix& outMatrix)
{
  mC2RFftPlan3D = fftwf_plan_dft_c2r_3d(outMatrix.getDimensionSizes().nz,
                                        outMatrix.getDimensionSizes().ny,
                                        outMatrix.getDimensionSizes().nx,
                                        reinterpret_cast<fftwf_complex*>(mData),
                                        outMatrix.getData(),
                                        kFftMeasureFlag);

  if (!mC2RFftPlan3D)
  {
    throw std::runtime_error(kErrFmtCreateC2RFftPlan3D);
  }
}//end of createC2RFftPlan3D
//----------------------------------------------------------------------------------------------------------------------

/**
 * Create an FFTW plan for 1D Real-to-Complex in the x dimension.
 */
void FftwComplexMatrix::createR2CFftPlan1DX(RealMatrix& inMatrix)
{
  // the FFTW uses here 32b interface although it is internally 64b, it doesn't mind since the size of 1
  // domain will never be bigger than 2^31 - however it it not a clear solution :)
  const int nx  = static_cast<int> (inMatrix.getDimensionSizes().nx);
  const int ny  = static_cast<int> (inMatrix.getDimensionSizes().ny);
  const int nz  = static_cast<int> (inMatrix.getDimensionSizes().nz);
  const int nxR = ((nx / 2) + 1);

  // 1D FFT definition - over the x axis
  const int  rank = 1;
  fftw_iodim dims[1];

  dims[0].is = 1;
  dims[0].n  = nx;
  dims[0].os = 1;

  // GNU Compiler + FFTW does it all at once
  #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    // How FFTs we need to perform - Z * Y
    const int  howManyRank = 2;
    fftw_iodim howManyDims[2];

    // z dim
    howManyDims[0].is = nx * ny;
    howManyDims[0].n  = nz;
    howManyDims[0].os = nxR * ny;

    // y dim
    howManyDims[1].is = nx;
    howManyDims[1].n  = ny;
    howManyDims[1].os = nxR;
  #endif

  // Intel Compiler + MKL does it slab by slab
  #if (defined(__INTEL_COMPILER))
    const int  howManyRank = 1;
    fftw_iodim howManyDims[1];

    // y dim
    howManyDims[0].is = nx;
    howManyDims[0].n  = ny;
    howManyDims[0].os = nxR;
  #endif

  mR2CFftPlan1DX = fftwf_plan_guru_dft_r2c(rank,                                    // 1D FFT rank
                                           dims,                                    // 1D FFT dimensions of x
                                           howManyRank,                             // how many in y and z
                                           howManyDims,                             // Dims and strides in y and z
                                           inMatrix.getData(),                      // input data
                                           reinterpret_cast<fftwf_complex*>(mData), // output data
                                           kFftMeasureFlag);                        // flags

  if (!mR2CFftPlan1DX)
  {
    throw std::runtime_error(kErrFmtCreateR2CFftPlan1DX);
  }
}// end of createR2CFftPlan1DX
//----------------------------------------------------------------------------------------------------------------------


/**
 * Create an FFTW plan for 1D Real-to-Complex in the y dimension.
 */
void FftwComplexMatrix::createR2CFftPlan1DY(RealMatrix& inMatrix)
{
  // the FFTW uses here 32b interface although it is internally 64b, it doesn't mind
  // since the size of 1 domain will never be bigger than 2^31 - however it it not a clear solution :)
  const int nx  = static_cast<int>(inMatrix.getDimensionSizes().nx);
  const int ny  = static_cast<int>(inMatrix.getDimensionSizes().ny);
  const int nz  = static_cast<int>(inMatrix.getDimensionSizes().nz);
  const int nyR = ((ny / 2) + 1);

  // 1D FFT definition - over the ny axis
  const int  rank = 1;
  fftw_iodim dims[1];

  dims[0].is = nx;
  dims[0].n  = ny;
  dims[0].os = nx;

  // GNU Compiler + FFTW does it all at once
  #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    // How FFTs we need to perform - Z * X
    const int  howManyRank = 2;
    fftw_iodim howManyDims[2];

    // z dim
    howManyDims[0].is = nx * ny;
    howManyDims[0].n  = nz;
    howManyDims[0].os = nx * nyR ;

    // x dim
    howManyDims[1].is = 1;
    howManyDims[1].n  = nx;
    howManyDims[1].os = 1;
  #endif

   // Intel Compiler + MKL does it slab by slab
  #if (defined(__INTEL_COMPILER))
    const int  howManyRank = 1;
    fftw_iodim howManyDims[1];

    // x dim
    howManyDims[0].is = 1;
    howManyDims[0].n  = nx;
    howManyDims[0].os = 1;
  #endif

  mR2CFftPlan1DY = fftwf_plan_guru_dft_r2c(rank,                                    // 1D FFT rank
                                           dims,                                    // 1D FFT dimensions of y
                                           howManyRank,                             // how many in x and z
                                           howManyDims,                             // Dims and strides in x and z
                                           inMatrix.getData(),                      // input data
                                           reinterpret_cast<fftwf_complex*>(mData), // output data
                                           kFftMeasureFlag);                        // flags

  if (!mR2CFftPlan1DY)
  {
    throw std::runtime_error(kErrFmtCreateR2CFftPlan1DY);
  }
}// end of createR2CFftPlan1DY
//----------------------------------------------------------------------------------------------------------------------

/**
 * Create an FFTW plan for 1D Real-to-Complex in the z dimension.
 */
void FftwComplexMatrix::createR2CFftPlan1DZ(RealMatrix& inMatrix)
{
  // the FFTW uses here 32b interface although it is internally 64b, it doesn't mind
  // since the size of 1 domain will never be bigger than 2^31 - however it it not a clear solution :)
  const int nx = static_cast<int> (inMatrix.getDimensionSizes().nx);
  const int ny = static_cast<int> (inMatrix.getDimensionSizes().ny);
  const int nz = static_cast<int> (inMatrix.getDimensionSizes().nz);

  // 1D FFT definition - over the ny axis
  const int  rank = 1;
  fftw_iodim dims[1];

  dims[0].is = nx * ny;
  dims[0].n  = nz;
  dims[0].os = nx * ny;

  // GNU Compiler + FFTW
  #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    // How FFTs we need to perform - Y * X
    const int  howManyRank = 2;
    fftw_iodim howManyDims[2];

    // y dim
    howManyDims[0].is = nx;
    howManyDims[0].n  = ny;
    howManyDims[0].os = nx;

    // x dim
    howManyDims[1].is = 1;
    howManyDims[1].n  = nx;
    howManyDims[1].os = 1;
  #endif

  // Intel Compiler + MKL does it slab by slab
  #if (defined(__INTEL_COMPILER))
    const int  howManyRank = 1;
    fftw_iodim howManyDims[1];

    // x dim
    howManyDims[0].is = 1;
    howManyDims[0].n  = nx;
    howManyDims[0].os = 1;
  #endif

  mR2CFftPlan1DZ = fftwf_plan_guru_dft_r2c(rank,                                    // 1D FFT rank
                                           dims,                                    // 1D FFT dimensions of z
                                           howManyRank,                             // how many in x and y
                                           howManyDims,                             // Dims and strides in x and y
                                           inMatrix.getData(),                      // input data
                                           reinterpret_cast<fftwf_complex*>(mData), // output data
                                           kFftMeasureFlag);                        // flags

  if (!mR2CFftPlan1DZ)
  {
    throw std::runtime_error(kErrFmtCreateR2CFftPlan1DZ);
  }
}// end of createR2CFftPlan1DZ
//----------------------------------------------------------------------------------------------------------------------

/**
 * Create FFTW plan for Complex-to-Real in the x dimension.
 */
void FftwComplexMatrix::createC2RFftPlan1DX(RealMatrix& outMatrix)
{
  // the FFTW uses here 32b interface although it is internally 64b, it doesn't mind
  // since the size of 1 domain will never be bigger than 2^31 - however it it not a clear solution :)
  const int nx  = static_cast<int> (outMatrix.getDimensionSizes().nx);
  const int ny  = static_cast<int> (outMatrix.getDimensionSizes().ny);
  const int nz  = static_cast<int> (outMatrix.getDimensionSizes().nz);
  const int nxR = ((nx / 2) + 1);

  // 1D FFT definition - over the x axis
  const int  rank = 1;
  fftw_iodim dims[1];

  dims[0].is = 1;
  dims[0].n  = nx;
  dims[0].os = 1;

  // GNU Compiler + FFTW
  #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    // How FFTs we need to perform - Z * Y
    const int  howManyRank = 2;
    fftw_iodim howManyDims[2];

    // z dim
    howManyDims[0].is = nxR * ny;
    howManyDims[0].n  = nz;
    howManyDims[0].os = nx * ny;

    // y dim
    howManyDims[1].is = nxR;
    howManyDims[1].n  = ny;
    howManyDims[1].os = nx;
  #endif

  // Intel Compiler + MKL does it slab by slab
  #if (defined(__INTEL_COMPILER))
    const int  howManyRank = 1;
    fftw_iodim howManyDims[1];

    // y dim
    howManyDims[0].is = nxR;
    howManyDims[0].n  = ny;
    howManyDims[0].os = nx;
  #endif

  mC2RFftPlan1DX = fftwf_plan_guru_dft_c2r(rank,                                    // 1D FFT rank
                                           dims,                                    // 1D FFT dimensions of x
                                           howManyRank,                             // how many in y and z
                                           howManyDims,                             // Dims and strides in y and z
                                           reinterpret_cast<fftwf_complex*>(mData), // input data
                                           outMatrix.getData(),                     // output data
                                           kFftMeasureFlag);                        // flags

  if (!mC2RFftPlan1DX)
  {
    throw std::runtime_error(kErrFmtCreateC2RFftPlan1DX);
  }
}// end of createC2RFftPlan1DX
//----------------------------------------------------------------------------------------------------------------------

 /**
  * Create FFTW plan for Complex-to-Real in the y dimension.
  */
void FftwComplexMatrix::createC2RFftPlan1DY(RealMatrix& outMatrix)
{
  // the FFTW uses here 32b interface although it is internally 64b, it doesn't mind
  // since the size of 1 domain will never be bigger than 2^31 - however it it not a clear solution :)
  const int nx  = static_cast<int> (outMatrix.getDimensionSizes().nx);
  const int ny  = static_cast<int> (outMatrix.getDimensionSizes().ny);
  const int nz  = static_cast<int> (outMatrix.getDimensionSizes().nz);
  const int nyR = ((ny / 2) + 1);

  // 1D FFT definition - over the y axis
  const int  rank = 1;

  fftw_iodim dims[1];
  dims[0].is = nx;
  dims[0].n  = ny;
  dims[0].os = nx;

  // GNU Compiler + FFTW
  #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    // How FFTs we need to perform - Z * X
    const int  howManyRank = 2;
    fftw_iodim howManyDims[2];

    // z dim
    howManyDims[0].is = nx * nyR;
    howManyDims[0].n  = nz;
    howManyDims[0].os = nx * ny;

    // x dim
    howManyDims[1].is = 1;
    howManyDims[1].n  = nx;
    howManyDims[1].os = 1;
  #endif

  // Intel Compiler + MKL does it slab by slab
  #if (defined(__INTEL_COMPILER))
    const int  howManyRank = 1;
    fftw_iodim howManyDims[1];

    // x dim
    howManyDims[0].is = 1;
    howManyDims[0].n  = nx;
    howManyDims[0].os = 1;
  #endif

  mC2RFftPlan1DY = fftwf_plan_guru_dft_c2r(rank,                                    // 1D FFT rank
                                           dims,                                    // 1D FFT dimensions of y
                                           howManyRank,                             // how many in x and z
                                           howManyDims,                             // Dims and strides in x and z
                                           reinterpret_cast<fftwf_complex*>(mData), // input data
                                           outMatrix.getData(),                     // output data
                                           kFftMeasureFlag);                        // flags

  if (!mC2RFftPlan1DY)
  {
    throw std::runtime_error(kErrFmtCreateC2RFftPlan1DY);
  }
}// end of createC2RFftPlan1DY
//----------------------------------------------------------------------------------------------------------------------

/**
 * Create FFTW plan for Complex-to-Real in the z dimension.
 */
void FftwComplexMatrix::createC2RFftPlan1DZ(RealMatrix& outMatrix)
{
  // the FFTW uses here 32b interface although it is internally 64b, it doesn't mind
  // since the size of 1 domain will never be bigger than 2^31 - however it it not a clear solution :)
  const int nx   = static_cast<int> (outMatrix.getDimensionSizes().nx);
  const int ny   = static_cast<int> (outMatrix.getDimensionSizes().ny);
  const int nz   = static_cast<int> (outMatrix.getDimensionSizes().nz);

  // 1D FFT definition - over the z axis
  const int  rank = 1;
  fftw_iodim dims[1];

  dims[0].is = nx * ny;
  dims[0].n  = nz;
  dims[0].os = nx * ny;

  // GNU Compiler + FFTW
  #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    // How FFTs we need to perform - Y * X
    const int  howManyRank = 2;
    fftw_iodim howManyDims[2];

    // y dim
    howManyDims[0].is = nx;
    howManyDims[0].n  = ny;
    howManyDims[0].os = nx;

    // x dim
    howManyDims[1].is = 1;
    howManyDims[1].n  = nx;
    howManyDims[1].os = 1;
  #endif

  // Intel Compiler + MKL does it slab by slab
  #if (defined(__INTEL_COMPILER))
    const int howManyRank = 1;
    fftw_iodim howManyDims[1];

    // x dim
    howManyDims[0].is = 1;
    howManyDims[0].n  = nx;
    howManyDims[0].os = 1;
  #endif

  mC2RFftPlan1DZ = fftwf_plan_guru_dft_c2r(rank,                                    // 1D FFT rank
                                           dims,                                    // 1D FFT dimensions of z
                                           howManyRank,                             // how many in x and y
                                           howManyDims,                             // Dims and strides in x and y
                                           reinterpret_cast<fftwf_complex*>(mData), // input data
                                           outMatrix.getData(),                     // output data
                                           kFftMeasureFlag);                        // flags

  if (!mC2RFftPlan1DZ)
  {
    throw std::runtime_error(kErrFmtCreateC2RFftPlan1DX);
  }
}// end of createC2RFftPlan1DZ
//----------------------------------------------------------------------------------------------------------------------

/**
 * Computer forward out-of place 3D Real-to-Complex FFT.
 */
void FftwComplexMatrix::computeR2CFft3D(RealMatrix& inMatrix)
{
  if (mR2CFftPlan3D)
  {
    fftwf_execute_dft_r2c(mR2CFftPlan3D, inMatrix.getData(), reinterpret_cast<fftwf_complex*>(mData));
  }
  else //error
  {
    throw std::runtime_error(kErrFmtExecuteR2CFftPlan3D);
  }
}// end of computeR2CFft3D
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute inverse out-of-place 3D Complex to Real FFT.
 */
void FftwComplexMatrix::computeC2RFft3D(RealMatrix & outMatrix)
{
  if (mC2RFftPlan3D)
  {
    fftwf_execute_dft_c2r(mC2RFftPlan3D, reinterpret_cast<fftwf_complex*>(mData), outMatrix.getData());
  }
  else // error
  {
    throw std::runtime_error(kErrFmtExecuteC2RFftPlan3D);
  }
}// end of computeC2RFft3D
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute 1D out-of-place Real-to-Complex FFT in the x dimension.
 */
void FftwComplexMatrix::computeR2CFft1DX(RealMatrix& inMatrix)
{
  if (mR2CFftPlan1DX)
  {
    // GNU Compiler + FFTW
    #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
      fftwf_execute_dft_r2c(mR2CFftPlan1DX,
                            inMatrix.getData(),
                            reinterpret_cast<fftwf_complex*>(mData));
    #endif

    // Intel Compiler + MKL
    #if (defined(__INTEL_COMPILER))
      const DimensionSizes dims = Parameters::getInstance().getFullDimensionSizes();
      for (size_t slab_id = 0; slab_id < dims.nz; slab_id++)
      {
        fftwf_execute_dft_r2c(mR2CFftPlan1DX,
                              &inMatrix.getData()[slab_id * dims.nx * dims.ny],
                              (fftwf_complex *) &mData[slab_id * 2 * (dims.nx / 2 + 1) * dims.ny]);
      }
    #endif
  }
  else //error
  {
    throw std::runtime_error(kErrFmtExecuteR2CFftPlan1DX);
  }
}// end of computeR2CFft1DX
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute 1D out-of-place Real-to-Complex FFT in the Y dimension.
 */
void FftwComplexMatrix::computeR2CFft1DY(RealMatrix& inMatrix)
{
  if (mR2CFftPlan1DY)
  {
    // GNU Compiler + FFTW
    #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
      fftwf_execute_dft_r2c(mR2CFftPlan1DY,
                            inMatrix.getData(),
                            reinterpret_cast<fftwf_complex*>(mData));
    #endif

    // Intel Compiler + MKL
    #if (defined(__INTEL_COMPILER))
      const DimensionSizes dims = Parameters::getInstance().getFullDimensionSizes();
      for (size_t slab_id = 0; slab_id < dims.nz; slab_id++)
      {
        fftwf_execute_dft_r2c(mR2CFftPlan1DY,
                              &inMatrix.getData()[slab_id * dims.nx * dims.ny],
                              (fftwf_complex *) &mData[slab_id * dims.nx * 2 * (dims.ny / 2 + 1)]);
      }
    #endif
  }
  else //error
  {
    throw std::runtime_error(kErrFmtExecuteR2CFftPlan1DY);
  }
}// end of computeR2CFft1DY
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute 1D out-of-place Real-to-Complex FFT in the z dimension.
 */
void FftwComplexMatrix::computeR2CFft1DZ(RealMatrix& inMatrix)
{
  if (mR2CFftPlan1DZ)
  {
    // GNU Compiler + FFTW
    #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
      fftwf_execute_dft_r2c(mR2CFftPlan1DZ,
                            inMatrix.getData(),
                            reinterpret_cast<fftwf_complex*>(mData));
    #endif

    // Intel Compiler + MKL
    #if (defined(__INTEL_COMPILER))
      const DimensionSizes dims = Parameters::getInstance().getFullDimensionSizes();
      for (size_t slab_id = 0; slab_id < dims.ny; slab_id++)
      {
        fftwf_execute_dft_r2c(mR2CFftPlan1DZ,
                              &inMatrix.getData()[slab_id * dims.nx],
                              (fftwf_complex *) &mData[slab_id * 2 * dims.nx]);
      }
    #endif
  }
  else //error
  {
    throw std::runtime_error(kErrFmtExecuteR2CFftPlan1DZ);
  }
}// end of computeR2CFft1DZ
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute 1D out-of-place Complex-to-Real FFT in the x dimension.
 */
void FftwComplexMatrix::computeC2RFft1DX(RealMatrix& outMatrix)
{
  if (mC2RFftPlan1DX)
  {
    // GNU Compiler + FFTW
    #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
      fftwf_execute_dft_c2r(mC2RFftPlan1DX,
                            reinterpret_cast<fftwf_complex*>(mData),
                            outMatrix.getData());
    #endif

    // Intel Compiler + MKL
    #if (defined(__INTEL_COMPILER))
      const DimensionSizes dims = Parameters::getInstance().getFullDimensionSizes();
      for (size_t slab_id = 0; slab_id < dims.nz; slab_id++)
      {
        fftwf_execute_dft_c2r(mC2RFftPlan1DX,
                              (fftwf_complex *) &mData[slab_id * 2 * (dims.nx / 2 + 1) * dims.ny],
                              &outMatrix.getData()[slab_id * dims.nx * dims.ny]);
      }
    #endif
  }
  else //error
  {
    throw std::runtime_error(kErrFmtExecuteC2RFftPlan1DX);
  }
}// end of computeR2CFft1DX
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute 1D out-of-place Complex-to-Real FFT in the y dimension.
 */
void FftwComplexMatrix::computeC2RFft1DY(RealMatrix& outMatrix)
{
  if (mC2RFftPlan1DY)
  {
    // GNU Compiler + FFTW
    #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
      fftwf_execute_dft_c2r(mC2RFftPlan1DY,
                            reinterpret_cast<fftwf_complex*>(mData),
                            outMatrix.getData());
    #endif

    // Intel Compiler + MKL
    #if (defined(__INTEL_COMPILER))
      const DimensionSizes dims = Parameters::getInstance().getFullDimensionSizes();
      for (size_t slab_id = 0; slab_id < dims.nz; slab_id++)
      {
        fftwf_execute_dft_c2r(mC2RFftPlan1DY,
                              (fftwf_complex *) &mData[slab_id * dims.nx * 2 * (dims.ny / 2 + 1)],
                              &outMatrix.getData()[slab_id * dims.nx * dims.ny]);
      }
    #endif
  }
  else //error
  {
    throw std::runtime_error(kErrFmtExecuteC2RFftPlan1DY);
  }
}// end of computeR2CFft1DY
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute 1D out-of-place Complex-to-Real FFT in the Z dimension
 */
void FftwComplexMatrix::computeC2RFft1DZ(RealMatrix& outMatrix)
{
  if (mC2RFftPlan1DZ)
  {
    // GNU Compiler + FFTW
    #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
      fftwf_execute_dft_c2r(mC2RFftPlan1DZ,
                            reinterpret_cast<fftwf_complex*>(mData),
                            outMatrix.getData());
    #endif

    // Intel Compiler + MKL
    #if (defined(__INTEL_COMPILER))
      const DimensionSizes dims = Parameters::getInstance().getFullDimensionSizes();
      for (size_t slab_id = 0; slab_id < dims.ny; slab_id++)
      {
        fftwf_execute_dft_c2r(mC2RFftPlan1DZ,
                              (fftwf_complex *) &mData[slab_id * 2 * dims.nx ],
                              &outMatrix.getData()[slab_id * dims.nx]);
      }
    #endif
  }
  else //error
  {
    throw std::runtime_error(kErrFmtExecuteC2RFftPlan1DZ);
  }
}// end of computeR2CFft1DZ
//----------------------------------------------------------------------------------------------------------------------

/**
 * Export wisdom to the file.
 */
void FftwComplexMatrix::exportWisdom()
{
  #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    int success = fftwf_export_wisdom_to_filename(getWisdomFileName().c_str());
    if (success == 0)
    {
      throw std::runtime_error(kErrFmtFftWisdomNotExported);
    }
  #endif
}// end of exportWisdom
//----------------------------------------------------------------------------------------------------------------------

/**
 * Import wisdom from the file.
 */
void FftwComplexMatrix::importWisdom()
{
  #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    int success = fftwf_import_wisdom_from_filename(getWisdomFileName().c_str());
    if (success == 0)
    {
      // print out a warning!
      throw std::runtime_error(ErrFmtFftWisdomNotImported);
    }
  #endif
}// end of importWisdom
//----------------------------------------------------------------------------------------------------------------------

/**
 * Delete stored wisdom (delete the file).
 */
void FftwComplexMatrix::deleteStoredWisdom()
{
  std::remove(getWisdomFileName().c_str());
}// end of deleteStoredWisdom
//----------------------------------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------ Protected methods -------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Allocate Memory using fftwf_malloc function to ensure correct alignment.
 */
void FftwComplexMatrix::allocateMemory()
{
  /* No memory allocated before this function*/
  mData = (float *) fftwf_malloc(mCapacity * sizeof (float));

  if (!mData)
  {
    throw std::bad_alloc();
  }

  // first touch
  #pragma omp parallel for simd schedule(static)
  for (size_t i = 0; i < mCapacity; i++)
  {
    mData[i] = 0.0f;
  }
}// end of allocateMemory
//----------------------------------------------------------------------------------------------------------------------

 /**
  * Free memory using fftwf_free.
  */
void FftwComplexMatrix::freeMemory()
{
  if (mData) fftwf_free(mData);
  mData = nullptr;
 }// end of freeMemory
 //---------------------------------------------------------------------------------------------------------------------

/**
 * Get Wisdom file name (derive it form the checkpoint filename).
 */
std::string FftwComplexMatrix::getWisdomFileName()
{
  std::string fileName = Parameters::getInstance().getCheckpointFileName();
  fileName.erase(fileName.find_last_of("."), std::string::npos);

  fileName.append(".");
  fileName.append(kFftWisdomFileExtension);

  return fileName;
}// end of getWisdomFileName
//----------------------------------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------ Private methods ---------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//
