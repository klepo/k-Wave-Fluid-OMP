/**
 * @file        VelocityMatrix.cpp
 * @author      Jiri Jaros
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The implementation file containing the particle velocity matrix.
 *
 * @version     kspaceFirstOrder3D 2.16
 *
 * @date        28 July      2011, 11:37 (created) \n
 *              28 August    2017, 14:45 (revised)
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


#include <MatrixClasses/VelocityMatrix.h>
#include <MatrixClasses/FftwComplexMatrix.h>

//--------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------- Constants -----------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//


//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Public methods ---------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Compute acoustic velocity for initial pressure problem.
 */
void VelocityMatrix::computeInitialVelocity(const RealMatrix&  dtRho0Sgxyz,
                                            FftwComplexMatrix& fftTemp)
{
  fftTemp.computeC2RFft3D(*this);

  const float divider = 1.0f / (2.0f * static_cast<float>(mSize));

  #pragma omp parallel for schedule (static) firstprivate(divider)
  for (size_t i = 0; i < mSize; i++)
  {
    mData[i] *= dtRho0Sgxyz[i] * divider;
  }
}// end of computeInitialVelocity
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute velocity for the initial pressure problem, homogeneous medium, uniform grid.
 */
void VelocityMatrix::computeInitialVelocityHomogeneousUniform(const float dtRho0Sgxyz,
                                                              FftwComplexMatrix& fftTemp)
{
  fftTemp.computeC2RFft3D(*this);

  const float divider = 1.0f / (2.0f * static_cast<float>(mSize)) * dtRho0Sgxyz;

  #pragma omp parallel for schedule (static) firstprivate(divider)
  for (size_t i = 0; i < mSize; i++)
  {
    mData[i] *= divider;
  }
}// end of computeInitialVelocityHomogeneousUniform
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute acoustic velocity for initial pressure problem, homogenous medium, nonuniform grid, x direction.
 */
void VelocityMatrix::computeInitialVelocityXHomogeneousNonuniform(const float        dtRho0Sgx,
                                                                  const RealMatrix&  dxudxnSgx,
                                                                  FftwComplexMatrix& fftTemp)
{
  fftTemp.computeC2RFft3D(*this);

  const float divider = 1.0f / (2.0f * static_cast<float>(mSize)) * dtRho0Sgx;

  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    size_t i = z * mDimensionSizes.ny * mDimensionSizes.nx;
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        mData[i] *=  divider * dxudxnSgx[x];
        i++;
      } // x
    } // y
  } // z
}// end of computeInitialVelocityXHomogeneousNonuniform
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute acoustic velocity for initial pressure problem, homogenous medium, nonuniform grid, x direction.
 */
void VelocityMatrix::computeInitialVelocityYHomogeneousNonuniform(const float        dtRho0Sgy,
                                                                  const RealMatrix&  dyudynSgy,
                                                                  FftwComplexMatrix& fftTemp)
{
  fftTemp.computeC2RFft3D(*this);

  const float divider = 1.0f / (2.0f * static_cast<float>(mSize)) * dtRho0Sgy;

  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    register size_t i = z * mDimensionSizes.ny * mDimensionSizes.nx;
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      const float eDyudynSgy = dyudynSgy[y] * divider;
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        mData[i] *= eDyudynSgy;
        i++;
      } // x
    } // y
  } // z
}// end of computeInitialVelocityYHomogeneousNonuniform
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute acoustic velocity for initial pressure problem, homogenous medium, nonuniform grid, z direction.
 */
void VelocityMatrix::computeInitialVelocityZHomogeneousNonuniform(const float dtRho0Sgz,
                                                                  const  RealMatrix& dzudznSgz,
                                                                  FftwComplexMatrix& fftTemp)
{
  fftTemp.computeC2RFft3D(*this);

  const float divider = 1.0f / (2.0f * static_cast<float>(mSize)) * dtRho0Sgz;

  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    size_t i = z * mDimensionSizes.ny * mDimensionSizes.nx;
    const float eDzudznSgz = dzudznSgz[z] * divider;

    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        mData[i] *= eDzudznSgz;
        i++;
      } // x
    } // y
  } // z
 }// end of computeInitialVelocityZHomogeneousNonuniform
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute acoustic velocity for heterogeneous medium and a uniform grid, x direction.
 */
void VelocityMatrix::computeVelocityX(const RealMatrix& ifftX,
                                      const RealMatrix& dtRho0Sgx,
                                      const RealMatrix& pmlX)
{
  const float divider = 1.0f / static_cast<float>(mSize);

  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    register size_t i = z * mSlabSize;
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        register float eData = mData[i];

        //FFT_p.ElementMultiplyMatrices(dt_rho0);
        const float eIfftX = divider * ifftX[i] * dtRho0Sgx[i];

        //BSXElementRealMultiply_1D_X(pml);
        eData *= pmlX[x];

        //ElementSubMatrices(FFT_p);
        eData -= eIfftX;

        //BSXElementRealMultiply_1D_X(pml);
        mData[i] = eData * pmlX[x];

        i++;
      } // x
    } // y
  } // z
}// end of computeVelocityX
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute acoustic velocity for homogeneous medium and a uniform grid, x direction.
 */
void VelocityMatrix::computeVelocityXHomogeneousUniform(const RealMatrix& ifftX,
                                                        const float       dtRho0,
                                                        const RealMatrix& pmlX)
{
  const float divider = dtRho0 / static_cast<float>(mSize);

  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    register size_t i = z * mSlabSize;
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        register float eData = mData[i];

        //FFT_p.ElementMultiplyMatrices(dt_rho0);
        const float eIfftX = divider * ifftX[i];

        //BSXElementRealMultiply_1D_X(pml);
        eData *= pmlX[x];

        //ElementSubMatrices(FFT_p);
        eData -= eIfftX;

        //BSXElementRealMultiply_1D_X(pml);
        mData[i] = eData * pmlX[x];

        i++;
      } // x
    } // y
  } // z
}// end of computeVelocityXHomogeneousUniform
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute acoustic velocity for homogenous medium and nonuniform grid, x direction.
 */
void VelocityMatrix::computeVelocityXHomogeneousNonuniform(const RealMatrix& ifftX,
                                                           const float       dtRho0,
                                                           const RealMatrix& dxudxnSgx,
                                                           const RealMatrix& pmlX)
{
  const float divider = dtRho0 / static_cast<float>(mSize);

  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    register size_t i = z * mSlabSize;
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        register float eData = mData[i];

        //FFT_p.ElementMultiplyMatrices(dt_rho0);
        const float eIfftX = (divider * dxudxnSgx[x]) * ifftX[i];

        //BSXElementRealMultiply_1D_X(pml);
        eData *= pmlX[x];

        //ElementSubMatrices(FFT_p);
        eData -= eIfftX;

        //BSXElementRealMultiply_1D_X(pml);
        mData[i] = eData * pmlX[x];

        i++;
      } // x
    } // y
  } // z
}// end of computeVelocityXHomogeneousNonuniform
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute acoustic velocity for heterogeneous medium and a uniform grid, y direction.
 */
void VelocityMatrix::computeVelocityY(const RealMatrix& ifftY,
                                      const RealMatrix& dtRho0Sgy,
                                      const RealMatrix& pmlY)
{
  const float divider = 1.0f / static_cast<float>(mSize);

  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    size_t i = z * mSlabSize;
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      const float ePmlY = pmlY[y];
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        register float eData = mData[i];

        //FFT_p.ElementMultiplyMatrices(dt_rho0);
        const float eIfftY = divider * ifftY[i] * dtRho0Sgy[i];

        //BSXElementRealMultiply_1D_X(pml);
        eData *= ePmlY;

        //ElementSubMatrices(FFT_p);
        eData -= eIfftY;

        //BSXElementRealMultiply_1D_X(pml);
        mData[i] = eData * ePmlY;

        i++;
      } // x
    } // y
  } // z
}// end of computeVelocityY
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute acoustic velocity for homogeneous medium and a uniform grid, y direction.
 */
void VelocityMatrix::computeVelocityYHomogeneousUniform(const RealMatrix& ifftY,
                                                        const float       dtRho0,
                                                        const RealMatrix& pmlY)
{
  const float divider = dtRho0 / static_cast<float>(mSize);

  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    size_t i = z * mSlabSize;
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      const float ePmlY = pmlY[y];
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        register float eData = mData[i];

        //FFT_p.ElementMultiplyMatrices(dt_rho0);
        const float eIfftY = divider * ifftY[i];

        //BSXElementRealMultiply_1D_X(pml);
        eData *= ePmlY;

        //ElementSubMatrices(FFT_p);
        eData -= eIfftY;

        //BSXElementRealMultiply_1D_X(pml);
        mData[i] = eData * ePmlY;

        i++;
      } // x
    } // y
  } // z
}// end of computeVelocityYHomogeneousUniform
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute acoustic velocity for homogenous medium and non-uniform grid, y direction.
 */
void VelocityMatrix::computeVelocityYHomogeneousNonuniform(const RealMatrix& ifftY,
                                                           const float       dtRho0,
                                                           const RealMatrix& dyudynSgy,
                                                           const RealMatrix& pmlY)
{
  const float divider = dtRho0 / static_cast<float>(mSize);

  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    size_t i = z * mSlabSize;
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      const float ePmlY      = pmlY[y];
      const float eDyudynSgy = dyudynSgy[y];
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        register float eData = mData[i];

        //FFT_p.ElementMultiplyMatrices(dt_rho0);
        const float eIfftY = (divider * eDyudynSgy) * ifftY[i];

        //BSXElementRealMultiply_1D_X(pml);
        eData *= ePmlY;

        //ElementSubMatrices(FFT_p);
        eData -= eIfftY;

        //BSXElementRealMultiply_1D_X(pml);
        mData[i] = eData * ePmlY;

        i++;
      } // x
    } // y
  } // z
}// end of computeVelocityYHomogeneousNonuniform
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute acoustic velocity for heterogeneous medium and a uniform grid, z direction.
 */
void VelocityMatrix::computeVelocityZ(const RealMatrix& ifftZ,
                                      const RealMatrix& dtRho0Sgz,
                                      const RealMatrix& pmlZ)
{
  const float divider = 1.0f / static_cast<float>(mSize);

  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    size_t i = z * mSlabSize;
    const float ePmlZ = pmlZ[z];
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        register float eData = mData[i];

        //FFT_p.ElementMultiplyMatrices(dt_rho0);
        const float eIfftZ = divider * ifftZ[i] * dtRho0Sgz[i];

        //BSXElementRealMultiply_1D_X(pml);
        eData *= ePmlZ;

        //ElementSubMatrices(FFT_p);
        eData -= eIfftZ;

        //BSXElementRealMultiply_1D_X(pml);
        mData[i] = eData * ePmlZ;

        i++;
      } // x
    } // y
  } // z
}// end of computeVelocityZ
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute acoustic velocity for homogeneous medium and a uniform grid, z direction.
 */
void VelocityMatrix::computeVelocityZHomogeneousUniform(const RealMatrix& ifftZ,
                                                        const float       dtRho0,
                                                        const RealMatrix& pmlZ)
{
  const float divider = dtRho0 / static_cast<float>(mSize);

  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    size_t i = z* mSlabSize;
    const float ePmlZ = pmlZ[z];
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        register float eData = mData[i];

        //FFT_p.ElementMultiplyMatrices(dt_rho0);
        const float eIfftZ = divider * ifftZ[i];

        //BSXElementRealMultiply_1D_X(abc);
        eData *= ePmlZ;

        //ElementSubMatrices(FFT_p);
        eData -= eIfftZ;

        //BSXElementRealMultiply_1D_X(abc);

        mData[i] = eData * ePmlZ;

        i++;
      } // x
    } // y
  } // z
}// end of computeVelocityZHomogeneousUniform
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute acoustic velocity for homogenous medium and non-uniform grid, z direction.
 */
void VelocityMatrix::computeVelocityZHomogeneousNonuniform(const RealMatrix& ifftZ,
                                                           const float       dtRho0,
                                                           const RealMatrix& dzudznSgz,
                                                           const RealMatrix& pmlZ)
{
  const float divider = dtRho0 / static_cast<float>(mSize);

  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    size_t i = z * mSlabSize;
    const float ePmlZ = pmlZ[z];
    const float eDzudznSgz_data = dzudznSgz[z];

    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        register float eData = mData[i];

        //FFT_p.ElementMultiplyMatrices(dt_rho0);
        const float eIfftZ = (divider * eDzudznSgz_data) * ifftZ[i];

        //BSXElementRealMultiply_1D_X(abc);
        eData *= ePmlZ;

        //ElementSubMatrices(FFT_p);
        eData -= eIfftZ;

        //BSXElementRealMultiply_1D_X(abc);

        mData[i] = eData * ePmlZ;

        i++;
      } // x
    } // y
  } // z
}// end of computeVelocityZHomogeneousNonuniform
//----------------------------------------------------------------------------------------------------------------------

/**
 * Add transducer data source to velocity x component.
 */
void VelocityMatrix::addTransducerSource(const IndexMatrix& velocitySourceIndex,
                                         const RealMatrix&  transducerSourceInput,
                                         const IndexMatrix& delayMask,
                                         const size_t       timeIndex)
{
  const size_t sourceSize = velocitySourceIndex.size();
  #pragma omp parallel for schedule(static) firstprivate(timeIndex, sourceSize) if (velocitySourceIndex.size() > 16384)
  for (size_t i = 0; i < sourceSize; i++)
  {
    mData[velocitySourceIndex[i]] += transducerSourceInput[delayMask[i] + timeIndex];
  }
}// end of addTransducerSource
//----------------------------------------------------------------------------------------------------------------------

/**
 * Add in velocity source terms.
 */
void VelocityMatrix::addVelocitySource(const RealMatrix & velocitySourceInput,
                                       const IndexMatrix& velocitySourceIndex,
                                       const size_t       timeIndex,
                                       const size_t       velocitySourceMode,
                                       const size_t       velocitySourceMany)
{

  const size_t sourceSize = velocitySourceIndex.size();
  const size_t index2D    = (velocitySourceMany != 0) ? timeIndex * sourceSize : timeIndex;

  if (velocitySourceMode == 0)
  {
    #pragma omp parallel for if (sourceSize > 16384)
    for (size_t i = 0; i < sourceSize; i++)
    {
      const size_t signalIndex = (velocitySourceMany != 0) ? index2D + i : index2D;
      mData[velocitySourceIndex[i]] = velocitySourceInput[signalIndex];
    }
  }// end of Dirichlet

  if (velocitySourceMode == 1)
  {
    #pragma omp parallel for if (sourceSize > 16384)
    for (size_t i = 0; i < sourceSize; i++)
    {
      const size_t signalIndex = (velocitySourceMany != 0) ? index2D + i : index2D;
      mData[velocitySourceIndex[i]] += velocitySourceInput[signalIndex];
    }
  }// end of add
}// end of addVelocitySource
//----------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------ Protected methods -------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//


//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Private methods --------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//
