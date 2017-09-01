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
 *              01 September 2017, 17:17 (revised)
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

  const float  divider      = 1.0f / (2.0f * static_cast<float>(mSize));
  const float* pDtRho0Sgxyz = dtRho0Sgxyz.getData();

  #pragma omp parallel for simd schedule (static) aligned(pDtRho0Sgxyz)
  for (size_t i = 0; i < mSize; i++)
  {
    mData[i] *= pDtRho0Sgxyz[i] * divider;
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

  #pragma omp parallel for simd schedule (static)
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
  const float* pDxudxnSgx = dxudxnSgx.getData();

  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, mDimensionSizes);
        mData[i] *=  divider * pDxudxnSgx[x];
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
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      const float eDyudynSgy = dyudynSgy[y] * divider;
      #pragma omp simd
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, mDimensionSizes);
        mData[i] *= eDyudynSgy;
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
    const float eDzudznSgz = dzudznSgz[z] * divider;
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, mDimensionSizes);
        mData[i] *= eDzudznSgz;
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
  const float divider     = 1.0f / static_cast<float>(mSize);

  const float* dIfftX     = ifftX.getData();
  const float* dDtRho0Sgx = dtRho0Sgx.getData();
  const float* dPmlX      = pmlX.getData();

  #pragma omp for schedule(static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, mDimensionSizes);

        const float eIfftX = divider * dIfftX[i] * dDtRho0Sgx[i];
        const float ePmlX  = dPmlX[x];

        mData[i] = (mData[i] * ePmlX - eIfftX) * ePmlX;
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

  const float* eIfftX = ifftX.getData();
  const float* ePmlX  = pmlX.getData();

  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, mDimensionSizes);

        mData[i] = (mData[i] * ePmlX[x] - divider * eIfftX[i]) * ePmlX[x];
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
  const float divider     = dtRho0 / static_cast<float>(mSize);

  const float* eIfftX     = ifftX.getData();
  const float* eDxudxnSgx = dxudxnSgx.getData();
  const float* ePmlX      = pmlX.getData();


  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, mDimensionSizes);

        mData[i] = (mData[i] * ePmlX[x] - (divider * eDxudxnSgx[x] * eIfftX[i])) * ePmlX[x];
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

  const float* dIfftY = ifftY.getData();
  const float* dDtRho0Sgy = dtRho0Sgy.getData();

  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      const float ePmlY = pmlY[y];
      #pragma omp simd
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, mDimensionSizes);
        const float eIfftY = divider * dIfftY[i] * dDtRho0Sgy[i];

        mData[i] = (mData[i] * ePmlY - eIfftY) * ePmlY;
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

  const float* eIfftY = ifftY.getData();

  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      const float ePmlY = pmlY[y];
      #pragma omp simd
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, mDimensionSizes);

        mData[i] = (mData[i] * ePmlY - divider * eIfftY[i]) * ePmlY;
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
  const float  divider = dtRho0 / static_cast<float>(mSize);
  const float* eIfftY  = ifftY.getData();

  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      const float ePmlY      = pmlY[y];
      const float eDyudynSgy = dyudynSgy[y];
      #pragma omp simd
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, mDimensionSizes);

        mData[i] = (mData[i] * ePmlY - (divider * eDyudynSgy * eIfftY[i])) * ePmlY;
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

  const float* dIfftZ     = ifftZ.getData();
  const float* dDtRho0Sgz = dtRho0Sgz.getData();

  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    const float ePmlZ = pmlZ[z];
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, mDimensionSizes);
        const float eIfftZ = divider * dIfftZ[i] * dDtRho0Sgz[i];

        mData[i] = (mData[i] * ePmlZ - eIfftZ) * ePmlZ;
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
  const float* eIfftZ = ifftZ.getData();

  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    const float ePmlZ = pmlZ[z];
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, mDimensionSizes);

        mData[i] = (mData[i] * ePmlZ - divider * eIfftZ[i]) * ePmlZ;
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
  const float* eIfftZ = ifftZ.getData();

  #pragma omp for schedule (static)
  for (size_t z = 0; z < mDimensionSizes.nz; z++)
  {
    const float ePmlZ = pmlZ[z];
    const float eDzudznSgz = dzudznSgz[z];
    for (size_t y = 0; y < mDimensionSizes.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < mDimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, mDimensionSizes);

        mData[i] = (mData[i] * ePmlZ - (divider * eDzudznSgz * eIfftZ[i])) * ePmlZ;
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

#pragma omp declare simd
inline size_t VelocityMatrix::get1DIndex(const size_t          z,
                                         const size_t          y,
                                         const size_t          x,
                                         const DimensionSizes& dimensionSizes)
{
  return (z * dimensionSizes.ny + y) * dimensionSizes.nx + x;
}// end of get1DIndex
//----------------------------------------------------------------------------------------------------------------------
