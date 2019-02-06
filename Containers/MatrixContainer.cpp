/**
 * @file      MatrixContainer.cpp
 *
 * @author    Jiri Jaros\n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            jarosjir@fit.vutbr.cz
 *
 * @brief     The implementation file containing the matrix container.
 *
 * @version   kspaceFirstOrder3D 2.17
 *
 * @date      12 July      2012, 10:27 (created) \n
 *            06 February  2019, 16:09 (revised)
 *
 * @copyright Copyright (C) 2019 Jiri Jaros and Bradley Treeby.
 *
 * This file is part of the C++ extension of the [k-Wave Toolbox](http://www.k-wave.org).
 *
 * This file is part of the k-Wave. k-Wave is free software: you can redistribute it and/or modify it under the terms
 * of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with k-Wave.
 * If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).
 */

#include <stdexcept>

#include <Containers/MatrixContainer.h>
#include <Parameters/Parameters.h>
#include <Logger/Logger.h>

#include <MatrixClasses/RealMatrix.h>
#include <MatrixClasses/ComplexMatrix.h>
#include <MatrixClasses/FftwComplexMatrix.h>
#include <MatrixClasses/IndexMatrix.h>



using std::string;

//--------------------------------------------------------------------------------------------------------------------//
//-------------------------------------------------- Public methods --------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Constructor.
 */
MatrixContainer::MatrixContainer()
  : mContainer()
{

}// end of Constructor.
//----------------------------------------------------------------------------------------------------------------------


/**
 * Destructor.
 * No need for virtual destructor (no polymorphism).
 */
MatrixContainer::~MatrixContainer()
{
  mContainer.clear();
}// end of ~MatrixContainer
//----------------------------------------------------------------------------------------------------------------------

/**
 * This function creates the list of matrices being used in the simulation. It is done based on the
 * simulation parameters and the dimensionality. All matrices records are created here.
 */
void MatrixContainer::init()
{
  using MT = MatrixRecord::MatrixType;
  using MI = MatrixContainer::MatrixIdx;

  const Parameters& params = Parameters::getInstance();

  DimensionSizes fullDims = params.getFullDimensionSizes();
  DimensionSizes reducedDims = params.getReducedDimensionSizes();

  const bool is3DSimulation = params.isSimulation3D();

  constexpr bool kLoad         = true;
  constexpr bool kNoLoad       = false;
  constexpr bool kCheckpoint   = true;
  constexpr bool kNoCheckpoint = false;

  //--------------------------------------------- Allocate all matrices ----------------------------------------------//

  mContainer[MI::kKappa].set(MT::kReal, reducedDims, kNoLoad, kNoCheckpoint, kKappaRName);

  if (!params.getC0ScalarFlag())
  {
    mContainer[MI::kC2]   .set(MT::kReal, fullDims   ,   kLoad, kNoCheckpoint, kC0Name);
  }

  mContainer[MI::kP]      .set(MT::kReal, fullDims   , kNoLoad,   kCheckpoint, kPName);

  mContainer[MI::kRhoX]   .set(MT::kReal, fullDims   , kNoLoad,   kCheckpoint, kRhoXName);
  mContainer[MI::kRhoY]   .set(MT::kReal, fullDims   , kNoLoad,   kCheckpoint, kRhoYName);
  if (is3DSimulation)
  {
    mContainer[MI::kRhoZ] .set(MT::kReal, fullDims   , kNoLoad,   kCheckpoint, kRhoZName);
  }

  mContainer[MI::kUxSgx]  .set(MT::kReal, fullDims   , kNoLoad,   kCheckpoint, kUxSgxName);
  mContainer[MI::kUySgy]  .set(MT::kReal, fullDims   , kNoLoad,   kCheckpoint, kUySgyName);
  if (is3DSimulation)
  {
    mContainer[MI::kUzSgz].set(MT::kReal, fullDims   , kNoLoad,   kCheckpoint, kUzSgzName);
  }

  mContainer[MI::kDuxdx]  .set(MT::kReal, fullDims   , kNoLoad, kNoCheckpoint, kDuxdxName);
  mContainer[MI::kDuydy]  .set(MT::kReal, fullDims   , kNoLoad, kNoCheckpoint, kDuydyName);
  if (is3DSimulation)
  {
    mContainer[MI::kDuzdz].set(MT::kReal, fullDims   , kNoLoad, kNoCheckpoint, kDuzdzName);
  }

  if (!params.getRho0ScalarFlag())
  {
    mContainer[MI::kRho0]       .set(MT::kReal, fullDims,  kLoad, kNoCheckpoint, kRho0Name);
    mContainer[MI::kDtRho0Sgx]  .set(MT::kReal, fullDims,  kLoad, kNoCheckpoint, kRho0SgxName);
    mContainer[MI::kDtRho0Sgy]  .set(MT::kReal, fullDims,  kLoad, kNoCheckpoint, kRho0SgyName);
    if (is3DSimulation)
    {
      mContainer[MI::kDtRho0Sgz].set(MT::kReal, fullDims,  kLoad, kNoCheckpoint, kRho0SgzName);
    }
  }


  mContainer[MI::kDdxKShiftPosR] .set(MT::kComplex, DimensionSizes(reducedDims.nx, 1, 1),
                                      kLoad, kNoCheckpoint, kDdxKShiftPosRName);
  mContainer[MI::kDdyKShiftPos]  .set(MT::kComplex, DimensionSizes(1, reducedDims.ny, 1),
                                      kLoad, kNoCheckpoint, kDdyKShiftPosName);
  if (is3DSimulation)
  {
    mContainer[MI::kDdzKShiftPos].set(MT::kComplex, DimensionSizes(1, 1, reducedDims.nz),
                                      kLoad, kNoCheckpoint, kDdzKShiftPosName);
  }

  mContainer[MI::kDdxKShiftNegR] .set(MT::kComplex, DimensionSizes(reducedDims.nx ,1, 1),
                                     kLoad, kNoCheckpoint, kDdxKShiftNegRName);
  mContainer[MI::kDdyKShiftNeg]  .set(MT::kComplex, DimensionSizes(1, reducedDims.ny, 1),
                                     kLoad, kNoCheckpoint, kDdyKShiftNegName);
  if (is3DSimulation)
  {
    mContainer[MI::kDdzKShiftNeg].set(MT::kComplex, DimensionSizes(1, 1, reducedDims.nz),
                                      kLoad, kNoCheckpoint, kDdzKShiftNegName);
  }


  mContainer[MI::kPmlXSgx].set(MT::kReal, DimensionSizes(fullDims.nx, 1, 1),    kLoad, kNoCheckpoint, kPmlXSgxName);
  mContainer[MI::kPmlYSgy].set(MT::kReal, DimensionSizes(1, fullDims.ny, 1),    kLoad, kNoCheckpoint, kPmlYSgyName);
  if (is3DSimulation)
  {
    mContainer[MI::kPmlZSgz].set(MT::kReal, DimensionSizes(1, 1, fullDims.nz),    kLoad, kNoCheckpoint, kPmlZSgzName);
  }

  mContainer[MI::kPmlX]   .set(MT::kReal, DimensionSizes(fullDims.nx, 1, 1),    kLoad, kNoCheckpoint, kPmlXName);
  mContainer[MI::kPmlY]   .set(MT::kReal, DimensionSizes(1, fullDims.ny, 1),    kLoad, kNoCheckpoint, kPmlYName);
  if (is3DSimulation)
  {
    mContainer[MI::kPmlZ] .set(MT::kReal, DimensionSizes(1, 1, fullDims.nz),    kLoad, kNoCheckpoint, kPmlZName);
  }

  if (params.getNonLinearFlag())
  {
    if (! params.getBOnAScalarFlag())
    {
      mContainer[MI::kBOnA].set(MT::kReal, fullDims, kLoad, kNoCheckpoint, kBonAName);
    }
  }

  if (params.getAbsorbingFlag() != 0)
  {
    if (!((params.getC0ScalarFlag()) && (params.getAlphaCoeffScalarFlag())))
    {
      mContainer[MI::kAbsorbTau] .set(MT::kReal, fullDims   , kNoLoad, kNoCheckpoint, kAbsorbTauName);
      mContainer[MI::kAbsorbEta] .set(MT::kReal, fullDims   , kNoLoad, kNoCheckpoint, kAbsorbEtaName);
    }
    mContainer[MI::kAbsorbNabla1].set(MT::kReal, reducedDims, kNoLoad, kNoCheckpoint, kAbsorbNabla1RName);
    mContainer[MI::kAbsorbNabla2].set(MT::kReal, reducedDims, kNoLoad, kNoCheckpoint, kAbsorbNabla2RName);
  }

  // linear sensor mask
  if (params.getSensorMaskType() == Parameters::SensorMaskType::kIndex)
  {
    mContainer[MI::kSensorMaskIndex].set(MT::kIndex,
                                         DimensionSizes(params.getSensorMaskIndexSize(), 1, 1),
                                         kLoad, kNoCheckpoint, kSensorMaskIndexName);
  }

  // cuboid sensor mask
  if (params.getSensorMaskType() == Parameters::SensorMaskType::kCorners)
  {
    mContainer[MI::kSensorMaskCorners].set(MT::kIndex,
                                           DimensionSizes(6 ,params.getSensorMaskCornersSize(), 1),
                                           kLoad, kNoCheckpoint, kSensorMaskCornersName);
  }


  // if p0 source flag
  if (params.getInitialPressureSourceFlag() == 1)
  {
    mContainer[MI::kInitialPressureSourceInput].set(MT::kReal,fullDims, kLoad, kNoCheckpoint,
                                                    kInitialPressureSourceInputName);
  }


  // us_index
  if ((params.getTransducerSourceFlag() != 0) ||
      (params.getVelocityXSourceFlag() != 0)  ||
      (params.getVelocityYSourceFlag() != 0)  ||
      (params.getVelocityZSourceFlag() != 0))
  {
    mContainer[MI::kVelocitySourceIndex].set(MT::kIndex,
                                            DimensionSizes(1 ,1, params.getVelocitySourceIndexSize()),
                                            kLoad, kNoCheckpoint, kVelocitySourceIndexName);
  }

  //transducer source flag defined
  if (params.getTransducerSourceFlag() != 0)
  {
    mContainer[MI::kDelayMask]            .set(MT::kIndex,DimensionSizes(1 ,1, params.getVelocitySourceIndexSize())         ,
                                               kLoad, kNoCheckpoint, kDelayMaskName);
    mContainer[MI::kTransducerSourceInput].set(MT::kReal ,DimensionSizes(1 ,1, params.getTransducerSourceInputSize()),
                                               kLoad, kNoCheckpoint, kTransducerSourceInputName);
  }

  // p variables
  if (params.getPressureSourceFlag() != 0)
  {
    if (params.getPressureSourceMany() == 0)
    { // 1D case
      mContainer[MI::kPressureSourceInput].set(MT::kReal,
                                               DimensionSizes(1 ,1, params.getPressureSourceFlag()),
                                               kLoad, kNoCheckpoint, kPressureSourceInputName);
    }
    else
    { // 2D case
      mContainer[MI::kPressureSourceInput].set(MT::kReal,
                                               DimensionSizes(1 ,params.getPressureSourceIndexSize(),
                                               params.getPressureSourceFlag()),
                                               kLoad, kNoCheckpoint, kPressureSourceInputName);
    }

    mContainer[MI::kPressureSourceIndex].set(MT::kIndex,
                                            DimensionSizes(1 ,1, params.getPressureSourceIndexSize()),
                                            kLoad, kNoCheckpoint, kPressureSourceIndexName);
  }



  //-------------------------------------------- Velocity source flags -----------------------------------------------//
  if (params.getVelocityXSourceFlag() != 0)
  {
    if (params.getVelocitySourceMany() == 0)
    { // 1D
      mContainer[MI::kVelocityXSourceInput].set(MT::kReal,
                                                DimensionSizes(1 ,1, params.getVelocityXSourceFlag()),
                                                kLoad, kNoCheckpoint, kVelocityXSourceInputName);
    }
    else
    { // 2D
      mContainer[MI::kVelocityXSourceInput].set(MT::kReal,
                                                DimensionSizes(1 ,params.getVelocitySourceIndexSize(),
                                                params.getVelocityXSourceFlag()),
                                                kLoad, kNoCheckpoint, kVelocityXSourceInputName);
    }
  }// ux_source_input


  if (params.getVelocityYSourceFlag() != 0)
  {
    if (params.getVelocitySourceMany() == 0)
    { // 1D
      mContainer[MI::kVelocityYSourceInput].set(MT::kReal,
                                                DimensionSizes(1 ,1, params.getVelocityYSourceFlag()),
                                                kLoad, kNoCheckpoint, kVelocityYSourceInputName);
    }
    else
    { // 2D
      mContainer[MI::kVelocityYSourceInput].set(MT::kReal,
                                                DimensionSizes(1 ,params.getVelocitySourceIndexSize(),
                                                params.getVelocityYSourceFlag()),
                                                kLoad, kNoCheckpoint, kVelocityYSourceInputName);
    }
  }// uy_source_input

  if (is3DSimulation)
  {
    if (params.getVelocityZSourceFlag() != 0)
    {
      if (params.getVelocitySourceMany() == 0)
      { // 1D
        mContainer[MI::kVelocityZSourceInput].set(MT::kReal,
                                                  DimensionSizes(1 ,1, params.getVelocityZSourceFlag()),
                                                  kLoad, kNoCheckpoint, kVelocityZSourceInputName);
      }
      else
      { // 2D
        mContainer[MI::kVelocityZSourceInput].set(MT::kReal,
                                                  DimensionSizes(1 ,params.getVelocitySourceIndexSize(),
                                                  params.getVelocityZSourceFlag()),
                                                  kLoad, kNoCheckpoint, kVelocityZSourceInputName);
      }
    }// uz_source_input
  }

  /// Add sourceKappa
  if (((params.getVelocitySourceMode() == Parameters::SourceMode::kAdditive) ||
       (params.getPressureSourceMode() == Parameters::SourceMode::kAdditive)) &&
      (params.getPressureSourceFlag()  ||
       params.getVelocityXSourceFlag() || params.getVelocityYSourceFlag() || params.getVelocityZSourceFlag()))
  {
    mContainer[MI::kSourceKappa].set(MT::kReal, reducedDims, kNoLoad, kNoCheckpoint, kSourceKappaRName);
  }

  //------------------------------------------------ Nonlinear grid --------------------------------------------------//
  if (params.getNonUniformGridFlag()!= 0)
  {
    mContainer[MI::kDxudxn]   .set(MT::kReal, DimensionSizes(fullDims.nx, 1, 1), kLoad, kNoCheckpoint, kDxudxnName);
    mContainer[MI::kDyudyn]   .set(MT::kReal, DimensionSizes(1, fullDims.ny, 1), kLoad, kNoCheckpoint, kDyudynName);
    if (is3DSimulation)
    {
      mContainer[MI::kDzudzn] .set(MT::kReal, DimensionSizes(1 ,1, fullDims.nz), kLoad, kNoCheckpoint, kDzudznName);
    }

    mContainer[MI::kDxudxnSgx].set(MT::kReal, DimensionSizes(fullDims.nx, 1, 1), kLoad, kNoCheckpoint, kDxudxnSgxName);
    mContainer[MI::kDyudynSgy].set(MT::kReal, DimensionSizes(1, fullDims.ny, 1), kLoad, kNoCheckpoint, kDyudynSgyName);
    if (is3DSimulation)
    {
      mContainer[MI::kDzudznSgz].set(MT::kReal, DimensionSizes(1 ,1, fullDims.nz),
                                     kLoad, kNoCheckpoint, kDzudznSgzName);
    }
  }

  //-------------------------------------------- Non staggered velocity ----------------------------------------------//
  if (params.getStoreVelocityNonStaggeredRawFlag())
  {
    DimensionSizes shiftDims = fullDims;

    const size_t nxR = fullDims.nx / 2 + 1;
    const size_t nyR = fullDims.ny / 2 + 1;
    const size_t nzR = (is3DSimulation) ? fullDims.nz / 2 + 1 : 1;

    const size_t xCutSize = nxR         * fullDims.ny * fullDims.nz;
    const size_t yCutSize = fullDims.nx * nyR         * fullDims.nz;
    const size_t zCutSize = fullDims.nx * fullDims.ny * nzR;

    if ((xCutSize >= yCutSize) && (xCutSize >= zCutSize))
    { // X cut is the biggest
      shiftDims.nx = nxR;
    }
    else if ((yCutSize >= xCutSize) && (yCutSize >= zCutSize))
    { // Y cut is the biggest
      shiftDims.ny = nyR;
    }
    else if ((zCutSize >= xCutSize) && (zCutSize >= yCutSize))
    { // Z cut is the biggest
      shiftDims.nz = nzR;
    }
    else
    { //all are the same
      shiftDims.nx = nxR;
    }

    mContainer[MI::kTempFftwShift].set(MT::kFftw, shiftDims, kNoLoad, kNoCheckpoint, kFftwShiftTempName);

    // these three are necessary only for u_non_staggered calculation now
    mContainer[MI::kUxShifted].set(MT::kReal, fullDims, kNoLoad, kNoCheckpoint, kUxShiftedName);
    mContainer[MI::kUyShifted].set(MT::kReal, fullDims, kNoLoad, kNoCheckpoint, kUyShiftedName);
    if (is3DSimulation)
    {
      mContainer[MI::kUzShifted].set(MT::kReal, fullDims, kNoLoad, kNoCheckpoint, kUzShiftedName);
    }

    // shifts from the input file
    mContainer[MI::kXShiftNegR].set(MT::kComplex, DimensionSizes(nxR, 1  , 1  ), kLoad, kNoCheckpoint, kXShiftNegRName);
    mContainer[MI::kYShiftNegR].set(MT::kComplex, DimensionSizes(1  , nyR, 1  ), kLoad, kNoCheckpoint, kYShiftNegRName);
    if (is3DSimulation)
    {
      mContainer[MI::kZShiftNegR].set(MT::kComplex, DimensionSizes(1, 1  , nzR), kLoad, kNoCheckpoint, kZShiftNegRName);
    }
  }// u_non_staggered


  //----------------------------------------------- Temporary matrices -----------------------------------------------//
  // this matrix used to load alphaCoeff for absorbTau pre-calculation
  if ((params.getAbsorbingFlag() != 0) && (!params.getAlphaCoeffScalarFlag()))
  {
    mContainer[MI::kTemp1RealND].set(MT::kReal, fullDims ,   kLoad, kNoCheckpoint, kAlphaCoeffName);
  }
  else
  {
    mContainer[MI::kTemp1RealND].set(MT::kReal, fullDims , kNoLoad, kNoCheckpoint, kTemp1RealNDName);
  }

  mContainer[MI::kTemp2RealND]  .set(MT::kReal, fullDims   , kNoLoad, kNoCheckpoint, kTemp2RealNDName);
  if (is3DSimulation)
  {
    mContainer[MI::kTemp3RealND].set(MT::kReal, fullDims   , kNoLoad, kNoCheckpoint, kTemp3RealNDName);
  }

  mContainer[MI::kTempFftwX]    .set(MT::kFftw, reducedDims, kNoLoad, kNoCheckpoint, kFftwXTempName);
  mContainer[MI::kTempFftwY]    .set(MT::kFftw, reducedDims, kNoLoad, kNoCheckpoint, kFftwYTempName);
  if (is3DSimulation)
  {
    mContainer[MI::kTempFftwZ]  .set(MT::kFftw, reducedDims, kNoLoad, kNoCheckpoint, kFftwZTempName);
  }
}// end of init
//----------------------------------------------------------------------------------------------------------------------

/**
 * Create all matrix objects in the container.
 */
void MatrixContainer::createMatrices()
{
  using MatrixType = MatrixRecord::MatrixType;

  for (auto& it : mContainer)
  {
    if (it.second.matrixPtr != nullptr)
    {
      throw std::invalid_argument(Logger::formatMessage(kErrFmtRelocationError, it.second.matrixName.c_str()));
    }

    switch (it.second.matrixType)
    {
      case MatrixType::kReal:
      {
        it.second.matrixPtr = new RealMatrix(it.second.dimensionSizes);
        break;
      }

      case MatrixType::kComplex:
      {
        it.second.matrixPtr = new ComplexMatrix(it.second.dimensionSizes);
        break;
      }

      case MatrixType::kIndex:
      {
        it.second.matrixPtr = new IndexMatrix(it.second.dimensionSizes);
        break;
      }

      case MatrixType::kFftw:
      {
        it.second.matrixPtr = new FftwComplexMatrix(it.second.dimensionSizes);
        break;
      }

      default:
      {
        throw std::invalid_argument(Logger::formatMessage(kErrFmtBadMatrixType, it.second.matrixName.c_str()));
        break;
      }
    }// switch
  }// end for
}// end of createMatrices
//----------------------------------------------------------------------------------------------------------------------

/**
 * Free all matrix objects.
 */
void MatrixContainer::freeMatrices()
{
  for (auto& it : mContainer)
  {
    if (it.second.matrixPtr)
    {
      delete it.second.matrixPtr;
      it.second.matrixPtr = nullptr;
    }
  }
}// end of freeMatrices
//----------------------------------------------------------------------------------------------------------------------

/**
 * Load all marked matrices from the input HDF5 file.
 */
void MatrixContainer::loadDataFromInputFile()
{
  Hdf5File& inputFile = Parameters::getInstance().getInputFile();

  for (const auto& it : mContainer)
  {
    if (it.second.loadData)
    {
      it.second.matrixPtr->readData(inputFile, it.second.matrixName);
    }
  }
}// end of loadDataFromInputFile
//----------------------------------------------------------------------------------------------------------------------

/**
 * Load selected matrices from the checkpoint HDF5 file.
 */
void MatrixContainer::loadDataFromCheckpointFile()
{
  Hdf5File& checkpointFile = Parameters::getInstance().getCheckpointFile();

  for (const auto& it : mContainer)
  {
    if (it.second.checkpoint)
    {
      it.second.matrixPtr->readData(checkpointFile,it.second.matrixName);
    }
  }
}// end of loadDataFromCheckpointFile
//----------------------------------------------------------------------------------------------------------------------

/**
 * Store selected matrices into the checkpoint file.
 */
void MatrixContainer::storeDataIntoCheckpointFile()
{
  Hdf5File& checkpointFile = Parameters::getInstance().getCheckpointFile();
  auto compressionLevel    = Parameters::getInstance().getCompressionLevel();

  for (const auto& it : mContainer)
  {
    if (it.second.checkpoint)
    {
      it.second.matrixPtr->writeData(checkpointFile, it.second.matrixName, compressionLevel);
    }
  }
}// end of storeDataIntoCheckpointFile
//----------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------ Protected methods -------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Private methods --------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

