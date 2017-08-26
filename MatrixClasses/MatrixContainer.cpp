/**
 * @file        MatrixContainer.cpp
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The implementation file containing the matrix container.
 *
 * @version     kspaceFirstOrder3D 2.16
 * @date        12 July      2012, 10:27 (created) \n
 *              26 August    2017, 17:14 (revised)
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

#include <MatrixClasses/MatrixContainer.h>

#include <Parameters/Parameters.h>
#include <Utils/ErrorMessages.h>

//----------------------------------------------------------------------------//
//----------------------------- CONSTANTS ------------------------------------//
//----------------------------------------------------------------------------//




//============================================================================//
//                              TMatrixRecord                                 //
//============================================================================//

//----------------------------------------------------------------------------//
//--------------------------- Public methods ---------------------------------//
//----------------------------------------------------------------------------//


/**
 * Copy constructor of TMatrixRecord.
 * @param [in] src
 */
TMatrixRecord::TMatrixRecord(const TMatrixRecord& src) :
        MatrixPtr(src.MatrixPtr),
        MatrixDataType(src.MatrixDataType),
        dimensionSizes(src.dimensionSizes),
        LoadData(src.LoadData),
        Checkpoint(src.Checkpoint),
        HDF5MatrixName(src.HDF5MatrixName)
{

}// end of TMatrixRecord
//------------------------------------------------------------------------------


/**
 * operator = of TMatrixRecord
 * @param  [in] src
 * @return this
 */
TMatrixRecord& TMatrixRecord::operator = (const TMatrixRecord& src)
{
  if (this != &src)
  {
    MatrixPtr       = src.MatrixPtr;
    MatrixDataType  = src.MatrixDataType;
    dimensionSizes  = src.dimensionSizes;
    LoadData        = src.LoadData;
    Checkpoint      = src.Checkpoint;
    HDF5MatrixName  = src.HDF5MatrixName;
  }

  return *this;
}// end of operator =
//------------------------------------------------------------------------------



/**
 * Set all values for the record
 * @param [in] MatrixPtr        - Pointer to the MatrixClass object
 * @param [in] MatrixDataType   - Matrix data type
 * @param [in] DimensionSizes   - Dimension sizes
 * @param [in] LoadData         - Load data from file?
 * @param [in] Checkpoint       - Checkpoint this matrix?
 * @param [in] HDF5MatrixName   - HDF5 matrix name
 */
void TMatrixRecord::SetAllValues(BaseMatrix *          MatrixPtr,
                                 const TMatrixDataType  MatrixDataType,
                                 const DimensionSizes   dimensionSizes,
                                 const bool             LoadData,
                                 const bool             Checkpoint,
                                 const string           HDF5MatrixName)
{
  this->MatrixPtr       = MatrixPtr;
  this->MatrixDataType  = MatrixDataType;
  this->dimensionSizes  = dimensionSizes;
  this->LoadData        = LoadData;
  this->Checkpoint      = Checkpoint;
  this->HDF5MatrixName  = HDF5MatrixName;
}// end of SetAllValues
//------------------------------------------------------------------------------


//----------------------------------------------------------------------------//
//------------------------- Protected methods --------------------------------//
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//-------------------------- Private methods ---------------------------------//
//----------------------------------------------------------------------------//



//============================================================================//
//                              TMatrixContainer                              //
//============================================================================//

//----------------------------------------------------------------------------//
//--------------------------- Public methods ---------------------------------//
//----------------------------------------------------------------------------//


/**
 * Destructor of TMatrixContainer.
 */
TMatrixContainer::~TMatrixContainer()
{
  MatrixContainer.clear();
}// end of ~TMatrixContainer
//------------------------------------------------------------------------------


/**
 * Create all matrix objects in the container.
 * @throw errors cause an exception bad_alloc.
 */
void TMatrixContainer::CreateAllObjects()
{
  for (TMatrixRecordContainer::iterator it = MatrixContainer.begin(); it != MatrixContainer.end(); it++)
  {
    if (it->second.MatrixPtr != NULL)
    {
      PrintErrorAndThrowException(kErrFmtRelocationError,
                                  it->second.HDF5MatrixName,
                                 __FILE__,__LINE__);
    }

    switch (it->second.MatrixDataType)
    {
      case TMatrixRecord::mdtReal:
      {
        it->second.MatrixPtr = new RealMatrix(it->second.dimensionSizes);
        break;
      }

      case TMatrixRecord::mdtComplex:
      {
        it->second.MatrixPtr = new ComplexMatrix(it->second.dimensionSizes);
        break;
      }

      case TMatrixRecord::mdtIndex:
      {
        it->second.MatrixPtr = new IndexMatrix(it->second.dimensionSizes);
        break;
      }

      case TMatrixRecord::mdtFFTW:
      {
        it->second.MatrixPtr = new FftwComplexMatrix(it->second.dimensionSizes);
        break;
      }

      case TMatrixRecord::mdtUxyz:
      {
        it->second.MatrixPtr = new VelocityMatrix(it->second.dimensionSizes);
        break;
      }

      default:
      {
        PrintErrorAndThrowException(kErrFmtBadMatrixType,
                                    it->second.HDF5MatrixName,
                                    __FILE__,__LINE__);
        break;
      }
    }// switch
  }// end for
}// end of CreateAllObjects
//-----------------------------------------------------------------------------


/**
 * Load all marked matrices from the HDF5 file.
 * @param [in] HDF5_File - HDF5 file handle
 */
void TMatrixContainer::LoadDataFromInputHDF5File(Hdf5File & HDF5_File)
{
  for (TMatrixRecordContainer::iterator it = MatrixContainer.begin(); it != MatrixContainer.end(); it++)
  {
    if (it->second.LoadData)
    {
      it->second.MatrixPtr->readData(HDF5_File, it->second.HDF5MatrixName.c_str());
    }
  }
}// end of LoadMatricesDataFromDisk
//------------------------------------------------------------------------------


/**
 * Load selected matrices from checkpoint HDF5 file.
 * @param [in] HDF5_File - HDF5 file handle
 */
void TMatrixContainer::LoadDataFromCheckpointHDF5File(Hdf5File & HDF5_File)
{
  for (TMatrixRecordContainer::iterator it = MatrixContainer.begin(); it != MatrixContainer.end(); it++)
  {
    if (it->second.Checkpoint)
    {
      it->second.MatrixPtr->readData(HDF5_File, it->second.HDF5MatrixName.c_str());
    }
  }
}// end of LoadDataFromCheckpointHDF5File
//------------------------------------------------------------------------------

/**
 * Store selected matrices into the checkpoint file.
 * @param [in] HDF5_File
 */
void TMatrixContainer::StoreDataIntoCheckpointHDF5File(Hdf5File & HDF5_File)
{
  for (TMatrixRecordContainer::iterator it = MatrixContainer.begin(); it != MatrixContainer.end(); it++)
  {
    if (it->second.Checkpoint)
    {
      it->second.MatrixPtr->writeData(HDF5_File,
                                                it->second.HDF5MatrixName.c_str(),
                                                Parameters::getInstance().getCompressionLevel());
    }
  }
}// end of StoreDataIntoCheckpointHDF5File
//------------------------------------------------------------------------------


/**
 * Free all matrix objects.
 *
 */
void TMatrixContainer::FreeAllMatrices()
{
  for (TMatrixRecordContainer::iterator it = MatrixContainer.begin(); it != MatrixContainer.end(); it++)
  {
    if (it->second.MatrixPtr)
    {
      delete it->second.MatrixPtr;
      it->second.MatrixPtr = NULL;
    }
  }
}// end of FreeAllMatrices
//------------------------------------------------------------------------------




/**
 * This function defines common matrices in K-Wave.
 * All matrices records are created here.
 */
void TMatrixContainer::AddMatricesIntoContainer()
{

  Parameters& Params = Parameters::getInstance();

  DimensionSizes FullDims = Params.getFullDimensionSizes();
  DimensionSizes ReducedDims = Params.getReducedDimensionSizes();

  const bool LOAD         = true;
  const bool NOLOAD       = false;
  const bool CHECKPOINT   = true;
  const bool NOCHECKPOINT = false;


  //----------------------Allocate all matrices ----------------------------//

  MatrixContainer[kappa] .SetAllValues(NULL ,TMatrixRecord::mdtReal, ReducedDims, NOLOAD, NOCHECKPOINT, kKappaRName);
  if (!Params.getC0ScalarFlag())
  {
    MatrixContainer[c2]  .SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   ,   LOAD, NOCHECKPOINT, kC0Name);
  }
  MatrixContainer[p]     .SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   , NOLOAD,   CHECKPOINT, kPName);

  MatrixContainer[rhox]  .SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   , NOLOAD,   CHECKPOINT, kRhoXName);
  MatrixContainer[rhoy]  .SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   , NOLOAD,   CHECKPOINT, kRhoYName);
  MatrixContainer[rhoz]  .SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   , NOLOAD,   CHECKPOINT, kRhoZName);

  MatrixContainer[ux_sgx].SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   , NOLOAD,   CHECKPOINT, kUxSgxName);
  MatrixContainer[uy_sgy].SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   , NOLOAD,   CHECKPOINT, kUySgyName);
  MatrixContainer[uz_sgz].SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   , NOLOAD,   CHECKPOINT, kUzSgzName);

  MatrixContainer[duxdx] .SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   , NOLOAD, NOCHECKPOINT, kDuxdxName);
  MatrixContainer[duydy] .SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   , NOLOAD, NOCHECKPOINT, kDuydyName);
  MatrixContainer[duzdz] .SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   , NOLOAD, NOCHECKPOINT, kDuzdzName);

  if (!Params.getRho0ScalarFlag())
  {
    MatrixContainer[rho0]       .SetAllValues(NULL,TMatrixRecord::mdtReal, FullDims,  LOAD, NOCHECKPOINT, kRho0Name);
    MatrixContainer[dt_rho0_sgx].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDims,  LOAD, NOCHECKPOINT, kRho0SgxName);
    MatrixContainer[dt_rho0_sgy].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDims,  LOAD, NOCHECKPOINT, kRho0SgyName);
    MatrixContainer[dt_rho0_sgz].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDims,  LOAD, NOCHECKPOINT, kRho0SgzName);
  }


  MatrixContainer[ddx_k_shift_pos].SetAllValues(NULL,TMatrixRecord::mdtComplex, DimensionSizes(ReducedDims.nx, 1, 1), LOAD, NOCHECKPOINT, kDdxKShiftPosRName);
  MatrixContainer[ddy_k_shift_pos].SetAllValues(NULL,TMatrixRecord::mdtComplex, DimensionSizes(1, ReducedDims.ny, 1), LOAD, NOCHECKPOINT, kDdyKShiftPosName);
  MatrixContainer[ddz_k_shift_pos].SetAllValues(NULL,TMatrixRecord::mdtComplex, DimensionSizes(1, 1, ReducedDims.nz), LOAD, NOCHECKPOINT, kDdzKShiftPosName);

  MatrixContainer[ddx_k_shift_neg].SetAllValues(NULL,TMatrixRecord::mdtComplex, DimensionSizes(ReducedDims.nx ,1, 1), LOAD, NOCHECKPOINT, kDdxKShiftNegRName);
  MatrixContainer[ddy_k_shift_neg].SetAllValues(NULL,TMatrixRecord::mdtComplex, DimensionSizes(1, ReducedDims.ny, 1), LOAD, NOCHECKPOINT, kDdyKShiftNegName);
  MatrixContainer[ddz_k_shift_neg].SetAllValues(NULL,TMatrixRecord::mdtComplex, DimensionSizes(1, 1, ReducedDims.nz), LOAD, NOCHECKPOINT, kDdzKShiftNegName);


  MatrixContainer[pml_x_sgx] .SetAllValues(NULL,TMatrixRecord::mdtReal,         DimensionSizes(FullDims.nx, 1, 1),    LOAD, NOCHECKPOINT, kPmlXSgxName);
  MatrixContainer[pml_y_sgy] .SetAllValues(NULL,TMatrixRecord::mdtReal,         DimensionSizes(1, FullDims.ny, 1),    LOAD, NOCHECKPOINT, kPmlYSgyName);
  MatrixContainer[pml_z_sgz] .SetAllValues(NULL,TMatrixRecord::mdtReal,         DimensionSizes(1, 1, FullDims.nz),    LOAD, NOCHECKPOINT, kPmlZSgzName);

  MatrixContainer[pml_x]     .SetAllValues(NULL,TMatrixRecord::mdtReal,         DimensionSizes(FullDims.nx, 1, 1),    LOAD, NOCHECKPOINT, kPmlXName);
  MatrixContainer[pml_y]     .SetAllValues(NULL,TMatrixRecord::mdtReal,         DimensionSizes(1, FullDims.ny, 1),    LOAD, NOCHECKPOINT, kPmlYName);
  MatrixContainer[pml_z]     .SetAllValues(NULL,TMatrixRecord::mdtReal,         DimensionSizes(1, 1, FullDims.nz),    LOAD, NOCHECKPOINT, kPmlZName);

  if (Params.getNonLinearFlag())
  {
    if (! Params.getBOnAScalarFlag())
    {
      MatrixContainer[BonA]      .SetAllValues   (NULL,TMatrixRecord::mdtReal, FullDims ,   LOAD, NOCHECKPOINT, kBonAName);
    }
  }

  if (Params.getAbsorbingFlag() != 0)
  {
    if (!((Params.getC0ScalarFlag()) && (Params.getAlphaCoeffScalarFlag())))
    {
      MatrixContainer[absorb_tau].SetAllValues (NULL,TMatrixRecord::mdtReal, FullDims   , NOLOAD, NOCHECKPOINT, kAbsorbTauName);
      MatrixContainer[absorb_eta].SetAllValues (NULL,TMatrixRecord::mdtReal, FullDims   , NOLOAD, NOCHECKPOINT, kAbsorbEtaName);
    }
    MatrixContainer[absorb_nabla1].SetAllValues(NULL,TMatrixRecord::mdtReal, ReducedDims, NOLOAD, NOCHECKPOINT, kAbsorbNabla1RName);
    MatrixContainer[absorb_nabla2].SetAllValues(NULL,TMatrixRecord::mdtReal, ReducedDims, NOLOAD, NOCHECKPOINT, kAbsorbNabla2RName);
  }

  // linear sensor mask
  if (Params.getSensorMaskType() == Parameters::SensorMaskType::kIndex)
  {
    MatrixContainer[sensor_mask_index].SetAllValues(NULL,TMatrixRecord::mdtIndex,
                                                    DimensionSizes(Params.getSensorMaskIndexSize(), 1, 1),
                                                    LOAD, NOCHECKPOINT, kSensorMaskIndexName);
  }

  // cuboid sensor mask
  if (Params.getSensorMaskType() == Parameters::SensorMaskType::kCorners)
  {
    MatrixContainer[sensor_mask_corners].SetAllValues(NULL,TMatrixRecord::mdtIndex,
                                                      DimensionSizes(6 ,Params.getSensorMaskCornersSize(), 1),
                                                      LOAD, NOCHECKPOINT, kSensorMaskCornersName);
  }


  // if p0 source flag
  if (Params.getInitialPressureSourceFlag() == 1)
  {
    MatrixContainer[p0_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDims, LOAD, NOCHECKPOINT, kInitialPressureSourceInputName);
  }


  // us_index
  if ((Params.getTransducerSourceFlag() != 0) ||
      (Params.getVelocityXSourceFlag() != 0)         ||
      (Params.getVelocityYSourceFlag() != 0)         ||
      (Params.getVelocityZSourceFlag() != 0))
  {
    MatrixContainer[u_source_index].SetAllValues(NULL,TMatrixRecord::mdtIndex,
                                                 DimensionSizes(1 ,1, Params.getVelocitySourceIndexSize()),
                                                 LOAD, NOCHECKPOINT, kVelocitySourceIndexName);
  }

  //transducer source flag defined
  if (Params.getTransducerSourceFlag() != 0)
  {
    MatrixContainer[delay_mask]             .SetAllValues(NULL,TMatrixRecord::mdtIndex,DimensionSizes(1 ,1, Params.getVelocitySourceIndexSize())         ,
                                                          LOAD, NOCHECKPOINT, kDelayMaskName);
    MatrixContainer[transducer_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal ,DimensionSizes(1 ,1, Params.getTransducerSourceInputSize()),
                                                          LOAD, NOCHECKPOINT, kInitialPressureSourceInputName);
  }

  // p variables
  if (Params.getPressureSourceFlag() != 0)
  {
    if (Params.getPressureSourceMany() == 0)
    { // 1D case
      MatrixContainer[p_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal,
                                                   DimensionSizes(1 ,1, Params.getPressureSourceFlag()),
                                                   LOAD, NOCHECKPOINT, kPressureSourceInputName);
    }
    else
    { // 2D case
      MatrixContainer[p_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal,
                                                   DimensionSizes(1 ,Params.getPressureSourceIndexSize(),Params.getPressureSourceFlag()),
                                                   LOAD, NOCHECKPOINT, kPressureSourceInputName);
    }

    MatrixContainer[p_source_index].SetAllValues(NULL,TMatrixRecord::mdtIndex,
                                                 DimensionSizes(1 ,1, Params.getPressureSourceIndexSize()),
                                                 LOAD, NOCHECKPOINT, kPressureSourceIndexName);
  }



    //----------------------------uxyz source flags---------------------------//
  if (Params.getVelocityXSourceFlag() != 0)
  {
    if (Params.getVelocitySourceMany() == 0)
    { // 1D
      MatrixContainer[ux_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal,
                                                    DimensionSizes(1 ,1, Params.getVelocityXSourceFlag()),
                                                    LOAD, NOCHECKPOINT, kVelocityXSourceInputName);
    }
    else
    { // 2D
      MatrixContainer[ux_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal,
                                                    DimensionSizes(1 ,Params.getVelocitySourceIndexSize(),Params.getVelocityXSourceFlag()),
                                                    LOAD, NOCHECKPOINT, kVelocityXSourceInputName);
    }
  }// ux_source_input


  if (Params.getVelocityYSourceFlag() != 0)
  {
    if (Params.getVelocitySourceMany() == 0)
    { // 1D
      MatrixContainer[uy_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal,
                                                    DimensionSizes(1 ,1, Params.getVelocityYSourceFlag()),
                                                    LOAD, NOCHECKPOINT, kVelocityYSourceInputName);
    }
    else
    { // 2D
      MatrixContainer[uy_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal,
                                                    DimensionSizes(1 ,Params.getVelocitySourceIndexSize(),Params.getVelocityYSourceFlag()),
                                                    LOAD, NOCHECKPOINT, kVelocityYSourceInputName);
    }
  }// uy_source_input

  if (Params.getVelocityZSourceFlag() != 0)
  {
    if (Params.getVelocitySourceMany() == 0)
    { // 1D
      MatrixContainer[uz_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal,
                                                    DimensionSizes(1 ,1, Params.getVelocityZSourceFlag()),
                                                    LOAD, NOCHECKPOINT, kVelocityZSourceInputName);
    }
    else
    { // 2D
      MatrixContainer[uz_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal,
                                                    DimensionSizes(1 ,Params.getVelocitySourceIndexSize(),Params.getVelocityZSourceFlag()),
                                                    LOAD, NOCHECKPOINT, kVelocityZSourceInputName);
    }
  }// uz_source_input


  //-- Nonlinear grid
  if (Params.getNonUniformGridFlag()!= 0)
  {
    MatrixContainer[dxudxn]    .SetAllValues(NULL,TMatrixRecord::mdtReal, DimensionSizes(FullDims.nx, 1, 1), LOAD, NOCHECKPOINT, kDxudxnName);
    MatrixContainer[dyudyn]    .SetAllValues(NULL,TMatrixRecord::mdtReal, DimensionSizes(1, FullDims.ny, 1), LOAD, NOCHECKPOINT, kDyudynName);
    MatrixContainer[dzudzn]    .SetAllValues(NULL,TMatrixRecord::mdtReal, DimensionSizes(1 ,1, FullDims.nz), LOAD, NOCHECKPOINT, kDzudznName);

    MatrixContainer[dxudxn_sgx].SetAllValues(NULL,TMatrixRecord::mdtReal, DimensionSizes(FullDims.nx, 1, 1), LOAD, NOCHECKPOINT, kDxudxnSgxName);
    MatrixContainer[dyudyn_sgy].SetAllValues(NULL,TMatrixRecord::mdtReal, DimensionSizes(1, FullDims.ny, 1), LOAD, NOCHECKPOINT, kDyudynSgyName);
    MatrixContainer[dzudzn_sgz].SetAllValues(NULL,TMatrixRecord::mdtReal, DimensionSizes(1 ,1, FullDims.nz), LOAD, NOCHECKPOINT, kDzudznSgzName);
  }

  //-- u_non_staggered_raw
  if (Params.getStoreVelocityNonStaggeredRawFlag())
  {
    DimensionSizes ShiftDims = FullDims;

    size_t X_2 = FullDims.nx / 2 + 1;
    size_t Y_2 = FullDims.ny / 2 + 1;
    size_t Z_2 = FullDims.nz / 2 + 1;

    size_t XCutSize = X_2        * FullDims.ny * FullDims.nz;
    size_t YCutSize = FullDims.nx * Y_2        * FullDims.nz;
    size_t ZCutSize = FullDims.nx * FullDims.ny * Z_2;

    if ((XCutSize >= YCutSize) && (XCutSize >= ZCutSize))
    { // X cut is the biggest
      ShiftDims.nx = X_2;
    }
    else if ((YCutSize >= XCutSize) && (YCutSize >= ZCutSize))
    { // Y cut is the biggest
      ShiftDims.ny = Y_2;
    }
    else if ((ZCutSize >= XCutSize) && (ZCutSize >= YCutSize))
    { // Z cut is the biggest
      ShiftDims.nz = Z_2;
    }
    else
    { //all are the same
      ShiftDims.nx = X_2;
    }

    MatrixContainer[FFT_shift_temp].SetAllValues(NULL, TMatrixRecord::mdtFFTW, ShiftDims, NOLOAD, NOCHECKPOINT, kCufftShiftTempName);

    // these three are necessary only for u_non_staggered calculation now
    MatrixContainer[ux_shifted].SetAllValues(NULL, TMatrixRecord::mdtReal, FullDims, NOLOAD, NOCHECKPOINT, kUxShiftedName);
    MatrixContainer[uy_shifted].SetAllValues(NULL, TMatrixRecord::mdtReal, FullDims, NOLOAD, NOCHECKPOINT, kUyShiftedName);
    MatrixContainer[uz_shifted].SetAllValues(NULL, TMatrixRecord::mdtReal, FullDims, NOLOAD, NOCHECKPOINT, kUzShiftedName);

    // shifts from the input file
    MatrixContainer[x_shift_neg_r].SetAllValues(NULL, TMatrixRecord::mdtComplex, DimensionSizes(X_2, 1  , 1  ), LOAD, NOCHECKPOINT, kXShiftNegRName);
    MatrixContainer[y_shift_neg_r].SetAllValues(NULL, TMatrixRecord::mdtComplex, DimensionSizes(1  , Y_2, 1  ), LOAD, NOCHECKPOINT, kYShiftNegRName);
    MatrixContainer[z_shift_neg_r].SetAllValues(NULL, TMatrixRecord::mdtComplex, DimensionSizes(1  , 1  , Z_2), LOAD, NOCHECKPOINT, kZShiftNegRName);
  }// u_non_staggered


  //------------------------------------------------------------------------//
  //--------------------- Temporary matrices -------------------------------//
  //------------------------------------------------------------------------//
  // this matrix used to load alpha_coeff for absorb_tau pre-calculation

  if ((Params.getAbsorbingFlag() != 0) && (!Params.getAlphaCoeffScalarFlag()))
  {
    MatrixContainer[Temp_1_RS3D].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDims ,   LOAD, NOCHECKPOINT, kAlphaCoeffName);
  }
  else
  {
    MatrixContainer[Temp_1_RS3D].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDims , NOLOAD, NOCHECKPOINT, kTemp1Real3DName);
  }

  MatrixContainer[Temp_2_RS3D].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDims   , NOLOAD, NOCHECKPOINT, kTemp2Real3DName);
  MatrixContainer[Temp_3_RS3D].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDims   , NOLOAD, NOCHECKPOINT, kTemp3Real3DName);

  MatrixContainer[FFT_X_temp].SetAllValues(NULL, TMatrixRecord::mdtFFTW, ReducedDims, NOLOAD, NOCHECKPOINT, kCufftXTempName);
  MatrixContainer[FFT_Y_temp].SetAllValues(NULL, TMatrixRecord::mdtFFTW, ReducedDims, NOLOAD, NOCHECKPOINT, kCufftYTempName);
  MatrixContainer[FFT_Z_temp].SetAllValues(NULL, TMatrixRecord::mdtFFTW, ReducedDims, NOLOAD, NOCHECKPOINT, kCufftZTempName);
}// end of InitMatrixContainer
//------------------------------------------------------------------------------


//----------------------------------------------------------------------------//
//------------------------- Protected methods --------------------------------//
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
//-------------------------- Private methods ---------------------------------//
//----------------------------------------------------------------------------//


/**
 * Print error and and throw an exception.
 * @throw bad_alloc
 *
 * @param [in] FMT - format of error
 * @param [in] HDF5MatrixName - HDF5 dataset name
 * @param [in] File - File of error
 * @param [in] Line - Line of error
 *
 */
void TMatrixContainer::PrintErrorAndThrowException(const char * FMT,
                                                   const string HDF5MatrixName,
                                                   const char * File,
                                                   const int    Line)
{
  fprintf(stderr,FMT, HDF5MatrixName.c_str(), File, Line);
  throw bad_alloc();
}// end of PrintErrorAndAbort
//------------------------------------------------------------------------------



//============================================================================//
//                        TOutputStreamContainer                              //
//============================================================================//

//----------------------------------------------------------------------------//
//--------------------------- Public methods ---------------------------------//
//----------------------------------------------------------------------------//


/**
 * Destructor
 */
TOutputStreamContainer::~TOutputStreamContainer()
{
  OutputStreamContainer.clear();
}// end of Destructor
//------------------------------------------------------------------------------

/**
 * Add all streams in simulation in the container, set all streams records here!
 * Please note, the Matrixcontainer has to be populated before calling this routine.
 *
 * @param [in] MatrixContainer - matrix container to link the steams with
 *                               sampled matrices and sensor masks
 */
void TOutputStreamContainer::AddStreamsIntoContainer(TMatrixContainer & MatrixContainer)
{

  Parameters& Params = Parameters::getInstance();

  float * TempBufferX = MatrixContainer.GetMatrix<RealMatrix>(Temp_1_RS3D).getData();
  float * TempBufferY = MatrixContainer.GetMatrix<RealMatrix>(Temp_2_RS3D).getData();
  float * TempBufferZ = MatrixContainer.GetMatrix<RealMatrix>(Temp_3_RS3D).getData();

  //--------------------- Pressure ------------------/
  if (Params.getStorePressureRawFlag())
  {
    OutputStreamContainer[p_sensor_raw] = CreateNewOutputStream(MatrixContainer,
                                                                p,
                                                                kPressureRawName.c_str(),
                                                                TBaseOutputHDF5Stream::roNONE,
                                                                TempBufferX);
  }// IsStore_p_raw

  if (Params.getStorePressureRmsFlag())
  {
    OutputStreamContainer[p_sensor_rms] = CreateNewOutputStream(MatrixContainer,
                                                                p,
                                                                kPressureRmsName.c_str(),
                                                                TBaseOutputHDF5Stream::roRMS);
  }

  if (Params.getStorePressureMaxFlag())
  {
    OutputStreamContainer[p_sensor_max] = CreateNewOutputStream(MatrixContainer,
                                                                p,
                                                                kPressureMaxName.c_str(),
                                                                TBaseOutputHDF5Stream::roMAX);
  }

  if (Params.getStorePressureMinFlag())
  {
    OutputStreamContainer[p_sensor_min] = CreateNewOutputStream(MatrixContainer,
                                                                p,
                                                                kPressureMinName.c_str(),
                                                                TBaseOutputHDF5Stream::roMIN);
  }

  if (Params.getStorePressureMaxAllFlag())
  {
    OutputStreamContainer[p_sensor_max_all] =
            new TWholeDomainOutputHDF5Stream(Params.getOutputFile(),
                                             kPressureMaxAllName.c_str(),
                                             MatrixContainer.GetMatrix<RealMatrix>(p),
                                             TBaseOutputHDF5Stream::roMAX);
  }

  if (Params.getStorePressureMinAllFlag())
  {
    OutputStreamContainer[p_sensor_min_all] =
            new TWholeDomainOutputHDF5Stream(Params.getOutputFile(),
                                             kPressureMinAllName.c_str(),
                                             MatrixContainer.GetMatrix<RealMatrix>(p),
                                             TBaseOutputHDF5Stream::roMIN);
  }

  //--------------------- Velocity ------------------/
  if (Params.getStoreVelocityRawFlag())
  {
    OutputStreamContainer[ux_sensor_raw] = CreateNewOutputStream(MatrixContainer,
                                                                 ux_sgx,
                                                                 kUxName.c_str(),
                                                                 TBaseOutputHDF5Stream::roNONE,
                                                                 TempBufferX);
    OutputStreamContainer[uy_sensor_raw] = CreateNewOutputStream(MatrixContainer,
                                                                 uy_sgy,
                                                                 kUyName.c_str(),
                                                                 TBaseOutputHDF5Stream::roNONE,
                                                                 TempBufferY);
    OutputStreamContainer[uz_sensor_raw] = CreateNewOutputStream(MatrixContainer,
                                                                 uz_sgz,
                                                                 kUzName.c_str(),
                                                                 TBaseOutputHDF5Stream::roNONE,
                                                                 TempBufferZ);
  }

  if (Params.getStoreVelocityNonStaggeredRawFlag())
  {
    OutputStreamContainer[ux_shifted_sensor_raw] = CreateNewOutputStream(MatrixContainer,
                                                                         ux_shifted,
                                                                         kUxNonStaggeredName.c_str(),
                                                                         TBaseOutputHDF5Stream::roNONE,
                                                                         TempBufferX);
    OutputStreamContainer[uy_shifted_sensor_raw] = CreateNewOutputStream(MatrixContainer,
                                                                         uy_shifted,
                                                                         kUyNonStaggeredName.c_str(),
                                                                         TBaseOutputHDF5Stream::roNONE,
                                                                         TempBufferY);
    OutputStreamContainer[uz_shifted_sensor_raw] = CreateNewOutputStream(MatrixContainer,
                                                                         uz_shifted,
                                                                         kUzNonStaggeredName.c_str(),
                                                                         TBaseOutputHDF5Stream::roNONE,
                                                                         TempBufferZ);
  }

  if (Params.getStoreVelocityRmsFlag())
  {
    OutputStreamContainer[ux_sensor_rms] = CreateNewOutputStream(MatrixContainer,
                                                                 ux_sgx,
                                                                 kUxRmsName.c_str(),
                                                                 TBaseOutputHDF5Stream::roRMS);
    OutputStreamContainer[uy_sensor_rms] = CreateNewOutputStream(MatrixContainer,
                                                                 uy_sgy,
                                                                 kUyRmsName.c_str(),
                                                                 TBaseOutputHDF5Stream::roRMS);
    OutputStreamContainer[uz_sensor_rms] = CreateNewOutputStream(MatrixContainer,
                                                                 uz_sgz,
                                                                 kUzRmsName.c_str(),
                                                                 TBaseOutputHDF5Stream::roRMS);
  }

  if (Params.getStoreVelocityMaxFlag())
  {
    OutputStreamContainer[ux_sensor_max] = CreateNewOutputStream(MatrixContainer,
                                                                 ux_sgx,
                                                                 kUxMaxName.c_str(),
                                                                 TBaseOutputHDF5Stream::roMAX);
    OutputStreamContainer[uy_sensor_max] = CreateNewOutputStream(MatrixContainer,
                                                                 uy_sgy,
                                                                 kUyMaxName.c_str(),
                                                                 TBaseOutputHDF5Stream::roMAX);
    OutputStreamContainer[uz_sensor_max] = CreateNewOutputStream(MatrixContainer,
                                                                 uz_sgz,
                                                                 kUzMaxName.c_str(),
                                                                 TBaseOutputHDF5Stream::roMAX);
  }

  if (Params.getStoreVelocityMinFlag())
  {
    OutputStreamContainer[ux_sensor_min] = CreateNewOutputStream(MatrixContainer,
                                                                 ux_sgx,
                                                                 kUxMinName.c_str(),
                                                                 TBaseOutputHDF5Stream::roMIN);
    OutputStreamContainer[uy_sensor_min] = CreateNewOutputStream(MatrixContainer,
                                                                 uy_sgy,
                                                                 kUyMinName.c_str(),
                                                                 TBaseOutputHDF5Stream::roMIN);
    OutputStreamContainer[uz_sensor_min] = CreateNewOutputStream(MatrixContainer,
                                                                 uz_sgz,
                                                                 kUzMinName.c_str(),
                                                                 TBaseOutputHDF5Stream::roMIN);
  }

  if (Params.getStoreVelocityMaxAllFlag())
  {
    OutputStreamContainer[ux_sensor_max_all] =
            new TWholeDomainOutputHDF5Stream(Params.getOutputFile(),
                                             kUxMaxAllName.c_str(),
                                             MatrixContainer.GetMatrix<RealMatrix>(ux_sgx),
                                             TBaseOutputHDF5Stream::roMAX);
    OutputStreamContainer[uy_sensor_max_all] =
            new TWholeDomainOutputHDF5Stream(Params.getOutputFile(),
                                             kUyMaxAllName.c_str(),
                                             MatrixContainer.GetMatrix<RealMatrix>(uy_sgy),
                                             TBaseOutputHDF5Stream::roMAX);
    OutputStreamContainer[uz_sensor_max_all] =
            new TWholeDomainOutputHDF5Stream(Params.getOutputFile(),
                                             kUzMaxAllName.c_str(),
                                             MatrixContainer.GetMatrix<RealMatrix>(uz_sgz),
                                             TBaseOutputHDF5Stream::roMAX);
  }

  if (Params.getStoreVelocityMinAllFlag())
  {
    OutputStreamContainer[ux_sensor_min_all] =
            new TWholeDomainOutputHDF5Stream(Params.getOutputFile(),
                                             kUxMinAllName.c_str(),
                                             MatrixContainer.GetMatrix<RealMatrix>(ux_sgx),
                                             TBaseOutputHDF5Stream::roMIN);
    OutputStreamContainer[uy_sensor_min_all] =
            new TWholeDomainOutputHDF5Stream(Params.getOutputFile(),
                                             kUyMinAllName.c_str(),
                                             MatrixContainer.GetMatrix<RealMatrix>(uy_sgy),
                                             TBaseOutputHDF5Stream::roMIN);
    OutputStreamContainer[uz_sensor_min_all] =
            new TWholeDomainOutputHDF5Stream(Params.getOutputFile(),
                                             kUzMinAllName.c_str(),
                                             MatrixContainer.GetMatrix<RealMatrix>(uz_sgz),
                                             TBaseOutputHDF5Stream::roMIN);
  }
}// end of AddStreamsdIntoContainer
//------------------------------------------------------------------------------

/**
 * Create all streams.
 */
void TOutputStreamContainer::CreateStreams()
{
  for (TOutputStreamMap::iterator it = OutputStreamContainer.begin(); it != OutputStreamContainer.end(); it++)
  {
    if (it->second)
    {
      (it->second)->Create();
    }
  }
}// end of CreateStreams
//------------------------------------------------------------------------------

/**
 * Reopen all streams after restarting form checkpoint.
 */
void TOutputStreamContainer::ReopenStreams()
{
  for (TOutputStreamMap::iterator it = OutputStreamContainer.begin(); it != OutputStreamContainer.end(); it++)
  {
    if (it->second)
    {
      (it->second)->Reopen();
    }
  }
}// end of ReopenStreams
//------------------------------------------------------------------------------


/**
 * Sample all streams.
 */
void TOutputStreamContainer::SampleStreams()
{
  for (TOutputStreamMap::iterator it = OutputStreamContainer.begin(); it != OutputStreamContainer.end(); it++)
  {
    if (it->second)
    {
      (it->second)->Sample();
    }
  }
}// end of SampleStreams
//------------------------------------------------------------------------------


/**
 * Checkpoint streams without post-processing (flush to the file).
 */
void TOutputStreamContainer::CheckpointStreams()
{
  for (TOutputStreamMap::iterator it = OutputStreamContainer.begin(); it != OutputStreamContainer.end(); it++)
  {
    if (it->second)
    {
      (it->second)->Checkpoint();
    }
  }
}// end of CheckpointStreams
//------------------------------------------------------------------------------

/**
 * /// Post-process all streams and flush them to the file.
 */
void TOutputStreamContainer::PostProcessStreams()
{
  for (TOutputStreamMap::iterator it = OutputStreamContainer.begin(); it != OutputStreamContainer.end(); it++)
  {
    if (it->second)
    {
      (it->second)->PostProcess();
    }
  }
}// end of CheckpointStreams
//------------------------------------------------------------------------------


/**
 * Close all streams (apply post-processing if necessary, flush data and close).
 */
void TOutputStreamContainer::CloseStreams()
{
  for (TOutputStreamMap::iterator it = OutputStreamContainer.begin(); it != OutputStreamContainer.end(); it++)
  {
    if (it->second)
    {
      (it->second)->Close();
    }
  }
}// end of CloseStreams
//------------------------------------------------------------------------------

/**
 *  Free all streams- destroy them.
 */
void TOutputStreamContainer::FreeAllStreams()
{
  for (TOutputStreamMap::iterator it = OutputStreamContainer.begin(); it != OutputStreamContainer.end(); it++)
  {
    if (it->second)
    {
      delete it->second;
    }
  }
  OutputStreamContainer.clear();
}// end of FreeAllStreams
//------------------------------------------------------------------------------


//----------------------------------------------------------------------------//
//--------------------------- Protected methods ------------------------------//
//----------------------------------------------------------------------------//


/**
 * Create a new output stream.
 * @param [in] MatrixContainer  - name of the HDF5 dataset or group
 * @param [in] SampledMatrixID  - code id of the matrix
 * @param [in] HDF5_DatasetName - name of the HDF5 dataset or group
 * @param [in] ReductionOp      - reduction operator
 * @param [in] BufferToReuse   - buffer to reuse
 * @return new output stream with defined links.
 */
TBaseOutputHDF5Stream * TOutputStreamContainer::CreateNewOutputStream(TMatrixContainer & MatrixContainer,
                                                                      const TMatrixID    SampledMatrixID,
                                                                      const char *       HDF5_DatasetName,
                                                                      const TBaseOutputHDF5Stream::TReductionOperator  ReductionOp,
                                                                      float *            BufferToReuse)
{
  Parameters& Params = Parameters::getInstance();

  TBaseOutputHDF5Stream * Stream = NULL;

  if (Params.getSensorMaskType() == Parameters::SensorMaskType::kIndex)
  {
    Stream = new TIndexOutputHDF5Stream(Params.getOutputFile(),
                                        HDF5_DatasetName,
                                        MatrixContainer.GetMatrix<RealMatrix>(SampledMatrixID),
                                        MatrixContainer.GetMatrix<IndexMatrix>(sensor_mask_index),
                                        ReductionOp,
                                        BufferToReuse);
  }
  else
  {
    Stream = new TCuboidOutputHDF5Stream(Params.getOutputFile(),
                                         HDF5_DatasetName,
                                         MatrixContainer.GetMatrix<RealMatrix>(SampledMatrixID),
                                         MatrixContainer.GetMatrix<IndexMatrix>(sensor_mask_corners),
                                         ReductionOp,
                                         BufferToReuse);
  }

  return Stream;
}// end of CreateNewOutputStream
//------------------------------------------------------------------------------

//----------------------------------------------------------------------------//
//--------------------------- Private methods --------------------------------//
//----------------------------------------------------------------------------//
