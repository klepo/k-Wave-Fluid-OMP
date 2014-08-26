/**
 * @file        MatrixContainer.cpp
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The implementation file containing the matrix container
 *
 * @version     kspaceFirstOrder3D 2.15
 * @date        12 July     2012, 10:27  (created) \n
 *              26 August   2014, 16:40  (revised)
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
#include <MatrixClasses/OutputHDF5Stream.h>

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
 * Copy constructor of TMatrixRecord
 * @param [in] src
 */
TMatrixRecord::TMatrixRecord(const TMatrixRecord& src) :
        MatrixPtr(src.MatrixPtr),
        MatrixDataType(src.MatrixDataType),
        DimensionSizes(src.DimensionSizes),
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
    DimensionSizes  = src.DimensionSizes;
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
void TMatrixRecord::SetAllValues(TBaseMatrix *          MatrixPtr,
                                 const TMatrixDataType  MatrixDataType,
                                 const TDimensionSizes  DimensionSizes,
                                 const bool             LoadData,
                                 const bool             Checkpoint,
                                 const string           HDF5MatrixName)
{
  this->MatrixPtr        = MatrixPtr;
  this->MatrixDataType   = MatrixDataType;
  this->DimensionSizes   = DimensionSizes;
  this->LoadData         = LoadData;
  this->Checkpoint       = Checkpoint;
  this->HDF5MatrixName   = HDF5MatrixName;
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
 * Destructor of TMatrixContainer
 */
TMatrixContainer::~TMatrixContainer()
{
  MatrixContainer.clear();
}// end of ~TMatrixContainer
//------------------------------------------------------------------------------


/**
 * Create all matrix objects in the container.
 * \throw errors cause an exception bad_alloc.
 */
void TMatrixContainer::CreateAllObjects()
{
  for (TMatrixRecordContainer::iterator it = MatrixContainer.begin(); it != MatrixContainer.end(); it++)
  {
    if (it->second.MatrixPtr != NULL)
    {
      PrintErrorAndThrowException(MatrixContainer_ERR_FMT_ReloactaionError,
                                  it->second.HDF5MatrixName,
                                 __FILE__,__LINE__);
    }

    switch (it->second.MatrixDataType)
    {
      case TMatrixRecord::mdtReal:
      {
        it->second.MatrixPtr = new TRealMatrix(it->second.DimensionSizes);
        break;
      }

      case TMatrixRecord::mdtComplex:
      {
        it->second.MatrixPtr = new TComplexMatrix(it->second.DimensionSizes);
        break;
      }

      case TMatrixRecord::mdtIndex:
      {
        it->second.MatrixPtr = new TLongMatrix(it->second.DimensionSizes);
        break;
      }

      case TMatrixRecord::mdtFFTW:
      {
        it->second.MatrixPtr = new TFFTWComplexMatrix(it->second.DimensionSizes);
        break;
      }

      case TMatrixRecord::mdtUxyz:
      {
        it->second.MatrixPtr = new Tuxyz_sgxyzMatrix(it->second.DimensionSizes);
        break;
      }

      default:
      {
        PrintErrorAndThrowException(MatrixContainer_ERR_FMT_RecordUnknownDistributionType,
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
void TMatrixContainer::LoadDataFromInputHDF5File(THDF5_File & HDF5_File)
{
  for (TMatrixRecordContainer::iterator it = MatrixContainer.begin(); it != MatrixContainer.end(); it++)
  {
    if (it->second.LoadData)
    {
      it->second.MatrixPtr->ReadDataFromHDF5File(HDF5_File, it->second.HDF5MatrixName.c_str());
    }
  }
}// end of LoadMatricesDataFromDisk
//------------------------------------------------------------------------------


/**
 * Load selected matrices from checkpoint HDF5 file
 * @param [in] HDF5_File - HDF5 file handle
 */
void TMatrixContainer::LoadDataFromCheckpointHDF5File(THDF5_File & HDF5_File)
{
  for (TMatrixRecordContainer::iterator it = MatrixContainer.begin(); it != MatrixContainer.end(); it++)
  {
    if (it->second.Checkpoint)
    {
      it->second.MatrixPtr->ReadDataFromHDF5File(HDF5_File, it->second.HDF5MatrixName.c_str());
    }
  }
}// end of LoadDataFromCheckpointHDF5File
//------------------------------------------------------------------------------

/**
 * Store selected matrices into the checkpoint file
 * @param [in] HDF5_File
 */
void TMatrixContainer::StoreDataIntoCheckpointHDF5File(THDF5_File & HDF5_File)
{
  for (TMatrixRecordContainer::iterator it = MatrixContainer.begin(); it != MatrixContainer.end(); it++)
  {
    if (it->second.Checkpoint)
    {
      it->second.MatrixPtr->WriteDataToHDF5File(HDF5_File, it->second.HDF5MatrixName.c_str(),TParameters::GetInstance()->GetCompressionLevel());

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

  TParameters * Params = TParameters::GetInstance();

  TDimensionSizes FullDims = Params->GetFullDimensionSizes();
  TDimensionSizes ReducedDims = Params->GetReducedDimensionSizes();

  const bool LOAD         = true;
  const bool NOLOAD       = false;
  const bool CHECKPOINT   = true;
  const bool NOCHECKPOINT = false;


  //----------------------Allocate all matrices ----------------------------//

  MatrixContainer[kappa] .SetAllValues(NULL ,TMatrixRecord::mdtReal, ReducedDims, NOLOAD, NOCHECKPOINT, kappa_r_Name);
  if (!Params->Get_c0_scalar_flag())
  {
    MatrixContainer[c2]  .SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   ,   LOAD, NOCHECKPOINT, c0_Name);
  }
  MatrixContainer[p]     .SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   , NOLOAD,   CHECKPOINT, p_Name);

  MatrixContainer[rhox]  .SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   , NOLOAD,   CHECKPOINT, rhox_Name);
  MatrixContainer[rhoy]  .SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   , NOLOAD,   CHECKPOINT, rhoy_Name);
  MatrixContainer[rhoz]  .SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   , NOLOAD,   CHECKPOINT, rhoz_Name);

  MatrixContainer[ux_sgx].SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   , NOLOAD,   CHECKPOINT, ux_sgx_Name);
  MatrixContainer[uy_sgy].SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   , NOLOAD,   CHECKPOINT, uy_sgy_Name);
  MatrixContainer[uz_sgz].SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   , NOLOAD,   CHECKPOINT, uz_sgz_Name);

  MatrixContainer[duxdx] .SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   , NOLOAD, NOCHECKPOINT, duxdx_Name);
  MatrixContainer[duydy] .SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   , NOLOAD, NOCHECKPOINT, duydy_Name);
  MatrixContainer[duzdz] .SetAllValues(NULL ,TMatrixRecord::mdtReal, FullDims   , NOLOAD, NOCHECKPOINT, duzdz_Name);

  if (!Params->Get_rho0_scalar_flag())
  {
    MatrixContainer[rho0]       .SetAllValues(NULL,TMatrixRecord::mdtReal, FullDims,  LOAD, NOCHECKPOINT, rho0_Name);
    MatrixContainer[dt_rho0_sgx].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDims,  LOAD, NOCHECKPOINT, rho0_sgx_Name);
    MatrixContainer[dt_rho0_sgy].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDims,  LOAD, NOCHECKPOINT, rho0_sgy_Name);
    MatrixContainer[dt_rho0_sgz].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDims,  LOAD, NOCHECKPOINT, rho0_sgz_Name);
  }


  MatrixContainer[ddx_k_shift_pos].SetAllValues(NULL,TMatrixRecord::mdtComplex, TDimensionSizes(ReducedDims.X, 1, 1), LOAD, NOCHECKPOINT, ddx_k_shift_pos_r_Name);
  MatrixContainer[ddy_k_shift_pos].SetAllValues(NULL,TMatrixRecord::mdtComplex, TDimensionSizes(1, ReducedDims.Y, 1), LOAD, NOCHECKPOINT, ddy_k_shift_pos_Name);
  MatrixContainer[ddz_k_shift_pos].SetAllValues(NULL,TMatrixRecord::mdtComplex, TDimensionSizes(1, 1, ReducedDims.Z), LOAD, NOCHECKPOINT, ddz_k_shift_pos_Name);

  MatrixContainer[ddx_k_shift_neg].SetAllValues(NULL,TMatrixRecord::mdtComplex, TDimensionSizes(ReducedDims.X ,1, 1), LOAD, NOCHECKPOINT, ddx_k_shift_neg_r_Name);
  MatrixContainer[ddy_k_shift_neg].SetAllValues(NULL,TMatrixRecord::mdtComplex, TDimensionSizes(1, ReducedDims.Y, 1), LOAD, NOCHECKPOINT, ddy_k_shift_neg_Name);
  MatrixContainer[ddz_k_shift_neg].SetAllValues(NULL,TMatrixRecord::mdtComplex, TDimensionSizes(1, 1, ReducedDims.Z), LOAD, NOCHECKPOINT, ddz_k_shift_neg_Name);


  MatrixContainer[pml_x_sgx] .SetAllValues(NULL,TMatrixRecord::mdtReal,         TDimensionSizes(FullDims.X, 1, 1),    LOAD, NOCHECKPOINT, pml_x_sgx_Name);
  MatrixContainer[pml_y_sgy] .SetAllValues(NULL,TMatrixRecord::mdtReal,         TDimensionSizes(1, FullDims.Y, 1),    LOAD, NOCHECKPOINT, pml_y_sgy_Name);
  MatrixContainer[pml_z_sgz] .SetAllValues(NULL,TMatrixRecord::mdtReal,         TDimensionSizes(1, 1, FullDims.Z),    LOAD, NOCHECKPOINT, pml_z_sgz_Name);

  MatrixContainer[pml_x]     .SetAllValues(NULL,TMatrixRecord::mdtReal,         TDimensionSizes(FullDims.X, 1, 1),    LOAD, NOCHECKPOINT, pml_x_Name);
  MatrixContainer[pml_y]     .SetAllValues(NULL,TMatrixRecord::mdtReal,         TDimensionSizes(1, FullDims.Y, 1),    LOAD, NOCHECKPOINT, pml_y_Name);
  MatrixContainer[pml_z]     .SetAllValues(NULL,TMatrixRecord::mdtReal,         TDimensionSizes(1, 1, FullDims.Z),    LOAD, NOCHECKPOINT, pml_z_Name);

  if (Params->Get_nonlinear_flag())
  {
    if (! Params->Get_BonA_scalar_flag())
    {
      MatrixContainer[BonA]      .SetAllValues   (NULL,TMatrixRecord::mdtReal, FullDims ,   LOAD, NOCHECKPOINT, BonA_Name);
    }
  }

  if (Params->Get_absorbing_flag() != 0)
  {
    if (!((Params->Get_c0_scalar_flag()) && (Params->Get_alpha_coeff_scallar_flag())))
    {
      MatrixContainer[absorb_tau].SetAllValues   (NULL,TMatrixRecord::mdtReal, FullDims , NOLOAD, NOCHECKPOINT, absorb_tau_Name);
      MatrixContainer[absorb_eta].SetAllValues   (NULL,TMatrixRecord::mdtReal, FullDims , NOLOAD, NOCHECKPOINT, absorb_eta_Name);
    }
    MatrixContainer[absorb_nabla1].SetAllValues(NULL,TMatrixRecord::mdtReal, ReducedDims, NOLOAD, NOCHECKPOINT, absorb_nabla1_r_Name);
    MatrixContainer[absorb_nabla2].SetAllValues(NULL,TMatrixRecord::mdtReal, ReducedDims, NOLOAD, NOCHECKPOINT, absorb_nabla2_r_Name);
  }

  // linear sensor mask
  if (Params->Get_sensor_mask_type() == TParameters::smt_index)
  {
    MatrixContainer[sensor_mask_index].SetAllValues(NULL,TMatrixRecord::mdtIndex,
                                                    TDimensionSizes(Params->Get_sensor_mask_index_size(), 1, 1),
                                                    LOAD, NOCHECKPOINT, sensor_mask_index_Name);
  }

  // cuboid sensor mask
  if (Params->Get_sensor_mask_type() == TParameters::smt_corners)
  {
    MatrixContainer[sensor_mask_corners].SetAllValues(NULL,TMatrixRecord::mdtIndex,
                                                      TDimensionSizes(6 ,Params->Get_sensor_mask_corners_size(), 1),
                                                      LOAD, NOCHECKPOINT, sensor_mask_corners_Name);
  }


  // if p0 source flag
  if (Params->Get_p0_source_flag() == 1)
  {
    MatrixContainer[p0_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDims, LOAD, NOCHECKPOINT, p0_source_input_Name);
  }


  // us_index
  if ((Params->Get_transducer_source_flag() != 0) ||
      (Params->Get_ux_source_flag() != 0)         ||
      (Params->Get_uy_source_flag() != 0)         ||
      (Params->Get_uz_source_flag() != 0))
  {
    MatrixContainer[u_source_index].SetAllValues(NULL,TMatrixRecord::mdtIndex,
                                                 TDimensionSizes(1 ,1, Params->Get_u_source_index_size()),
                                                 LOAD, NOCHECKPOINT, u_source_index_Name);
  }

  //transducer source flag defined
  if (Params->Get_transducer_source_flag() != 0)
  {
    MatrixContainer[delay_mask]             .SetAllValues(NULL,TMatrixRecord::mdtIndex,TDimensionSizes(1 ,1, Params->Get_u_source_index_size())         ,
                                                          LOAD, NOCHECKPOINT, delay_mask_Name);
    MatrixContainer[transducer_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal ,TDimensionSizes(1 ,1, Params->Get_transducer_source_input_size()),
                                                          LOAD, NOCHECKPOINT, transducer_source_input_Name);
  }

  // p variables
  if (Params->Get_p_source_flag() != 0)
  {
    if (Params->Get_p_source_many() == 0)
    { // 1D case
      MatrixContainer[p_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal,
                                                   TDimensionSizes(1 ,1, Params->Get_p_source_flag()),
                                                   LOAD, NOCHECKPOINT, p_source_input_Name);
    }
    else
    { // 2D case
      MatrixContainer[p_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal,
                                                   TDimensionSizes(1 ,Params->Get_p_source_index_size(),Params->Get_p_source_flag()),
                                                   LOAD, NOCHECKPOINT, p_source_input_Name);
    }

    MatrixContainer[p_source_index].SetAllValues(NULL,TMatrixRecord::mdtIndex,
                                                 TDimensionSizes(1 ,1, Params->Get_p_source_index_size()),
                                                 LOAD, NOCHECKPOINT, p_source_index_Name);
  }



    //----------------------------uxyz source flags---------------------------//
  if (Params->Get_ux_source_flag() != 0)
  {
    if (Params->Get_u_source_many() == 0)
    { // 1D
      MatrixContainer[ux_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal,
                                                    TDimensionSizes(1 ,1, Params->Get_ux_source_flag()),
                                                    LOAD, NOCHECKPOINT, ux_source_input_Name);
    }
    else
    { // 2D
      MatrixContainer[ux_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal,
                                                    TDimensionSizes(1 ,Params->Get_u_source_index_size(),Params->Get_ux_source_flag()),
                                                    LOAD, NOCHECKPOINT, ux_source_input_Name);
    }
  }// ux_source_input


  if (Params->Get_uy_source_flag() != 0)
  {
    if (Params->Get_u_source_many() == 0)
    { // 1D
      MatrixContainer[uy_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal,
                                                    TDimensionSizes(1 ,1, Params->Get_uy_source_flag()),
                                                    LOAD, NOCHECKPOINT, uy_source_input_Name);
    }
    else
    { // 2D
      MatrixContainer[uy_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal,
                                                    TDimensionSizes(1 ,Params->Get_u_source_index_size(),Params->Get_uy_source_flag()),
                                                    LOAD, NOCHECKPOINT, uy_source_input_Name);
    }
  }// uy_source_input

  if (Params->Get_uz_source_flag() != 0)
  {
    if (Params->Get_u_source_many() == 0)
    { // 1D
      MatrixContainer[uz_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal,
                                                    TDimensionSizes(1 ,1, Params->Get_uz_source_flag()),
                                                    LOAD, NOCHECKPOINT, uz_source_input_Name);
    }
    else
    { // 2D
      MatrixContainer[uz_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal,
                                                    TDimensionSizes(1 ,Params->Get_u_source_index_size(),Params->Get_uz_source_flag()),
                                                    LOAD, NOCHECKPOINT, uz_source_input_Name);
    }
  }// uz_source_input




  //-- Nonlinear grid
  if (Params->Get_nonuniform_grid_flag()!= 0)
  {
    MatrixContainer[dxudxn]    .SetAllValues(NULL,TMatrixRecord::mdtReal, TDimensionSizes(FullDims.X, 1, 1), LOAD, NOCHECKPOINT, dxudxn_Name);
    MatrixContainer[dyudyn]    .SetAllValues(NULL,TMatrixRecord::mdtReal, TDimensionSizes(1, FullDims.Y, 1), LOAD, NOCHECKPOINT, dyudyn_Name);
    MatrixContainer[dzudzn]    .SetAllValues(NULL,TMatrixRecord::mdtReal, TDimensionSizes(1 ,1, FullDims.Z), LOAD, NOCHECKPOINT, dzudzn_Name);

    MatrixContainer[dxudxn_sgx].SetAllValues(NULL,TMatrixRecord::mdtReal, TDimensionSizes(FullDims.X, 1, 1), LOAD, NOCHECKPOINT, dxudxn_sgx_Name);
    MatrixContainer[dyudyn_sgy].SetAllValues(NULL,TMatrixRecord::mdtReal, TDimensionSizes(1, FullDims.Y, 1), LOAD, NOCHECKPOINT, dyudyn_sgy_Name);
    MatrixContainer[dzudzn_sgz].SetAllValues(NULL,TMatrixRecord::mdtReal, TDimensionSizes(1 ,1, FullDims.Z), LOAD, NOCHECKPOINT, dzudzn_sgz_Name);
  }

  //------------------------------------------------------------------------//
  //--------------------- Temporary matrices -------------------------------//
  //------------------------------------------------------------------------//
  // this matrix used to load alpha_coeff for absorb_tau pre-calculation

  if ((Params->Get_absorbing_flag() != 0) && (!Params->Get_alpha_coeff_scallar_flag()))
  {
    MatrixContainer[Temp_1_RS3D].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDims ,   LOAD, NOCHECKPOINT, alpha_coeff_Name);
  }
  else
  {
    MatrixContainer[Temp_1_RS3D].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDims , NOLOAD, NOCHECKPOINT, "");
  }

  MatrixContainer[Temp_2_RS3D].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDims   , NOLOAD, NOCHECKPOINT, "");
  MatrixContainer[Temp_3_RS3D].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDims   , NOLOAD, NOCHECKPOINT, "");

  MatrixContainer[FFT_X_temp].SetAllValues(NULL, TMatrixRecord::mdtFFTW, ReducedDims, NOLOAD, NOCHECKPOINT, "");
  MatrixContainer[FFT_Y_temp].SetAllValues(NULL, TMatrixRecord::mdtFFTW, ReducedDims, NOLOAD, NOCHECKPOINT, "");
  MatrixContainer[FFT_Z_temp].SetAllValues(NULL, TMatrixRecord::mdtFFTW, ReducedDims, NOLOAD, NOCHECKPOINT, "");


  //-- u_non_staggered_raw
  if (Params->IsStore_u_non_staggered_raw())
  {
    TDimensionSizes ShiftDims = FullDims;

    if ((FullDims.X > FullDims.Y) && (FullDims.X > FullDims.Z))
    { // X is the longest
      ShiftDims.X = FullDims.X / 2 + 1;
    }
    else if ((FullDims.Y > FullDims.X) && (FullDims.Y > FullDims.Z))
    { // Y is the longest
      ShiftDims.Y = FullDims.Y / 2 + 1;
    }
    else if ((FullDims.Z > FullDims.X) && (FullDims.Z > FullDims.Y))
    { // Z is the longest
      ShiftDims.Z = FullDims.Z / 2 + 1;
    }
    else
    { //all are the same
      ShiftDims.X = FullDims.X / 2 + 1;
    }

    MatrixContainer[FFT_shift_temp].SetAllValues(NULL, TMatrixRecord::mdtFFTW, ShiftDims, NOLOAD, NOCHECKPOINT, "");

    // these three are necessary only for u_non_staggered calculation now
    MatrixContainer[ux_shifted].SetAllValues(NULL, TMatrixRecord::mdtReal, FullDims, NOLOAD, NOCHECKPOINT, "");
    MatrixContainer[uy_shifted].SetAllValues(NULL, TMatrixRecord::mdtReal, FullDims, NOLOAD, NOCHECKPOINT, "");
    MatrixContainer[uz_shifted].SetAllValues(NULL, TMatrixRecord::mdtReal, FullDims, NOLOAD, NOCHECKPOINT, "");
  }// u_non_staggered

}// end of InitMatrixContainer
//------------------------------------------------------------------------------


//----------------------------------------------------------------------------//
//------------------------- Protected methods --------------------------------//
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
//-------------------------- Private methods ---------------------------------//
//----------------------------------------------------------------------------//


/**
 * Print error and and throw an exception
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
 * Please note, the Matrixcontainer has to be populated before calling this routine
 *
 * @param [in] MatrixContainer - matrix container to link the steams with
 *                               sampled matrices and sensor masks
 */
void TOutputStreamContainer::AddStreamsIntoContainer(TMatrixContainer & MatrixContainer)
{

  TParameters * Params = TParameters::GetInstance();

  float * TempBufferX = MatrixContainer.GetRealMatrix(Temp_1_RS3D).GetRawData();
  float * TempBufferY = MatrixContainer.GetRealMatrix(Temp_2_RS3D).GetRawData();
  float * TempBufferZ = MatrixContainer.GetRealMatrix(Temp_3_RS3D).GetRawData();


  //--------------------- Pressure ------------------/
  if (Params->IsStore_p_raw())
  {
    OutputStreamContainer[p_sensor_raw] = CreateNewOutputStream(MatrixContainer,
                                                                p,
                                                                p_Name,
                                                                TBaseOutputHDF5Stream::roNONE,
                                                                TempBufferX);
  }// IsStore_p_raw

  if (Params->IsStore_p_rms())
  {
    OutputStreamContainer[p_sensor_rms] = CreateNewOutputStream(MatrixContainer,
                                                                p,
                                                                p_rms_Name,
                                                                TBaseOutputHDF5Stream::roRMS);
  }

  if (Params->IsStore_p_max())
  {
    OutputStreamContainer[p_sensor_max] = CreateNewOutputStream(MatrixContainer,
                                                                p,
                                                                p_max_Name,
                                                                TBaseOutputHDF5Stream::roMAX);
  }

  if (Params->IsStore_p_min())
  {
    OutputStreamContainer[p_sensor_min] = CreateNewOutputStream(MatrixContainer,
                                                                p,
                                                                p_min_Name,
                                                                TBaseOutputHDF5Stream::roMIN);
  }

  if (Params->IsStore_p_max_all())
  {
    OutputStreamContainer[p_sensor_max_all] =
            new TWholeDomainOutputHDF5Stream(Params->HDF5_OutputFile,
                                             p_max_all_Name,
                                             MatrixContainer.GetRealMatrix(p),
                                             TBaseOutputHDF5Stream::roMAX);
  }

  if (Params->IsStore_p_min_all())
  {
    OutputStreamContainer[p_sensor_min_all] =
            new TWholeDomainOutputHDF5Stream(Params->HDF5_OutputFile,
                                             p_min_all_Name,
                                             MatrixContainer.GetRealMatrix(p),
                                             TBaseOutputHDF5Stream::roMIN);
  }

  //--------------------- Velocity ------------------/
  if (Params->IsStore_u_raw())
  {
    OutputStreamContainer[ux_sensor_raw] = CreateNewOutputStream(MatrixContainer,
                                                                 ux_sgx,
                                                                 ux_Name,
                                                                 TBaseOutputHDF5Stream::roNONE,
                                                                 TempBufferX);
    OutputStreamContainer[uy_sensor_raw] = CreateNewOutputStream(MatrixContainer,
                                                                 uy_sgy,
                                                                 uy_Name,
                                                                 TBaseOutputHDF5Stream::roNONE,
                                                                 TempBufferY);
    OutputStreamContainer[uz_sensor_raw] = CreateNewOutputStream(MatrixContainer,
                                                                 uz_sgz,
                                                                 uz_Name,
                                                                 TBaseOutputHDF5Stream::roNONE,
                                                                 TempBufferZ);
  }

  if (Params->IsStore_u_rms())
  {
    OutputStreamContainer[ux_sensor_rms] = CreateNewOutputStream(MatrixContainer,
                                                                 ux_sgx,
                                                                 ux_rms_Name,
                                                                 TBaseOutputHDF5Stream::roRMS);
    OutputStreamContainer[uy_sensor_rms] = CreateNewOutputStream(MatrixContainer,
                                                                 uy_sgy,
                                                                 uy_rms_Name,
                                                                 TBaseOutputHDF5Stream::roRMS);
    OutputStreamContainer[uz_sensor_rms] = CreateNewOutputStream(MatrixContainer,
                                                                 uz_sgz,
                                                                 uz_rms_Name,
                                                                 TBaseOutputHDF5Stream::roRMS);
  }

  if (Params->IsStore_u_max())
  {
    OutputStreamContainer[ux_sensor_max] = CreateNewOutputStream(MatrixContainer,
                                                                 ux_sgx,
                                                                 ux_max_Name,
                                                                 TBaseOutputHDF5Stream::roMAX);
    OutputStreamContainer[uy_sensor_max] = CreateNewOutputStream(MatrixContainer,
                                                                 uy_sgy,
                                                                 uy_max_Name,
                                                                 TBaseOutputHDF5Stream::roMAX);
    OutputStreamContainer[uz_sensor_max] = CreateNewOutputStream(MatrixContainer,
                                                                 uz_sgz,
                                                                 uz_max_Name,
                                                                 TBaseOutputHDF5Stream::roMAX);
  }

  if (Params->IsStore_u_min())
  {
    OutputStreamContainer[ux_sensor_min] = CreateNewOutputStream(MatrixContainer,
                                                                 ux_sgx,
                                                                 ux_min_Name,
                                                                 TBaseOutputHDF5Stream::roMIN);
    OutputStreamContainer[uy_sensor_min] = CreateNewOutputStream(MatrixContainer,
                                                                 uy_sgy,
                                                                 uy_min_Name,
                                                                 TBaseOutputHDF5Stream::roMIN);
    OutputStreamContainer[uz_sensor_min] = CreateNewOutputStream(MatrixContainer,
                                                                 uz_sgz,
                                                                 uz_min_Name,
                                                                 TBaseOutputHDF5Stream::roMIN);
  }

  if (Params->IsStore_u_max_all())
  {
    OutputStreamContainer[ux_sensor_max_all] =
            new TWholeDomainOutputHDF5Stream(Params->HDF5_OutputFile,
                                             ux_max_all_Name,
                                             MatrixContainer.GetRealMatrix(ux_sgx),
                                             TBaseOutputHDF5Stream::roMAX);
    OutputStreamContainer[uy_sensor_max_all] =
            new TWholeDomainOutputHDF5Stream(Params->HDF5_OutputFile,
                                             uy_max_all_Name,
                                             MatrixContainer.GetRealMatrix(uy_sgy),
                                             TBaseOutputHDF5Stream::roMAX);
    OutputStreamContainer[uz_sensor_max_all] =
            new TWholeDomainOutputHDF5Stream(Params->HDF5_OutputFile,
                                             uz_max_all_Name,
                                             MatrixContainer.GetRealMatrix(uz_sgz),
                                             TBaseOutputHDF5Stream::roMAX);
  }

  if (Params->IsStore_u_min_all())
  {
    OutputStreamContainer[ux_sensor_min_all] =
            new TWholeDomainOutputHDF5Stream(Params->HDF5_OutputFile,
                                             ux_min_all_Name,
                                             MatrixContainer.GetRealMatrix(ux_sgx),
                                             TBaseOutputHDF5Stream::roMIN);
    OutputStreamContainer[uy_sensor_min_all] =
            new TWholeDomainOutputHDF5Stream(Params->HDF5_OutputFile,
                                             uy_min_all_Name,
                                             MatrixContainer.GetRealMatrix(uy_sgy),
                                             TBaseOutputHDF5Stream::roMIN);
    OutputStreamContainer[uz_sensor_min_all] =
            new TWholeDomainOutputHDF5Stream(Params->HDF5_OutputFile,
                                             uz_min_all_Name,
                                             MatrixContainer.GetRealMatrix(uz_sgz),
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
}// end of CreateStreams
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
}// end of CloseStreams
//------------------------------------------------------------------------------


/**
 * Checkpoint streams without post-processing (flush to the file)
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
 * /// Post-process all streams and flush them to the file
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
 * Create a new output stream
 * @param [in] MatrixContainer  - name of the HDF5 dataset or group
 * @param [in] SampledMatrixID  - code id of the matrix
 * @param [in] HDF5_DatasetName - name of the HDF5 dataset or group
 * @param [in] ReductionOp      - reduction operator
 * @param [in] BufferToReuse   - buffer to reuse
 * @return - new output stream with defined links
 */
TBaseOutputHDF5Stream * TOutputStreamContainer::CreateNewOutputStream(TMatrixContainer & MatrixContainer,
                                                                      const TMatrixID       SampledMatrixID,
                                                                      const char *          HDF5_DatasetName,
                                                                      const TBaseOutputHDF5Stream::TReductionOperator  ReductionOp,
                                                                      float *               BufferToReuse)
{

  TParameters * Params = TParameters::GetInstance();

  TBaseOutputHDF5Stream * Stream = NULL;

  if (Params->Get_sensor_mask_type() == TParameters::smt_index)
  {

    Stream = new TIndexOutputHDF5Stream(Params->HDF5_OutputFile,
                                      HDF5_DatasetName,
                                      MatrixContainer.GetRealMatrix(SampledMatrixID),
                                      MatrixContainer.GetLongMatrix(sensor_mask_index),
                                      ReductionOp,
                                      BufferToReuse);
  }
  else
  {
    Stream = new TCuboidOutputHDF5Stream(Params->HDF5_OutputFile,
                                      HDF5_DatasetName,
                                      MatrixContainer.GetRealMatrix(SampledMatrixID),
                                      MatrixContainer.GetLongMatrix(sensor_mask_corners),
                                      ReductionOp,
                                      BufferToReuse);
  }

  return Stream;
}// end of CreateNewOutputStream





//----------------------------------------------------------------------------//
//--------------------------- Private methods --------------------------------//
//----------------------------------------------------------------------------//