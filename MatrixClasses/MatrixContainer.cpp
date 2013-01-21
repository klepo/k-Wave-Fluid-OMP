/**
 * @file        MatrixContainer.cpp
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au
 * 
 * @brief       The implementation file containing the matrix container
 * 
 * @version     kspaceFirstOrder3D 2.13
 * @date        12 July 2012, 10:27  (created) \n
 *              14 September 2012, 14:20(revised)
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
 * Copy constructor of TMatrixRecord
 * @param [in] src
 */
TMatrixRecord::TMatrixRecord(const TMatrixRecord& src) :
    MatrixPtr(src.MatrixPtr), MatrixDataType(src.MatrixDataType),             
    DimensionSizes(src.DimensionSizes), LoadData(src.LoadData),  
    HDF5MatrixName(src.HDF5MatrixName)
{
    

}// end of TMatrixRecord
//------------------------------------------------------------------------------


/**
 * operator = of TMatrixRecord
 * @param  [in] src
 * @return this
 */
TMatrixRecord& TMatrixRecord::operator = (const TMatrixRecord& src){

    if (this != &src){
        MatrixPtr = src.MatrixPtr;           
        MatrixDataType  = src.MatrixDataType;            
        DimensionSizes  = src.DimensionSizes;
        LoadData        = src.LoadData;
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
 * @param [in] HDF5MatrixName   - HDF5 matrix name
 */
void TMatrixRecord::SetAllValues(TBaseMatrix *          MatrixPtr,                                                                   
                                 const TMatrixDataType  MatrixDataType,                                                                   
                                 const TDimensionSizes  DimensionSizes,
                                 const bool             LoadData, 
                                 const string           HDF5MatrixName){
    
    this->MatrixPtr        = MatrixPtr;    
    this->MatrixDataType   = MatrixDataType;        
    this->DimensionSizes   = DimensionSizes;    
    this->LoadData         = LoadData;
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
TMatrixContainer::~TMatrixContainer(){
    
    MatrixContainer.clear();    
    
}// end of ~TMatrixContainer
//------------------------------------------------------------------------------

 
/**
 * Create all matrix objects in the container.
 * \throw errors cause an exception bad_alloc.
 */
void TMatrixContainer::CreateAllObjects(){
    
    
    for (TMatrixRecordContainer::iterator it = MatrixContainer.begin(); it != MatrixContainer.end(); it++){
                
        if (it->second.MatrixPtr != NULL) {
            PrintErrorAndThrowException(MatrixContainer_ERR_FMT_ReloactaionError, it->second.HDF5MatrixName,
                               __FILE__,__LINE__);                    
        }                       
        
        switch (it->second.MatrixDataType){                        
            
            case TMatrixRecord::mdtReal:  {
                 it->second.MatrixPtr =  new TRealMatrix(it->second.DimensionSizes);                    
                 break;
            }
            
            case TMatrixRecord::mdtComplex:  {
                 it->second.MatrixPtr =  new TComplexMatrix(it->second.DimensionSizes);                    
                 break;
            }
            
            case TMatrixRecord::mdtIndex:  {
                 it->second.MatrixPtr =  new TLongMatrix(it->second.DimensionSizes);                    
                 break;
            }
            case TMatrixRecord::mdtFFTW:  {
                 it->second.MatrixPtr =  new TFFTWComplexMatrix(it->second.DimensionSizes);                    
                 break;
            }
            
            case TMatrixRecord::mdtUxyz:  {
                 it->second.MatrixPtr =  new Tuxyz_sgxyzMatrix(it->second.DimensionSizes);                    
                 break;
            }
                       
            default:{
                PrintErrorAndThrowException(MatrixContainer_ERR_FMT_RecordUnknownDistributionType, it->second.HDF5MatrixName,
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
void TMatrixContainer::LoadMatricesDataFromDisk(THDF5_File & HDF5_File){
   
    for (TMatrixRecordContainer::iterator it = MatrixContainer.begin(); it != MatrixContainer.end(); it++){
      
        if (it->second.LoadData) {
                    
            it->second.MatrixPtr->ReadDataFromHDF5File(HDF5_File, it->second.HDF5MatrixName.c_str());        
        
        }
    }       
    
}// end of LoadMatricesDataFromDisk
//------------------------------------------------------------------------------

  
/**
 * Free all matrix objects.
 * 
 */
void TMatrixContainer::FreeAllMatrices(){
     for (TMatrixRecordContainer::iterator it = MatrixContainer.begin(); it != MatrixContainer.end(); it++){
        if (it->second.MatrixPtr){
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
void TMatrixContainer::AddMatricesIntoContainer(){
    
    TParameters * Params = TParameters::GetInstance();
        
    TDimensionSizes FullDim = Params->GetFullDimensionSizes();
    TDimensionSizes ReducedDim = Params->GetReducedDimensionSizes();
    
    //----------------------Allocate all matrices ----------------------------//
                                             
    MatrixContainer[kappa] .SetAllValues(NULL ,TMatrixRecord::mdtReal  , ReducedDim, false, kappa_r_Name);
    if (!Params->Get_c0_scalar_flag()) {
        MatrixContainer[c2]    .SetAllValues(NULL ,TMatrixRecord::mdtReal  , FullDim   , true,  c0_Name);
    }
    MatrixContainer[p]     .SetAllValues(NULL ,TMatrixRecord::mdtReal  , FullDim   , false, p_Name);
     
    MatrixContainer[rhox]  .SetAllValues(NULL ,TMatrixRecord::mdtReal  , FullDim   , false, rhox_Name);
    MatrixContainer[rhoy]  .SetAllValues(NULL ,TMatrixRecord::mdtReal  , FullDim   , false, rhoy_Name);
    MatrixContainer[rhoz]  .SetAllValues(NULL ,TMatrixRecord::mdtReal  , FullDim   , false, rhoz_Name);    
    
    MatrixContainer[ux_sgx].SetAllValues(NULL ,TMatrixRecord::mdtReal  , FullDim   , false, ux_sgx_Name);
    MatrixContainer[uy_sgy].SetAllValues(NULL ,TMatrixRecord::mdtReal  , FullDim   , false, uy_sgy_Name);
    MatrixContainer[uz_sgz].SetAllValues(NULL ,TMatrixRecord::mdtReal  , FullDim   , false, uz_sgz_Name);
    
    MatrixContainer[duxdx] .SetAllValues(NULL ,TMatrixRecord::mdtReal  , FullDim   , false, duxdx_Name);    
    MatrixContainer[duydy] .SetAllValues(NULL ,TMatrixRecord::mdtReal  , FullDim   , false, duydy_Name);
    MatrixContainer[duzdz] .SetAllValues(NULL ,TMatrixRecord::mdtReal  , FullDim   , false, duzdz_Name);
    
    if (!Params->Get_rho0_scalar_flag()) {
       MatrixContainer[rho0]       .SetAllValues(NULL,TMatrixRecord::mdtReal, FullDim , true, rho0_Name);
       MatrixContainer[dt_rho0_sgx].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDim , true, rho0_sgx_Name);
       MatrixContainer[dt_rho0_sgy].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDim , true, rho0_sgy_Name);
       MatrixContainer[dt_rho0_sgz].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDim , true, rho0_sgz_Name);
    }   
    

    MatrixContainer[ddx_k_shift_pos].SetAllValues(NULL,TMatrixRecord::mdtComplex, TDimensionSizes(ReducedDim.X, 1, 1), true, ddx_k_shift_pos_r_Name);
    MatrixContainer[ddy_k_shift_pos].SetAllValues(NULL,TMatrixRecord::mdtComplex, TDimensionSizes(1, ReducedDim.Y, 1), true, ddy_k_shift_pos_Name);
    MatrixContainer[ddz_k_shift_pos].SetAllValues(NULL,TMatrixRecord::mdtComplex, TDimensionSizes(1, 1, ReducedDim.Z), true, ddz_k_shift_pos_Name);

    MatrixContainer[ddx_k_shift_neg].SetAllValues(NULL,TMatrixRecord::mdtComplex, TDimensionSizes(ReducedDim.X ,1, 1), true, ddx_k_shift_neg_r_Name);
    MatrixContainer[ddy_k_shift_neg].SetAllValues(NULL,TMatrixRecord::mdtComplex, TDimensionSizes(1, ReducedDim.Y, 1), true, ddy_k_shift_neg_Name);
    MatrixContainer[ddz_k_shift_neg].SetAllValues(NULL,TMatrixRecord::mdtComplex, TDimensionSizes(1, 1, ReducedDim.Z), true, ddz_k_shift_neg_Name);

    
    MatrixContainer[pml_x_sgx] .SetAllValues(NULL,TMatrixRecord::mdtReal,         TDimensionSizes(FullDim.X, 1, 1), true, pml_x_sgx_Name);
    MatrixContainer[pml_y_sgy] .SetAllValues(NULL,TMatrixRecord::mdtReal,         TDimensionSizes(1, FullDim.Y, 1), true, pml_y_sgy_Name);
    MatrixContainer[pml_z_sgz] .SetAllValues(NULL,TMatrixRecord::mdtReal,         TDimensionSizes(1, 1, FullDim.Z), true, pml_z_sgz_Name);

    MatrixContainer[pml_x]     .SetAllValues(NULL,TMatrixRecord::mdtReal,         TDimensionSizes(FullDim.X, 1, 1), true, pml_x_Name);
    MatrixContainer[pml_y]     .SetAllValues(NULL,TMatrixRecord::mdtReal,         TDimensionSizes(1, FullDim.Y, 1), true, pml_y_Name);
    MatrixContainer[pml_z]     .SetAllValues(NULL,TMatrixRecord::mdtReal,         TDimensionSizes(1, 1, FullDim.Z), true, pml_z_Name);
    
    if (Params->Get_nonlinear_flag()) {
        if (! Params->Get_BonA_scalar_flag()) {
            MatrixContainer[BonA]      .SetAllValues   (NULL,TMatrixRecord::mdtReal, FullDim   , true, BonA_Name);                                
        }
    }

    if (Params->Get_absorbing_flag() != 0) {
        if (!((Params->Get_c0_scalar_flag()) && (Params->Get_alpha_coeff_scallar_flag()))) {
            MatrixContainer[absorb_tau].SetAllValues   (NULL,TMatrixRecord::mdtReal, FullDim   , false, absorb_tau_Name);                                
            MatrixContainer[absorb_eta].SetAllValues   (NULL,TMatrixRecord::mdtReal, FullDim   , false, absorb_eta_Name);                                
        }
        MatrixContainer[absorb_nabla1].SetAllValues(NULL,TMatrixRecord::mdtReal, ReducedDim, false, absorb_nabla1_r_Name);                                
        MatrixContainer[absorb_nabla2].SetAllValues(NULL,TMatrixRecord::mdtReal, ReducedDim, false, absorb_nabla2_r_Name
        );                                
    }
        
    //-- index matrices --//                
    MatrixContainer[sensor_mask_ind].SetAllValues(NULL,TMatrixRecord::mdtIndex, TDimensionSizes(1 ,1, Params->Get_sensor_mask_index_size()), true, sensor_mask_index_Name);                                
   
    
    
    // if p0 source flag 
    if (Params->Get_p0_source_flag() == 1){
       MatrixContainer[p0_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDim, true, p0_source_input_Name);                                
    }
    
    
    // us_index    
    if ((Params->Get_transducer_source_flag() != 0) || 
        (Params->Get_ux_source_flag() != 0)         ||
        (Params->Get_uy_source_flag() != 0)         ||
        (Params->Get_uz_source_flag() != 0)            ){
                
                MatrixContainer[u_source_index].SetAllValues(NULL,TMatrixRecord::mdtIndex,TDimensionSizes(1 ,1, Params->Get_u_source_index_size()), true, u_source_index_Name);                                           
    }
                
                    
        // -- transducer source flag defined
    if (Params->Get_transducer_source_flag() != 0) {
            
        MatrixContainer[delay_mask]             .SetAllValues(NULL,TMatrixRecord::mdtIndex,TDimensionSizes(1 ,1, Params->Get_u_source_index_size())          , true, delay_mask_Name);                                                   
        MatrixContainer[transducer_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal ,TDimensionSizes(1 ,1, Params->Get_transducer_source_input_size()), true, transducer_source_input_Name);                                           
     
    }

    
      //-- p variables --// 
    if (Params->Get_p_source_flag() != 0){
                
           if (Params->Get_p_source_many() == 0) {    // 1D case               
                
               MatrixContainer[p_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal, TDimensionSizes(1 ,1, Params->Get_p_source_flag()), true, p_source_input_Name);                                           
               
           } else {                                             // 2D case

               MatrixContainer[p_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal, TDimensionSizes(1 ,Params->Get_p_source_index_size(), Params->Get_p_source_flag()), true, p_source_input_Name);                                           
           }   
           
    
           MatrixContainer[p_source_index].SetAllValues(NULL,TMatrixRecord::mdtIndex, TDimensionSizes(1 ,1, Params->Get_p_source_index_size()), true, p_source_index_Name);                                                      
    }

    
    
    //----------------------------uxyz source flags---------------------------//
    if (Params->Get_ux_source_flag() != 0) {
        if (Params->Get_u_source_many() == 0) { // 1D
            
            MatrixContainer[ux_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal, TDimensionSizes(1 ,1, Params->Get_ux_source_flag()), true, ux_source_input_Name);                                           
        }
        else {                                        // 2D
            MatrixContainer[ux_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal, TDimensionSizes(1 ,Params->Get_u_source_index_size(),Params->Get_ux_source_flag()), true, ux_source_input_Name);                                           
            
        }                                 
    }// ux_source_input
    
    
    if (Params->Get_uy_source_flag() != 0) {
        if (Params->Get_u_source_many() == 0) { // 1D
            
            MatrixContainer[uy_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal, TDimensionSizes(1 ,1, Params->Get_uy_source_flag()), true, uy_source_input_Name);                                           
        }
        else {                                        // 2D
            MatrixContainer[uy_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal, TDimensionSizes(1 ,Params->Get_u_source_index_size(),Params->Get_uy_source_flag()), true, uy_source_input_Name);                                           
            
        }                                 
                        
    }// uy_source_input
    
    if (Params->Get_uz_source_flag() != 0) {
        if (Params->Get_u_source_many() == 0) { // 1D
            
            MatrixContainer[uz_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal,TDimensionSizes(1 ,1, Params->Get_uz_source_flag()), true, uz_source_input_Name); 
        }
        else {                                        // 2D
            MatrixContainer[uz_source_input].SetAllValues(NULL,TMatrixRecord::mdtReal,TDimensionSizes(1 ,Params->Get_u_source_index_size(),Params->Get_uz_source_flag()), true, uz_source_input_Name);
            
        }                                     
        
        
    }// uz_source_input
    
    
    
    
    //-- Non linear grid
    if (Params->Get_nonuniform_grid_flag()!= 0) {            
        MatrixContainer[dxudxn].SetAllValues(NULL,TMatrixRecord::mdtReal, TDimensionSizes(FullDim.X, 1, 1), true, dxudxn_Name);                                            
        MatrixContainer[dyudyn].SetAllValues(NULL,TMatrixRecord::mdtReal, TDimensionSizes(1, FullDim.Y, 1), true, dyudyn_Name);                                            
        MatrixContainer[dzudzn].SetAllValues(NULL,TMatrixRecord::mdtReal, TDimensionSizes(1 ,1, FullDim.Z), true, dzudzn_Name);                                            

        
        MatrixContainer[dxudxn_sgx].SetAllValues(NULL,TMatrixRecord::mdtReal, TDimensionSizes(FullDim.X, 1, 1), true, dxudxn_sgx_Name);                                            
        MatrixContainer[dyudyn_sgy].SetAllValues(NULL,TMatrixRecord::mdtReal, TDimensionSizes(1, FullDim.Y, 1), true, dyudyn_sgy_Name);
        MatrixContainer[dzudzn_sgz].SetAllValues(NULL,TMatrixRecord::mdtReal, TDimensionSizes(1 ,1, FullDim.Z), true, dzudzn_sgz_Name);
        
    }

    
    //--Temporary matrices --//        
    // this matrix used to load alpha_coeff for absorb_tau pre-calculation    
    
    if ((Params->Get_absorbing_flag() != 0) && (!Params->Get_alpha_coeff_scallar_flag())){                
        MatrixContainer[Temp_1_RS3D].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDim, true, alpha_coeff_Name);         
    }else{
        MatrixContainer[Temp_1_RS3D].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDim, false, "");
    }

        
    MatrixContainer[Temp_2_RS3D].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDim, false, "");                                
    MatrixContainer[Temp_3_RS3D].SetAllValues(NULL,TMatrixRecord::mdtReal, FullDim, false, "");                                
    
    MatrixContainer[FFT_X_temp].SetAllValues(NULL, TMatrixRecord::mdtFFTW, ReducedDim, false, "");                                
    MatrixContainer[FFT_Y_temp].SetAllValues(NULL, TMatrixRecord::mdtFFTW, ReducedDim, false, "");                                
    MatrixContainer[FFT_Z_temp].SetAllValues(NULL, TMatrixRecord::mdtFFTW, ReducedDim, false, "");                                
    
    //--- output buffers --//
    
    
    TDimensionSizes SensorDims(Params->Get_sensor_mask_index_size(), 1, 1);
            
            
    if (Params->IsStore_p_rms()){
       MatrixContainer[p_sensor_rms].SetAllValues(NULL,TMatrixRecord::mdtReal, SensorDims, false, p_rms_Name);                                
    }
        
    if (Params->IsStore_p_max()){
       MatrixContainer[p_sensor_max].SetAllValues(NULL,TMatrixRecord::mdtReal, SensorDims, false, p_max_Name);                                
    }
                
    
    if (Params->IsStore_u_rms()){
       MatrixContainer[ux_sensor_rms].SetAllValues(NULL,TMatrixRecord::mdtReal, SensorDims, false, ux_rms_Name);
       MatrixContainer[uy_sensor_rms].SetAllValues(NULL,TMatrixRecord::mdtReal, SensorDims, false, uy_rms_Name);
       MatrixContainer[uz_sensor_rms].SetAllValues(NULL,TMatrixRecord::mdtReal, SensorDims, false, uz_rms_Name);
    }
        
    if (Params->IsStore_u_max()){
       MatrixContainer[ux_sensor_max].SetAllValues(NULL,TMatrixRecord::mdtReal, SensorDims, false, ux_max_Name);
       MatrixContainer[uy_sensor_max].SetAllValues(NULL,TMatrixRecord::mdtReal, SensorDims, false, uy_max_Name);                                
       MatrixContainer[uz_sensor_max].SetAllValues(NULL,TMatrixRecord::mdtReal, SensorDims, false, uz_max_Name);                                       
    }
            
    if (Params->IsStore_I_avg()){
       MatrixContainer[Ix_sensor_avg].SetAllValues(NULL,TMatrixRecord::mdtReal, SensorDims, false, Ix_avg_Name);                                
       MatrixContainer[Iy_sensor_avg].SetAllValues(NULL,TMatrixRecord::mdtReal, SensorDims, false, Iy_avg_Name);                                
       MatrixContainer[Iz_sensor_avg].SetAllValues(NULL,TMatrixRecord::mdtReal, SensorDims, false, Iz_avg_Name);                                
       
    }
        
    if (Params->IsStore_I_max()){
       MatrixContainer[Ix_sensor_max].SetAllValues(NULL,TMatrixRecord::mdtReal, SensorDims, false, Ix_max_Name);
       MatrixContainer[Iy_sensor_max].SetAllValues(NULL,TMatrixRecord::mdtReal, SensorDims, false, Iy_max_Name);
       MatrixContainer[Iz_sensor_max].SetAllValues(NULL,TMatrixRecord::mdtReal, SensorDims, false, Iz_max_Name);
    }
    
    
    // in case of intensity create buffer for time staggered values
    if ((Params->IsStore_I_avg()) || (Params->IsStore_I_max())){        
       MatrixContainer[p_sensor_i_1_raw] .SetAllValues(NULL,TMatrixRecord::mdtReal, SensorDims, false, "");
       MatrixContainer[ux_sensor_i_1_agr_2].SetAllValues(NULL,TMatrixRecord::mdtReal, SensorDims, false, "");
       MatrixContainer[uy_sensor_i_1_agr_2].SetAllValues(NULL,TMatrixRecord::mdtReal, SensorDims, false, "");
       MatrixContainer[uz_sensor_i_1_agr_2].SetAllValues(NULL,TMatrixRecord::mdtReal, SensorDims, false, "");        
    }
    
    
    
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
 * @param [in] File  File of error
 * @param [in] Line  Line of error
 * 
 */
void TMatrixContainer::PrintErrorAndThrowException(const char * FMT, const string HDF5MatrixName,           
                                          const char * File, const int Line){
    
     fprintf(stderr,FMT, HDF5MatrixName.c_str(), File, Line);                                
     throw bad_alloc();
    
}// end of PrintErrorAndAbort
//------------------------------------------------------------------------------    
    


