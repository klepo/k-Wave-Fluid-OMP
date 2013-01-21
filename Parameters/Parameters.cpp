/**
 * @file        Parameters.cpp
 * @author      Jiri Jaros
 *              CECS, ANU, Australia
 *              jiri.jaros@anu.edu.au
 * @brief       The implementation file containing parameters of the simulation
 * 
 * @version     kspaceFirstOrder3D 2.13
 * @date        9 August 2012, 1:39      (created) \n
 *              14 September 2012, 14:20 (revised)
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

#include <omp.h>
#include <iostream>
#include <string.h>
#include <sstream>
#include <exception>
#include <stdexcept>

#include <Parameters/Parameters.h>

#include <Utils/MatrixNames.h>
#include <Utils/ErrorMessages.h>

using namespace std;

//----------------------------------------------------------------------------//
//                              Constants                                     //
//----------------------------------------------------------------------------//




//----------------------------------------------------------------------------//
//                              Definitions                                   //
//----------------------------------------------------------------------------//





bool TParameters::ParametersInstanceFlag = false;

TParameters* TParameters::ParametersSingleInstance = NULL;


//----------------------------------------------------------------------------//
//                              Implementation                                //
//                              public methods                                //
//----------------------------------------------------------------------------//

/**
 * Get instance of singleton class.
 */
TParameters* TParameters::GetInstance(){

      
  if(!ParametersInstanceFlag)
    {        
        ParametersSingleInstance = new TParameters();
        ParametersInstanceFlag = true;
        return ParametersSingleInstance;
    }
    else
    {
        return ParametersSingleInstance;
    }
    
}// end of Create
//------------------------------------------------------------------------------

/**
 * Parse command line
 * @param [in] argc
 * @param [in] argv
 */
void TParameters::ParseCommandLine(int argc, char** argv){

    CommandLinesParameters.ParseCommandLine(argc, argv);
    
    if (CommandLinesParameters.IsVersion()){
        return;
    }
    
    ReadScalarsFromHDF5InputFile(HDF5_InputFile);
   
    if (CommandLinesParameters.IsBenchmarkFlag())
        Nt = CommandLinesParameters.GetBenchmarkTimeStepsCount();
   
    if ((Nt <= CommandLinesParameters.GetStartTimeIndex()) || 
        ( 0 > CommandLinesParameters.GetStartTimeIndex()) ){
        fprintf(stderr,Parameters_ERR_FMT_Illegal_StartTime_value, (long) 1, Nt);        
        CommandLinesParameters.PrintUsageAndExit();
    }
   
}// end of ParseCommandLine
//------------------------------------------------------------------------------


/**
 * Read scalar values from the input HDF5 file.
 * 
 * @param [in] HDF5_InputFile - Handle to an opened input file
 */
void TParameters::ReadScalarsFromHDF5InputFile(THDF5_File & HDF5_InputFile){
        
        
    TDimensionSizes ScalarSizes(1,1,1);
    
    if (!HDF5_InputFile.IsOpened()) {
    
        // Open file
        try{
          HDF5_InputFile.Open(CommandLinesParameters.GetInputFileName().c_str());      
        } catch (ios::failure e){
            fprintf(stderr,"%s",e.what());
            PrintUsageAndExit();
        }
    
    }    
    
    HDF5_FileHeader.ReadHeaderFromInputFile(HDF5_InputFile);
    
    if (HDF5_FileHeader.GetFileType() != THDF5_FileHeader::hdf5_ft_input) {
       char ErrorMessage[256] = "";
       sprintf(ErrorMessage,Parameters_ERR_FMT_IncorrectInputFileFormat,GetInputFileName().c_str());
       throw ios::failure(ErrorMessage);
    }

    if (!HDF5_FileHeader.CheckMajorFileVersion()) {
       char ErrorMessage[256] = "";
       sprintf(ErrorMessage,Parameters_ERR_FMT_IncorrectMajorHDF5FileVersion,GetInputFileName().c_str(),
               HDF5_FileHeader.GetSupportedHDF5_MajorVersion().c_str());
       throw ios::failure(ErrorMessage);
    }

    if (!HDF5_FileHeader.CheckMinorFileVersion()) {
       char ErrorMessage[256] = "";
       sprintf(ErrorMessage,Parameters_ERR_FMT_IncorrectMinorHDF5FileVersion,GetInputFileName().c_str(),
               HDF5_FileHeader.GetSupportedHDF5_MinorVersion().c_str());
       throw ios::failure(ErrorMessage);
    }
    
    
    
    HDF5_InputFile.ReadCompleteDataset(Nt_Name,  ScalarSizes, &Nt);
    
    HDF5_InputFile.ReadCompleteDataset(dt_Name,  ScalarSizes, &dt);
    HDF5_InputFile.ReadCompleteDataset(dx_Name,  ScalarSizes, &dx);
    HDF5_InputFile.ReadCompleteDataset(dy_Name,  ScalarSizes, &dy);
    HDF5_InputFile.ReadCompleteDataset(dz_Name,  ScalarSizes, &dz);
    
    HDF5_InputFile.ReadCompleteDataset(c_ref_Name,  ScalarSizes, &c_ref);
    
    HDF5_InputFile.ReadCompleteDataset(pml_x_size_Name, ScalarSizes, &pml_x_size);
    HDF5_InputFile.ReadCompleteDataset(pml_y_size_Name, ScalarSizes, &pml_y_size);
    HDF5_InputFile.ReadCompleteDataset(pml_z_size_Name, ScalarSizes, &pml_z_size);
        
    HDF5_InputFile.ReadCompleteDataset(pml_x_alpha_Name, ScalarSizes, &pml_x_alpha);
    HDF5_InputFile.ReadCompleteDataset(pml_y_alpha_Name, ScalarSizes, &pml_y_alpha);
    HDF5_InputFile.ReadCompleteDataset(pml_z_alpha_Name, ScalarSizes, &pml_z_alpha);

    
    
    
    long X, Y, Z;
    HDF5_InputFile.ReadCompleteDataset(Nx_Name  , ScalarSizes, &X);    
    HDF5_InputFile.ReadCompleteDataset(Ny_Name  , ScalarSizes, &Y);
    HDF5_InputFile.ReadCompleteDataset(Nz_Name  , ScalarSizes, &Z);
    
    FullDimensionSizes.X = X;
    FullDimensionSizes.Y = Y;
    FullDimensionSizes.Z = Z;

    ReducedDimensionSizes.X = ((X/2) + 1);
    ReducedDimensionSizes.Y = Y;
    ReducedDimensionSizes.Z = Z; 
    

    sensor_mask_ind_size = HDF5_InputFile.GetDatasetElementCount(sensor_mask_index_Name);
    
    // flags
    HDF5_InputFile.ReadCompleteDataset(ux_source_flag_Name        , ScalarSizes, &ux_source_flag);
    HDF5_InputFile.ReadCompleteDataset(uy_source_flag_Name        , ScalarSizes, &uy_source_flag);
    HDF5_InputFile.ReadCompleteDataset(uz_source_flag_Name        , ScalarSizes, &uz_source_flag);
    HDF5_InputFile.ReadCompleteDataset(transducer_source_flag_Name, ScalarSizes, &transducer_source_flag);

    HDF5_InputFile.ReadCompleteDataset(p_source_flag_Name         , ScalarSizes, &p_source_flag);
    HDF5_InputFile.ReadCompleteDataset(p0_source_flag_Name        , ScalarSizes, &p0_source_flag);
    
    HDF5_InputFile.ReadCompleteDataset(nonuniform_grid_flag_Name  , ScalarSizes, &nonuniform_grid_flag);
    HDF5_InputFile.ReadCompleteDataset(absorbing_flag_Name        , ScalarSizes, &absorbing_flag);
    HDF5_InputFile.ReadCompleteDataset(nonlinear_flag_Name        , ScalarSizes, &nonlinear_flag);    
    
    

    //--- Vector sizes ---//
    if (transducer_source_flag == 0) transducer_source_input_size = 0;
    else {    
       transducer_source_input_size = HDF5_InputFile.GetDatasetElementCount(transducer_source_input_Name);
    }

    if ((transducer_source_flag > 0) || (ux_source_flag > 0) || (uy_source_flag > 0) || (uz_source_flag > 0)){        
        
        u_source_index_size = HDF5_InputFile.GetDatasetElementCount(u_source_index_Name);                       
    }


    //-- uxyz_source_flags --//
    if ((ux_source_flag > 0) || (uy_source_flag > 0) || (uz_source_flag > 0)){

        HDF5_InputFile.ReadCompleteDataset(u_source_many_Name, ScalarSizes, &u_source_many);
        HDF5_InputFile.ReadCompleteDataset(u_source_mode_Name, ScalarSizes, &u_source_mode);
    } else{
        u_source_many = 0;
        u_source_mode = 0;
    }    


    //-- p_source_flag --//
    if (p_source_flag != 0) {
        HDF5_InputFile.ReadCompleteDataset(p_source_many_Name, ScalarSizes, &p_source_many);                
        HDF5_InputFile.ReadCompleteDataset(p_source_mode_Name, ScalarSizes, &p_source_mode);                
        
        p_source_index_size = HDF5_InputFile.GetDatasetElementCount(p_source_index_Name);                                       

    } else{        
        p_source_mode = 0;
        p_source_many = 0;
        p_source_index_size = 0;
    }
    
    
    // absorb flag
    if (absorbing_flag != 0) {
        HDF5_InputFile.ReadCompleteDataset(alpha_power_Name, ScalarSizes, &alpha_power); 
        if (alpha_power == 1.0f) {
            fprintf(stderr,"%s", Parameters_ERR_FMT_Illegal_alpha_power_value);
            PrintUsageAndExit();
        }    
           
        alpha_coeff_scalar_flag = HDF5_InputFile.GetDatasetDimensionSizes(alpha_coeff_Name) == ScalarSizes;
        if (alpha_coeff_scalar_flag) HDF5_InputFile.ReadCompleteDataset(alpha_coeff_Name, ScalarSizes, &alpha_coeff_scalar);                
    }
    
       
    c0_scalar_flag = HDF5_InputFile.GetDatasetDimensionSizes(c0_Name) == ScalarSizes;
    if (c0_scalar_flag) HDF5_InputFile.ReadCompleteDataset(c0_Name, ScalarSizes, &c0_scalar);                
    
    if (nonlinear_flag) {
        BonA_scalar_flag = HDF5_InputFile.GetDatasetDimensionSizes(BonA_Name) == ScalarSizes;
        if (BonA_scalar_flag) HDF5_InputFile.ReadCompleteDataset(BonA_Name, ScalarSizes, &BonA_scalar);                    
    }
        
    rho0_scalar_flag = HDF5_InputFile.GetDatasetDimensionSizes(rho0_Name) == ScalarSizes;
    if ( rho0_scalar_flag) {
        HDF5_InputFile.ReadCompleteDataset(rho0_Name, ScalarSizes, &rho0_scalar);                    
        HDF5_InputFile.ReadCompleteDataset(rho0_sgx_Name, ScalarSizes, &rho0_sgx_scalar);                    
        HDF5_InputFile.ReadCompleteDataset(rho0_sgy_Name, ScalarSizes, &rho0_sgy_scalar);                    
        HDF5_InputFile.ReadCompleteDataset(rho0_sgz_Name, ScalarSizes, &rho0_sgz_scalar);                    
    }
    
    
         
}// end of ReadScalarsFromMatlabInputFile
//------------------------------------------------------------------------------


/** 
 * Save scalars into the output HDF5 file.
 * @param [in] HDF5_OutputFile - Handle to an opened output file where to store
 */
void TParameters::SaveScalarsToHDF5File(THDF5_File & HDF5_OutputFile){
    
      
    // Write dimension sizes
    HDF5_OutputFile.WriteScalarValue(Nx_Name  ,(long) FullDimensionSizes.X);    
    HDF5_OutputFile.WriteScalarValue(Ny_Name  ,(long) FullDimensionSizes.Y);
    HDF5_OutputFile.WriteScalarValue(Nz_Name  ,(long) FullDimensionSizes.Z);
    
    HDF5_OutputFile.WriteScalarValue(Nt_Name  ,(long) Nt);
    
    HDF5_OutputFile.WriteScalarValue(dt_Name  , dt);
    HDF5_OutputFile.WriteScalarValue(dx_Name  , dx);
    HDF5_OutputFile.WriteScalarValue(dy_Name  , dy);
    HDF5_OutputFile.WriteScalarValue(dz_Name  , dz);
    
    HDF5_OutputFile.WriteScalarValue(c_ref_Name, c_ref);
    
    
    HDF5_OutputFile.WriteScalarValue(pml_x_size_Name  , pml_x_size);
    HDF5_OutputFile.WriteScalarValue(pml_y_size_Name  , pml_y_size);
    HDF5_OutputFile.WriteScalarValue(pml_z_size_Name  , pml_z_size);
    
    HDF5_OutputFile.WriteScalarValue(pml_x_alpha_Name, pml_x_alpha);
    HDF5_OutputFile.WriteScalarValue(pml_y_alpha_Name, pml_y_alpha);
    HDF5_OutputFile.WriteScalarValue(pml_z_alpha_Name, pml_z_alpha);

    
    
    HDF5_OutputFile.WriteScalarValue(ux_source_flag_Name   , ux_source_flag);
    HDF5_OutputFile.WriteScalarValue(uy_source_flag_Name   , uy_source_flag);
    HDF5_OutputFile.WriteScalarValue(uz_source_flag_Name   , uz_source_flag);
    HDF5_OutputFile.WriteScalarValue(transducer_source_flag_Name   , transducer_source_flag);
    
    HDF5_OutputFile.WriteScalarValue(p_source_flag_Name   , p_source_flag);
    HDF5_OutputFile.WriteScalarValue(p0_source_flag_Name  , p0_source_flag);
    
    HDF5_OutputFile.WriteScalarValue(nonuniform_grid_flag_Name, nonuniform_grid_flag);
    HDF5_OutputFile.WriteScalarValue(absorbing_flag_Name      , absorbing_flag);
    HDF5_OutputFile.WriteScalarValue(nonlinear_flag_Name      , nonlinear_flag);    
    

    //-- uxyz_source_flags --//
    if ((ux_source_flag > 0) || (uy_source_flag > 0) || (uz_source_flag > 0)){

        HDF5_OutputFile.WriteScalarValue(u_source_many_Name   , u_source_many);
        HDF5_OutputFile.WriteScalarValue(u_source_mode_Name   , u_source_mode);                
    }

    //-- p_source_flag --//
    if (p_source_flag != 0) {
        
        HDF5_OutputFile.WriteScalarValue(p_source_many_Name   , p_source_many);
        HDF5_OutputFile.WriteScalarValue(p_source_mode_Name   , p_source_mode);                
        
    }
        
    
    // absorb flag
    if (absorbing_flag != 0) {
        HDF5_OutputFile.WriteScalarValue(alpha_power_Name, alpha_power);                
        
    }

    
    
}// end of SaveScalarsToHDF5File




//----------------------------------------------------------------------------//
//                              Implementation                                //
//                            protected methods                               //
//----------------------------------------------------------------------------//


/**
 * Constructor 
 */
TParameters::TParameters() : 
        HDF5_InputFile(), HDF5_OutputFile(), HDF5_FileHeader(),
        CommandLinesParameters(),
        Nt(0),dt(0.0f),    
        dx(0.0f), dy(0.0f), dz(0.0f), 
        c_ref(0.0f), alpha_power(0.0f),
        FullDimensionSizes(0,0,0), ReducedDimensionSizes(0,0,0),
        sensor_mask_ind_size (0), u_source_index_size(0), p_source_index_size(0), transducer_source_input_size(0),        
        ux_source_flag(0), uy_source_flag(0), uz_source_flag(0),
        p_source_flag(0), p0_source_flag(0), transducer_source_flag(0),                  
        u_source_many(0), u_source_mode(0), p_source_mode(0), p_source_many(0),
        nonuniform_grid_flag(0), absorbing_flag(0), nonlinear_flag(0),    
        pml_x_size(0), pml_y_size(0), pml_z_size(0),    
        alpha_coeff_scalar_flag(false), alpha_coeff_scalar(0.0f),
        c0_scalar_flag(false), c0_scalar(0.0f),
        absorb_eta_scalar(0.0f), absorb_tau_scalar (0.0f),
        BonA_scalar_flag(false), BonA_scalar (0.0f),
        rho0_scalar_flag(false), rho0_scalar(0.0f), rho0_sgx_scalar(0.0f), rho0_sgy_scalar(0.0f), rho0_sgz_scalar(0.0f)
        
{                         
            
}// end of TFFT1DParameters
//------------------------------------------------------------------------------



//----------------------------------------------------------------------------//
//                              Implementation                                //
//                              private methods                               //
//----------------------------------------------------------------------------//

/**
 * print usage end exit
 */
void TParameters::PrintUsageAndExit(){
    
 CommandLinesParameters.PrintUsageAndExit();        
  
    
}// end of PrintUsage
//------------------------------------------------------------------------------

