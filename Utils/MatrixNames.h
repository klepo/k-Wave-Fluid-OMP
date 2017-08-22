/**
 * @file        MatrixNames.h
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The header file storing names of all variables.
 *
 * @version     kspaceFirstOrder3D 2.16
 * @date        14 September 2012, 14:33 (created) \n
 *              22 August    2017, 14:33 (revised)
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

#ifndef MATRIX_NAMES_H
#define MATRIX_NAMES_H

/**
 * @brief   Datatype for matrix names.
 * @details Datatype for matrix names.
 */
using MatrixName = const char*;


//--------------------------------------------------------------------------------------------------------------------//
//--------------------------------------------------- Constants ------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/// Nt variable name
MatrixName const kNtName         = "Nt";
/// t_index name
MatrixName const kTimeIndexName  = "t_index";
/// dt variable name
MatrixName const kDtName         = "dt";
/// dx variable name
MatrixName const kDxName         = "dx";
/// dy variable name
MatrixName const kDyName         = "dy";
/// dz variable name
MatrixName const kDzName         = "dz";

/// c_ref variable name
MatrixName const kCRefName       = "c_ref";
/// c0 variable name
MatrixName const kC0Name         = "c0";

/// alpha_power variable name
MatrixName const kAlphaPowerName = "alpha_power";
/// alpha_coeff variable name
MatrixName const kAlphaCoeffName = "alpha_coeff";

/// Nx variable name
MatrixName const kNxName         = "Nx";
/// Ny variable name
MatrixName const kNyName         = "Ny";
/// Nz variable name
MatrixName const kNzName         = "Nz";

/// x_shift_neg_r variable name
MatrixName const kXShiftNegRName = "x_shift_neg_r";
/// y_shift_neg_r variable name
MatrixName const kYShiftNegRName = "y_shift_neg_r";
/// z_shift_neg_r variable name
MatrixName const kZShiftNegRName = "z_shift_neg_r";

/// ux_shifted variable name
MatrixName const kUxShiftedName  = "ux_shifted";
/// uy_shifted variable name
MatrixName const kUyShiftedName  = "uy_shifted";
/// uz_shifted variable name
MatrixName const kUzShiftedName  = "uz_shifted";

/// pml_x_size variable name
MatrixName const kPmlXSizeName   = "pml_x_size";
/// pml_y_size variable name
MatrixName const kPmlYSizeName   = "pml_y_size";
/// pml_z_size variable name
MatrixName const kPmlZSizeName   = "pml_z_size";

/// pml_x_sgx variable name
MatrixName const kPmlXSgxName    = "pml_x_sgx";
/// pml_y_sgy variable name
MatrixName const kPmlYSgyName    = "pml_y_sgy";
/// pml_z_sgz variable name
MatrixName const kPmlZSgzName    = "pml_z_sgz";

// pml_x variable name
MatrixName const kPmlXName       = "pml_x";
/// pml_y variable name
MatrixName const kPmlYName       = "pml_y";
/// pml_z variable name
MatrixName const kPmlZName       = "pml_z";


/// pml_x_alpha variable name
MatrixName const kPmlXAlphaName    = "pml_x_alpha";
/// pml_y_alpha variable name
MatrixName const kPmlYAlphaName    = "pml_y_alpha";
/// pml_z_alpha variable name
MatrixName const kPmlZAlphaName    = "pml_z_alpha";

/// ux_source_flag variable name
MatrixName const kVelocityXSourceFlagName = "ux_source_flag";
/// uy_source_flag variable name
MatrixName const kVelocityYSourceFlagName = "uy_source_flag";
/// uz_source_flag variable name
MatrixName const kVelocityZSourceFlagName = "uz_source_flag";

/// u_source_many variable name
MatrixName const kVelocitySourceManyName  = "u_source_many";
/// p_source_many variable name
MatrixName const kPressureSourceManyName  = "p_source_many";

/// p_source_flag variable name
MatrixName const kPressureSourceFlagName        = "p_source_flag";
/// p0_source_flag variable name
MatrixName const kInitialPressureSourceFlagName = "p0_source_flag";

/// u_source_mode variable name
MatrixName const kVelocitySourceModeName  = "u_source_mode";
/// p_source_mode variable name
MatrixName const kPressureSourceModeName  = "p_source_mode";

/// p_source_input variable name
MatrixName const kPressureSourceInputName = "p_source_input";
/// p_source_index variable name
MatrixName const kPressureSourceIndexName = "p_source_index";

/// u_source_index variable name
MatrixName const kVelocitySourceIndexName  = "u_source_index";
/// ux_source_input variable name
MatrixName const kVelocityXSourceInputName = "ux_source_input";
/// uy_source_input variable name
MatrixName const kVelocityYSourceInputName = "uy_source_input";
/// uz_source_input variable name
MatrixName const kVelocityZSourceInputName = "uz_source_input";

/// nonuniform_grid_flag variable name
MatrixName const kNonUniformGridFlagName   = "nonuniform_grid_flag";
/// absorbing_flag variable name
MatrixName const kAbsorbingFlagName        = "absorbing_flag";
/// nonlinear_flag variable name
MatrixName const kNonLinearFlagName        = "nonlinear_flag";

/// transducer_source_flag variable name
MatrixName const kTransducerSourceFlagName = "transducer_source_flag";
/// sensor_mask_index variable name
MatrixName const kSensorMaskIndexName      = "sensor_mask_index";
/// sensor_mask_type variable name
MatrixName const kSensorMaskTypeName       = "sensor_mask_type";
/// sensor_mask_corners variable name
MatrixName const kSensorMaskCornersName    = "sensor_mask_corners";

/// transducer_source_input variable name
MatrixName const kTransducerSourceInputName = "transducer_source_input";

/// p0_source_input variable name
MatrixName const kInitialPressureSourceInputName = "p0_source_input";
/// delay_mask variable name
MatrixName const kDelayMaskName                  = "delay_mask";


/// kappa_r variable name
MatrixName const kKappaRName = "kappa_r";
/// BonA variable name
MatrixName const kBonAName   = "BonA";
/// p variable name
MatrixName const kPName      = "p";
/// rhox variable name
MatrixName const kRhoXName   = "rhox";
/// rhoy variable name
MatrixName const kRhoYName   = "rhoy";
/// rhoz variable name
MatrixName const kRhoZName   = "rhoz";

/// ux variable name
MatrixName const kUxName     = "ux";
/// uy variable name
MatrixName const kUyName     = "uy";
/// uz variable name
MatrixName const kUzName     = "uz";

/// ux_sgx variable name
MatrixName const kUxSgxName  = "ux_sgx";
/// uy_sgy variable name
MatrixName const kUySgyName  = "uy_sgy";
/// uz_sgz variable name
MatrixName const kUzSgzName  = "uz_sgz";

/// ux_non_staggered variable name
MatrixName const kUxNonStaggeredName = "ux_non_staggered";
/// uy_non_staggered variable name
MatrixName const kUyNonStaggeredName = "uy_non_staggered";
/// uz_non_staggered variable name
MatrixName const kUzNonStaggeredName = "uz_non_staggered";

/// duxdx variable name
MatrixName const kDuxdxName          = "duxdx";
/// duydy variable name
MatrixName const kDuydyName          = "duydy";
/// duzdz variable name
MatrixName const kDuzdzName          = "duzdz";

/// dxudxn variable name
MatrixName const kDxudxnName         = "dxudxn";
/// dyudyn variable name
MatrixName const kDyudynName         = "dyudyn";
/// dzudzn variable name
MatrixName const kDzudznName         = "dzudzn";

/// dxudxn_sgx variable name
MatrixName const kDxudxnSgxName      = "dxudxn_sgx";
/// dyudyn_sgy variable name
MatrixName const kDyudynSgyName      = "dyudyn_sgy";
/// dzudzn_sgz variable name
MatrixName const kDzudznSgzName      = "dzudzn_sgz";

/// ddx_k_shift_pos_r variable name
MatrixName const kDdxKShiftPosRName  = "ddx_k_shift_pos_r";
/// ddy_k_shift_pos variable name
MatrixName const kDdyKShiftPosName   = "ddy_k_shift_pos";
/// ddz_k_shift_pos variable name
MatrixName const kDdzKShiftPosName   = "ddz_k_shift_pos";

/// ddx_k_shift_neg_r variable name
MatrixName const kDdxKShiftNegRName  = "ddx_k_shift_neg_r";
/// ddy_k_shift_neg variable name
MatrixName const kDdyKShiftNegName   = "ddy_k_shift_neg";
/// ddz_k_shift_neg variable name
MatrixName const kDdzKShiftNegName   = "ddz_k_shift_neg";

/// rho0 variable name
MatrixName const kRho0Name           = "rho0";
/// rho0_sgx variable name
MatrixName const kRho0SgxName        = "rho0_sgx";
/// rho0_sgy variable name
MatrixName const kRho0SgyName        = "rho0_sgy";
/// rho0_sgz variable name
MatrixName const kRho0SgzName        = "rho0_sgz";

/// absorb_tau variable name
MatrixName const kAbsorbTauName      = "absorb_tau";
/// absorb_eta variable name
MatrixName const kAbsorbEtaName      = "absorb_eta";
/// absorb_nabla1_r variable name
MatrixName const kAbsorbNabla1RName  = "absorb_nabla1_r";
/// absorb_nabla2_r variable name
MatrixName const kAbsorbNabla2RName  = "absorb_nabla2_r";

/// p variable name in the output file
MatrixName const kPressureRawName    = "p";
/// p_rms variable name
MatrixName const kPressureRmsName    = "p_rms";
/// p_max variable name
MatrixName const kPressureMaxName    = "p_max";
/// p_min variable name
MatrixName const kPressureMinName    = "p_min";
/// p_max_all variable name
MatrixName const kPressureMaxAllName = "p_max_all";
/// p_min_all variable name
MatrixName const kPressureMinAllName = "p_min_all";
/// p_final variable name
MatrixName const kPressureFinalName  = "p_final";

/// ux_rms variable name
MatrixName const kUxRmsName = "ux_rms";
/// uy_rms variable name
MatrixName const kUyRmsName = "uy_rms";
/// uz_rms variable name
MatrixName const kUzRmsName = "uz_rms";

/// ux_max variable name
MatrixName const kUxMaxName = "ux_max";
/// uy_max variable name
MatrixName const kUyMaxName = "uy_max";
/// uz_max variable name
MatrixName const kUzMaxName = "uz_max";
/// ux_min variable name
MatrixName const kUxMinName = "ux_min";
/// uy_min variable name
MatrixName const kUyMinName = "uy_min";
/// uz_min variable name
MatrixName const kUzMinName = "uz_min";

/// ux_max_all variable name
MatrixName const kUxMaxAllName = "ux_max_all";
/// uy_max_all variable name
MatrixName const kUyMaxAllName = "uy_max_all";
/// uz_max_all variable name
MatrixName const kUzMaxAllName = "uz_max_all";
/// ux_min_all variable name
MatrixName const kUxMinAllName = "ux_min_all";
/// uy_min_all variable name
MatrixName const kUyMinAllName = "uy_min_all";
/// uz_min_all variable name
MatrixName const kUzMinAllName = "uz_min_all";

/// ux_final variable name
MatrixName const kUxFinalName = "ux_final";
/// uy_final variable name
MatrixName const kUyFinalName = "uy_final";
/// uz_final variable name
MatrixName const kUzFinalName = "uz_final";

/// Temp_1_RS3D variable name
MatrixName const kTemp1Real3DName = "Temp_1_RS3D";
/// Temp_2_RS3D variable name
MatrixName const kTemp2Real3DName = "Temp_2_RS3D";
/// Temp_3_RS3D variable name
MatrixName const kTemp3Real3DName = "Temp_3_RS3D";


/// CUFFT_shift_temp variable name
MatrixName const kCufftShiftTempName = "CUFFT_shift_temp";
/// CUFFT_X_temp variable name
MatrixName const kCufftXTempName     = "CUFFT_X_temp";
/// CUFFT_Y_temp variable name
MatrixName const kCufftYTempName     = "CUFFT_Y_temp";
/// CUFFT_Z_temp variable name
MatrixName const kCufftZTempName     = "CUFFT_Z_temp";

#endif	/* MATRIX_NAMES_H */
