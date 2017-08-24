/**
 * @file        Parameters.cpp
 * @author      Jiri Jaros
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The implementation file containing parameters of the simulation.
 *
 * @version     kspaceFirstOrder3D 2.16
 *
 * @date        09 August    2012, 13:39 (created) \n
 *              24 August    2017, 12:23 (revised)
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

#ifdef _OPENMP
  #include <omp.h>
#endif

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
TParameters* TParameters::GetInstance()
{
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
}// end of GetInstance
//------------------------------------------------------------------------------

/**
 * Parse command line.
 * @param [in] argc
 * @param [in] argv
 */
void TParameters::ParseCommandLine(int argc, char** argv)
{
  CommandLinesParameters.ParseCommandLine(argc, argv);

  if (CommandLinesParameters.IsVersion())
  {
    return;
  }

  ReadScalarsFromHDF5InputFile(HDF5_InputFile);

  if (CommandLinesParameters.IsBenchmarkFlag())
  {
    Nt = CommandLinesParameters.GetBenchmarkTimeStepsCount();
  }

  if ((Nt <= (size_t) CommandLinesParameters.GetStartTimeIndex()) ||
      ( 0 > CommandLinesParameters.GetStartTimeIndex()) )
  {
    fprintf(stderr,kErrFmtIllegalSamplingStartTimeStep, (size_t) 1, Nt);
    CommandLinesParameters.PrintUsageAndExit();
  }
}// end of ParseCommandLine
//------------------------------------------------------------------------------


/**
 * Read scalar values from the input HDF5 file.
 *
 * @param [in] HDF5_InputFile - Handle to an opened input file.
 * @throw ios:failure if the file cannot be open or is of a wrong type or version.
 */
void TParameters::ReadScalarsFromHDF5InputFile(Hdf5File & HDF5_InputFile)
{
  DimensionSizes ScalarSizes(1, 1, 1);

  if (!HDF5_InputFile.IsOpened())
  {
    // Open file
    try
    {
      HDF5_InputFile.Open(CommandLinesParameters.GetInputFileName().c_str());
    }
    catch (ios::failure e)
    {
      fprintf(stderr, "%s", e.what());
      PrintUsageAndExit();
    }
  }

  HDF5_FileHeader.readHeaderFromInputFile(HDF5_InputFile);

  // check file type
  if (HDF5_FileHeader.getFileType() != Hdf5FileHeader::FileType::kInput)
  {
    char ErrorMessage[256] = "";
    sprintf(ErrorMessage, kErrFmtBadInputFileFormat, GetInputFileName().c_str());
    throw ios::failure(ErrorMessage);
  }

  // check version
  if (!HDF5_FileHeader.checkMajorFileVersion())
  {
    char ErrorMessage[256] = "";
    sprintf(ErrorMessage, kErrFmtBadMajorFileVersion, GetInputFileName().c_str(),
            HDF5_FileHeader.getFileMajorVersion().c_str());
    throw ios::failure(ErrorMessage);
  }

  if (!HDF5_FileHeader.checkMinorFileVersion())
  {
    char ErrorMessage[256] = "";
    sprintf(ErrorMessage, kErrFmtBadMinorFileVersion, GetInputFileName().c_str(),
            HDF5_FileHeader.getFileMinorVersion().c_str());
    throw ios::failure(ErrorMessage);
  }

  const hid_t HDF5RootGroup = HDF5_InputFile.GetRootGroup();

  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kNtName, Nt);

  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kDtName, dt);
  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kDxName, dx);
  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kDyName, dy);
  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kDzName, dz);

  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kCRefName, c_ref);

  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kPmlXSizeName, pml_x_size);
  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kPmlYSizeName, pml_y_size);
  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kPmlZSizeName, pml_z_size);

  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kPmlXAlphaName, pml_x_alpha);
  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kPmlYAlphaName, pml_y_alpha);
  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kPmlZAlphaName, pml_z_alpha);

  size_t X, Y, Z;
  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kNxName, X);
  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kNyName, Y);
  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kNzName, Z);

  FullDimensionSizes.nx = X;
  FullDimensionSizes.ny = Y;
  FullDimensionSizes.nz = Z;

  ReducedDimensionSizes.nx = ((X / 2) + 1);
  ReducedDimensionSizes.ny = Y;
  ReducedDimensionSizes.nz = Z;

  // if the file is of version 1.0, there must be a sensor mask index (backward compatibility)
  if (HDF5_FileHeader.getFileVersion() == Hdf5FileHeader::FileVersion::kVersion10)
  {
    sensor_mask_ind_size = HDF5_InputFile.GetDatasetElementCount(HDF5RootGroup, kSensorMaskIndexName);

    //if -u_non_staggered_raw enabled, throw an error - not supported
    if (IsStore_u_non_staggered_raw())
    {
      throw ios::failure(kErrFmtNonStaggeredVelocityNotSupportedFileVersion);
    }
  }

  // This is the current version 1.1
  if (HDF5_FileHeader.getFileVersion() == Hdf5FileHeader::FileVersion::kVersion11)
  {

    // read sensor mask type as a size_t value to enum
    size_t SensorMaskTypeNumericalue = 0;
    HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kSensorMaskTypeName, SensorMaskTypeNumericalue);

    // convert the size_t value to enum
    switch (SensorMaskTypeNumericalue)
    {
      case 0: sensor_mask_type = smt_index;
        break;
      case 1: sensor_mask_type = smt_corners;
        break;
      default:
      {
        throw ios::failure(kErrFmtBadSensorMaskType);
        break;
      }
    }//case

    // read the input mask size
    switch (sensor_mask_type)
    {
      case smt_index:
      {
        sensor_mask_ind_size = HDF5_InputFile.GetDatasetElementCount(HDF5RootGroup, kSensorMaskIndexName);
        break;
      }
      case smt_corners:
      {
        // mask dimensions are [6, N, 1] - I want to know N
        sensor_mask_corners_size = HDF5_InputFile.GetDatasetDimensionSizes(HDF5RootGroup, kSensorMaskCornersName).ny;
        break;
      }
    }// switch
  }// version 1.1


  // flags.
  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kVelocityXSourceFlagName, ux_source_flag);
  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kVelocityYSourceFlagName, uy_source_flag);
  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kVelocityZSourceFlagName, uz_source_flag);
  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kTransducerSourceFlagName, transducer_source_flag);

  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kPressureSourceFlagName, p_source_flag);
  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kInitialPressureSourceFlagName,p0_source_flag);

  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kNonUniformGridFlagName, nonuniform_grid_flag);
  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kAbsorbingFlagName, absorbing_flag);
  HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kNonLinearFlagName, nonlinear_flag);



  // Vector sizes.
  if (transducer_source_flag == 0)
  {
   transducer_source_input_size = 0;
  }
  else
  {
    transducer_source_input_size = HDF5_InputFile.GetDatasetElementCount(HDF5RootGroup, kInitialPressureSourceInputName);
  }

  if ((transducer_source_flag > 0) || (ux_source_flag > 0) || (uy_source_flag > 0) || (uz_source_flag > 0))
  {
    u_source_index_size = HDF5_InputFile.GetDatasetElementCount(HDF5RootGroup, kVelocitySourceIndexName);
  }


  // uxyz_source_flags
  if ((ux_source_flag > 0) || (uy_source_flag > 0) || (uz_source_flag > 0))
  {
    HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kVelocitySourceManyName, u_source_many);
    HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kVelocitySourceModeName, u_source_mode);
  }
  else
  {
    u_source_many = 0;
    u_source_mode = 0;
  }

  // p_source_flag
  if (p_source_flag != 0)
  {
    HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kPressureSourceManyName, p_source_many);
    HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kPressureSourceModeName, p_source_mode);

    p_source_index_size = HDF5_InputFile.GetDatasetElementCount(HDF5RootGroup, kPressureSourceIndexName);
  }
  else
  {
    p_source_mode = 0;
    p_source_many = 0;
    p_source_index_size = 0;
  }


  // absorb flag
  if (absorbing_flag != 0)
  {
    HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kAlphaPowerName, alpha_power);
    if (alpha_power == 1.0f)
    {
      fprintf(stderr, "%s", kErrFmtIllegalAlphaPowerValue);
      PrintUsageAndExit();
    }

    alpha_coeff_scalar_flag = HDF5_InputFile.GetDatasetDimensionSizes(HDF5RootGroup, kAlphaCoeffName) == ScalarSizes;
    if (alpha_coeff_scalar_flag)
    {
      HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kAlphaCoeffName, alpha_coeff_scalar);
    }
  }


  c0_scalar_flag = HDF5_InputFile.GetDatasetDimensionSizes(HDF5RootGroup, kC0Name) == ScalarSizes;
  if (c0_scalar_flag)
  {
    HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kC0Name, c0_scalar);
  }

  if (nonlinear_flag)
  {
    BonA_scalar_flag = HDF5_InputFile.GetDatasetDimensionSizes(HDF5RootGroup, kBonAName) == ScalarSizes;
    if (BonA_scalar_flag)
    {
      HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kBonAName, BonA_scalar);
    }
  }

  rho0_scalar_flag = HDF5_InputFile.GetDatasetDimensionSizes(HDF5RootGroup, kRho0Name) == ScalarSizes;
  if (rho0_scalar_flag)
  {
    HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kRho0Name, rho0_scalar);
    HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kRho0SgxName, rho0_sgx_scalar);
    HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kRho0SgyName, rho0_sgy_scalar);
    HDF5_InputFile.ReadScalarValue(HDF5RootGroup, kRho0SgzName, rho0_sgz_scalar);
  }
}// end of ReadScalarsFromHDF5InputFile
//------------------------------------------------------------------------------

/**
 * Save scalars into the output HDF5 file.
 * @param [in] HDF5_OutputFile - Handle to an opened output file where to store
 */
void TParameters::SaveScalarsToHDF5File(Hdf5File & HDF5_OutputFile)
{
  const hid_t HDF5RootGroup = HDF5_OutputFile.GetRootGroup();

  // Write dimension sizes
  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kNxName, FullDimensionSizes.nx);
  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kNyName, FullDimensionSizes.ny);
  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kNzName, FullDimensionSizes.nz);

  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kNtName, Nt);

  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kDtName, dt);
  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kDxName, dx);
  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kDyName, dy);
  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kDzName, dz);

  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kCRefName, c_ref);

  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kPmlXSizeName, pml_x_size);
  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kPmlYSizeName, pml_y_size);
  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kPmlZSizeName, pml_z_size);

  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kPmlXAlphaName, pml_x_alpha);
  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kPmlYAlphaName, pml_y_alpha);
  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kPmlZAlphaName, pml_z_alpha);

  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kVelocityXSourceFlagName, ux_source_flag);
  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kVelocityYSourceFlagName, uy_source_flag);
  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kVelocityZSourceFlagName, uz_source_flag);
  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kTransducerSourceFlagName, transducer_source_flag);

  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kPressureSourceFlagName, p_source_flag);
  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kInitialPressureSourceFlagName, p0_source_flag);

  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kNonUniformGridFlagName, nonuniform_grid_flag);
  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kAbsorbingFlagName, absorbing_flag);
  HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kNonLinearFlagName, nonlinear_flag);


  // uxyz_source_flags
  if ((ux_source_flag > 0) || (uy_source_flag > 0) || (uz_source_flag > 0))
  {
    HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kVelocitySourceManyName, u_source_many);
    HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kVelocitySourceModeName, u_source_mode);
  }

  // p_source_flag
  if (p_source_flag != 0)
  {
    HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kPressureSourceManyName, p_source_many);
    HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kPressureSourceModeName, p_source_mode);
  }

  // absorb flag
  if (absorbing_flag != 0)
  {
    HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kAlphaPowerName, alpha_power);
  }

  // if copy sensor mask, then copy the mask type
  if (IsCopySensorMask())
  {
    size_t SensorMaskTypeNumericValue = 0;

    switch (sensor_mask_type)
    {
      case smt_index: SensorMaskTypeNumericValue = 0;
        break;
      case smt_corners: SensorMaskTypeNumericValue = 1;
        break;
    }// switch

    HDF5_OutputFile.WriteScalarValue(HDF5RootGroup, kSensorMaskTypeName, SensorMaskTypeNumericValue);
  }
}// end of SaveScalarsToHDF5File
//------------------------------------------------------------------------------



//----------------------------------------------------------------------------//
//                              Implementation                                //
//                            protected methods                               //
//----------------------------------------------------------------------------//


/**
 * Constructor.
 */
TParameters::TParameters() :
        HDF5_InputFile(), HDF5_OutputFile(), HDF5_CheckpointFile(), HDF5_FileHeader(),
        CommandLinesParameters(),
        Nt(0), t_index(0), dt(0.0f),
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
 * Print usage end exit.
 */
void TParameters::PrintUsageAndExit()
{
  CommandLinesParameters.PrintUsageAndExit();
}// end of PrintUsage
//------------------------------------------------------------------------------

