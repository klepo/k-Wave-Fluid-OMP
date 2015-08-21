/**
 * @file        KSpaceFirstOrder3DSolver.cpp
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The implementation file containing the main class of the project
 *              responsible for the entire simulation.
 *
 * @version     kspaceFirstOrder3D 2.16
 * @date        12 July      2012, 10:27  (created)\n
 *              21 August    2015, 18:28  (revised)
 *
 * @section License
 * This file is part of the C++ extension of the k-Wave Toolbox (http://www.k-wave.org).\n
 * Copyright (C) 2014 Jiri Jaros and Bradley Treeby.
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

// Linux build
#ifdef __linux__
  #include <sys/resource.h>
  #include <cmath>
#endif

// Windows build
#ifdef _WIN64
  #define _USE_MATH_DEFINES
  #include <cmath>
  #include <Windows.h>
  #include <Psapi.h>
  #pragma comment(lib, "Psapi.lib")
#endif

#ifdef _OPENMP
  #include <omp.h>
#endif

#include <iostream>
#include <sstream>
#include <cstdio>
#include <limits>

#include <immintrin.h>
#include <time.h>

#include <KSpaceSolver/KSpaceFirstOrder3DSolver.h>

#include <Utils/ErrorMessages.h>

#include <MatrixClasses/FFTWComplexMatrix.h>
#include <MatrixClasses/MatrixContainer.h>

using namespace std;

//----------------------------------------------------------------------------//
//                              Constants                                     //
//----------------------------------------------------------------------------//




//----------------------------------------------------------------------------//
//                              Public methods                                //
//----------------------------------------------------------------------------//


/**
 * Constructor of the class.
 */
TKSpaceFirstOrder3DSolver::TKSpaceFirstOrder3DSolver():
        MatrixContainer(), OutputStreamContainer(),
        ActPercent(0),
        Parameters(NULL),
        TotalTime(), PreProcessingTime(), DataLoadTime (), SimulationTime(),
        PostProcessingTime(), IterationTime()
{
  TotalTime.Start();
  Parameters = TParameters::GetInstance();

  //Switch off HDF5 error messages
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);
}// end of TKSpaceFirstOrder3DSolver
//------------------------------------------------------------------------------


/**
 * Destructor of the class.
 */
TKSpaceFirstOrder3DSolver::~TKSpaceFirstOrder3DSolver()
{
  FreeMemory();
}// end of TKSpace3DSolver
//------------------------------------------------------------------------------

/**
 * The method allocates the matrix container, creates all matrices and
 * creates all output streams (however not allocating memory).
 */
void TKSpaceFirstOrder3DSolver::AllocateMemory()
{
  // create container, then all matrices
  MatrixContainer.AddMatricesIntoContainer();
  MatrixContainer.CreateAllObjects();

  // add output streams into container
  //@todo Think about moving under LoadInputData routine...
  OutputStreamContainer.AddStreamsIntoContainer(MatrixContainer);
}// end of AllocateMemory
//------------------------------------------------------------------------------

/**
 * The method frees all memory allocated by the class.
 */
void TKSpaceFirstOrder3DSolver::FreeMemory()
{
  MatrixContainer.FreeAllMatrices();
  OutputStreamContainer.FreeAllStreams();
}// end of FreeMemory
//------------------------------------------------------------------------------

/**
 * Load data from the input file provided by the Parameter class and creates
 * the output time series streams.
 */
void TKSpaceFirstOrder3DSolver::LoadInputData()
{
  // Start timer
  DataLoadTime.Start();

  // get handles
  THDF5_File& HDF5_InputFile      = Parameters->HDF5_InputFile; // file is opened (in Parameters)
  THDF5_File& HDF5_OutputFile     = Parameters->HDF5_OutputFile;
  THDF5_File& HDF5_CheckpointFile = Parameters->HDF5_CheckpointFile;

  // Load data from disk
  MatrixContainer.LoadDataFromInputHDF5File(HDF5_InputFile);

  // close the input file
  HDF5_InputFile.Close();

  // The simulation does not use checkpointing or this is the first turn
  bool RecoverFromPrevState = (Parameters->IsCheckpointEnabled() &&
                               THDF5_File::IsAccessible(Parameters->GetCheckpointFileName().c_str()));

  //-------------------- Read data from the checkpoint file -----------------//
  if (RecoverFromPrevState)
  {
    // Open checkpoint file
    HDF5_CheckpointFile.Open(Parameters->GetCheckpointFileName().c_str());

    // Check the checkpoint file
    CheckCheckpointFile();

    // read the actual value of t_index
    size_t new_t_index;
    HDF5_CheckpointFile.ReadScalarValue(HDF5_CheckpointFile.GetRootGroup(),
                                        t_index_Name,
                                        new_t_index);
    Parameters->Set_t_index(new_t_index);

    // Read necessary matrices from the checkpoint file
    MatrixContainer.LoadDataFromCheckpointHDF5File(HDF5_CheckpointFile);

    HDF5_CheckpointFile.Close();


    //------------- Read data from the output file -------------------------//

    // Reopen output file for RW access
    HDF5_OutputFile.Open(Parameters->GetOutputFileName().c_str(), H5F_ACC_RDWR);
    //Read file header of the output file
    Parameters->HDF5_FileHeader.ReadHeaderFromOutputFile(HDF5_OutputFile);
    // Check the checkpoint file
    CheckOutputFile();
    // Restore elapsed time
    RestoreCumulatedElapsedFromOutputFile();

    OutputStreamContainer.ReopenStreams();
  }
  else
  { //-------------------- First round of multi-leg simulation ---------------//
    // Create the output file
    HDF5_OutputFile.Create(Parameters->GetOutputFileName().c_str());


    // Create the steams, link them with the sampled matrices, however DO NOT allocate memory!
    OutputStreamContainer.CreateStreams();

  }

 // Stop timer
  DataLoadTime.Stop();
}// end of LoadInputData
//------------------------------------------------------------------------------


/**
 * This method computes k-space First Order 3D simulation.
 * It launches calculation on a given dataset going through
 * FFT initialization, pre-processing, main loop and post-processing phases.
 *
 */
void TKSpaceFirstOrder3DSolver::Compute()
{
  PreProcessingTime.Start();

  fprintf(stdout,"FFT plans creation.........."); fflush(stdout);

    InitializeFFTWPlans( );

  fprintf(stdout,"Done \n");
  fprintf(stdout,"Pre-processing phase........"); fflush(stdout);


    PreProcessingPhase( );

  fprintf(stdout,"Done \n");
  fprintf(stdout,"Current memory in use:%8ldMB\n", ShowMemoryUsageInMB());
  PreProcessingTime.Stop();

  fprintf(stdout,"Elapsed time:          %8.2fs\n",PreProcessingTime.GetElapsedTime());

  SimulationTime.Start();

    Compute_MainLoop();

  SimulationTime.Stop();

  PostProcessingTime.Start();

  if (IsCheckpointInterruption())
  { // Checkpoint
    fprintf(stdout,"-------------------------------------------------------------\n");
    fprintf(stdout,".............. Interrupted to checkpoint! ...................\n");
    fprintf(stdout,"Number of time steps completed:                    %10ld\n",  Parameters->Get_t_index());
    fprintf(stdout,"Elapsed time:                                       %8.2fs\n",SimulationTime.GetElapsedTime());
    fprintf(stdout,"-------------------------------------------------------------\n");
    fprintf(stdout,"Checkpoint in progress......"); fflush(stdout);

    SaveCheckpointData();
  }
  else
  { // Finish
    fprintf(stdout,"-------------------------------------------------------------\n");
    fprintf(stdout,"Elapsed time:                                       %8.2fs\n",SimulationTime.GetElapsedTime());
    fprintf(stdout,"-------------------------------------------------------------\n");
    fprintf(stdout,"Post-processing phase......."); fflush(stdout);

    PostPorcessing();

    // if checkpointing is enabled and the checkpoint file was created created in the past, delete it
    if (Parameters->IsCheckpointEnabled())
    {
      std::remove(Parameters->GetCheckpointFileName().c_str());
      TFFTWComplexMatrix::DeleteStoredWisdom();
    }
  }

  PostProcessingTime.Stop();

  fprintf(stdout,"Done \n");
  fprintf(stdout,"Elapsed time:          %8.2fs\n",PostProcessingTime.GetElapsedTime());

    WriteOutputDataInfo();

  Parameters->HDF5_OutputFile.Close();
}// end of Compute()
//------------------------------------------------------------------------------




/**
 * Print parameters of the simulation.
 * @param [in,out] file - where to print the parameters
 */
void TKSpaceFirstOrder3DSolver::PrintParametersOfSimulation(FILE * file)
{
  fprintf(file,"Domain dims:   [%4ld, %4ld,%4ld]\n",
                Parameters->GetFullDimensionSizes().X,
                Parameters->GetFullDimensionSizes().Y,
                Parameters->GetFullDimensionSizes().Z);

  fprintf(file,"Simulation time steps:  %8ld\n", Parameters->Get_Nt());
}// end of PrintParametersOfSimulation
//------------------------------------------------------------------------------



/**
 * Get peak memory usage.
 * @return Peak memory usage in MBs.
 *
 */
size_t TKSpaceFirstOrder3DSolver::ShowMemoryUsageInMB()
{
  // Linux build
  #ifdef __linux__
    struct rusage mem_usage;
    getrusage(RUSAGE_SELF, &mem_usage);

    return mem_usage.ru_maxrss >> 10;
  #endif

  // Windows build
  #ifdef _WIN64
    HANDLE hProcess;
    PROCESS_MEMORY_COUNTERS pmc;

    hProcess = OpenProcess(PROCESS_QUERY_INFORMATION | PROCESS_VM_READ,
                           FALSE,
                           GetCurrentProcessId());

    GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc));
    CloseHandle(hProcess);

    return pmc.PeakWorkingSetSize >> 20;
  #endif
}// end of ShowMemoryUsageInMB
//------------------------------------------------------------------------------


/**
 * Print Full code name and the license.
 * @param [in] file - file to print the data (stdout)
 */
void TKSpaceFirstOrder3DSolver::PrintFullNameCodeAndLicense(FILE * file)
{
  fprintf(file,"\n");
  fprintf(file,"+----------------------------------------------------+\n");
  fprintf(file,"| Build name:       kspaceFirstOrder3D v2.16         |\n");
  fprintf(file,"| Build date:       %*.*s                      |\n", 10,11,__DATE__);
  fprintf(file,"| Build time:       %*.*s                         |\n", 8,8,__TIME__);
  #if (defined (__KWAVE_GIT_HASH__))
    fprintf(file,"| Git hash: %s |\n",__KWAVE_GIT_HASH__);
  #endif
  fprintf(file,"|                                                    |\n");

  // OS detection
  #ifdef __linux__
    fprintf(file,"| Operating system: Linux x64                        |\n");
  #endif
  #ifdef _WIN64
    fprintf(file,"| Operating system: Windows x64                      |\n");
  #endif

  // Compiler detections
  #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    fprintf(file,"| Compiler name:    GNU C++ %.19s                    |\n", __VERSION__);
  #endif
  #ifdef __INTEL_COMPILER
    fprintf(file,"| Compiler name:    Intel C++ %d                   |\n", __INTEL_COMPILER);
  #endif

  // instruction set
  #if (defined (__AVX2__))
    fprintf(file,"| Instruction set:  Intel AVX 2                      |\n");
  #elif (defined (__AVX__))
    fprintf(file,"| Instruction set:  Intel AVX                        |\n");
  #elif (defined (__SSE4_2__))
    fprintf(file,"| Instruction set:  Intel SSE 4.2                    |\n");
  #elif (defined (__SSE4_1__))
    fprintf(file,"| Instruction set:  Intel SSE 4.1                    |\n");
  #elif (defined (__SSE3__))
    fprintf(file,"| Instruction set:  Intel SSE 3                      |\n");
  #elif (defined (__SSE2__))
    fprintf(file,"| Instruction set:  Intel SSE 2                      |\n");
  #endif

  fprintf(file,"|                                                    |\n");
  fprintf(file,"| Copyright (C) 2014 Jiri Jaros and Bradley Treeby   |\n");
  fprintf(file,"| http://www.k-wave.org                              |\n");
  fprintf(file,"+----------------------------------------------------+\n");
  fprintf(file,"\n");
}// end of GetFullCodeAndLincence
//------------------------------------------------------------------------------

/**
 * Set processor affinity.
 */
void TKSpaceFirstOrder3DSolver::SetProcessorAffinity()
{
  // Linux Build
  #ifdef __linux__
    //GNU compiler
    #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
      setenv("OMP_PROC_BIND","TRUE",1);
    #endif

    #ifdef __INTEL_COMPILER
      setenv("KMP_AFFINITY","none",1);
    #endif
  #endif

  // Windows build is always compiled by the Intel Compiler
  #ifdef _WIN64
    _putenv_s("KMP_AFFINITY","none");
  #endif
}//end of SetProcessorAffinity
//------------------------------------------------------------------------------


//----------------------------------------------------------------------------//
//                            Protected methods                               //
//----------------------------------------------------------------------------//

/**
 * Initialize FFTW plans.
 *
 */
void TKSpaceFirstOrder3DSolver::InitializeFFTWPlans()
{
  // initialization of FFTW library
  #ifdef _OPENMP
    fftwf_init_threads();
    fftwf_plan_with_nthreads(Parameters->GetNumberOfThreads());
  #endif

  // The simulation does not use checkpointing or this is the first turn
  bool RecoverFromPrevState = (Parameters->IsCheckpointEnabled() &&
                               THDF5_File::IsAccessible(Parameters->GetCheckpointFileName().c_str()));

  // import FFTW wisdom if it is here
  if (RecoverFromPrevState)
  {
     // try to find the wisdom in the file that has the same name as the checkpoint file (different extension)
    TFFTWComplexMatrix::ImportWisdom();
  }

  // create real to complex plans
  Get_FFT_X_temp().Create_FFT_Plan_3D_R2C(Get_p());
  Get_FFT_Y_temp().Create_FFT_Plan_3D_R2C(Get_p());
  Get_FFT_Z_temp().Create_FFT_Plan_3D_R2C(Get_p());

  // create real to complex plans
  Get_FFT_X_temp().Create_FFT_Plan_3D_C2R(Get_p());
  Get_FFT_Y_temp().Create_FFT_Plan_3D_C2R(Get_p());
  Get_FFT_Z_temp().Create_FFT_Plan_3D_C2R(Get_p());

  // if necessary, create 1D shift plans.
  // in this case, the matrix has a bit bigger dimensions to be able to store
  // shifted matrices.
  if (TParameters::GetInstance()->IsStore_u_non_staggered_raw())
  {
    // X shifts
    Get_FFT_shift_temp().Create_FFT_Plan_1DX_R2C(Get_p());
    Get_FFT_shift_temp().Create_FFT_Plan_1DX_C2R(Get_p());

    // Y shifts
    Get_FFT_shift_temp().Create_FFT_Plan_1DY_R2C(Get_p());
    Get_FFT_shift_temp().Create_FFT_Plan_1DY_C2R(Get_p());

    // Z shifts
    Get_FFT_shift_temp().Create_FFT_Plan_1DZ_R2C(Get_p());
    Get_FFT_shift_temp().Create_FFT_Plan_1DZ_C2R(Get_p());
  }// end u_non_staggered
}// end of InitializeFFTWPlans
//------------------------------------------------------------------------------



/**
 * Compute pre-processing phase \n
 * Initialize all indices, pre-compute constants such as c^2, rho0_sg* x dt
 * and create kappa, absorb_eta, absorb_tau, absorb_nabla1, absorb_nabla2 matrices.
 */
void TKSpaceFirstOrder3DSolver::PreProcessingPhase()
{
  // get the correct sensor mask and recompute indices
  if (Parameters->Get_sensor_mask_type() == TParameters::smt_index)
  {
    Get_sensor_mask_index().RecomputeIndicesToCPP();
  }

  if (Parameters->Get_sensor_mask_type() == TParameters::smt_corners)
  {
    Get_sensor_mask_corners().RecomputeIndicesToCPP();
  }


  if ((Parameters->Get_transducer_source_flag() != 0) ||
      (Parameters->Get_ux_source_flag() != 0)         ||
      (Parameters->Get_uy_source_flag() != 0)         ||
      (Parameters->Get_uz_source_flag() != 0)
     )
  {
    Get_u_source_index().RecomputeIndicesToCPP();
  }

  if (Parameters->Get_transducer_source_flag() != 0)
  {
    Get_delay_mask().RecomputeIndicesToCPP();
  }

  if (Parameters->Get_p_source_flag() != 0)
  {
    Get_p_source_index().RecomputeIndicesToCPP();
  }


  // compute dt / rho0_sg...
  if (Parameters->Get_rho0_scalar_flag())
  { // rho is scalar
    Parameters->Get_rho0_sgx_scalar() = Parameters->Get_dt() / Parameters->Get_rho0_sgx_scalar();
    Parameters->Get_rho0_sgy_scalar() = Parameters->Get_dt() / Parameters->Get_rho0_sgy_scalar();
    Parameters->Get_rho0_sgz_scalar() = Parameters->Get_dt() / Parameters->Get_rho0_sgz_scalar();
  }
  else
  { // non-uniform grid cannot be pre-calculated :-(
    // rho is matrix
    if (Parameters->Get_nonuniform_grid_flag())
    {
      Calculate_dt_rho0_non_uniform();
    }
    else
    {
      Get_dt_rho0_sgx().ScalarDividedBy(Parameters->Get_dt());
      Get_dt_rho0_sgy().ScalarDividedBy(Parameters->Get_dt());
      Get_dt_rho0_sgz().ScalarDividedBy(Parameters->Get_dt());
    }
  }

  // generate different matrices
  if (Parameters->Get_absorbing_flag() != 0)
  {
    Generate_kappa_absorb_nabla1_absorb_nabla2();
    Generate_absorb_tau_absorb_eta_matrix();
  }
  else
  {
    Generate_kappa();
  }

  // calculate c^2. It has to be after kappa gen... because of c modification
  Compute_c2();
}// end of PreProcessingPhase
//------------------------------------------------------------------------------


/**
 * Generate kappa matrix for lossless case.
 *
 */
void TKSpaceFirstOrder3DSolver::Generate_kappa()
{
  #pragma omp parallel
  {
    const float dx_sq_rec = 1.0f / (Parameters->Get_dx()*Parameters->Get_dx());
    const float dy_sq_rec = 1.0f / (Parameters->Get_dy()*Parameters->Get_dy());
    const float dz_sq_rec = 1.0f / (Parameters->Get_dz()*Parameters->Get_dz());

    const float c_ref_dt_pi = Parameters->Get_c_ref() * Parameters->Get_dt() * float(M_PI);

    const float Nx_rec   = 1.0f / (float) Parameters->GetFullDimensionSizes().X;
    const float Ny_rec   = 1.0f / (float) Parameters->GetFullDimensionSizes().Y;
    const float Nz_rec   = 1.0f / (float) Parameters->GetFullDimensionSizes().Z;


    const size_t X_Size  = Parameters->GetReducedDimensionSizes().X;
    const size_t Y_Size  = Parameters->GetReducedDimensionSizes().Y;
    const size_t Z_Size  = Parameters->GetReducedDimensionSizes().Z;

    float * kappa = Get_kappa().GetRawData();

    #pragma omp for schedule (static)
    for (size_t z = 0; z < Z_Size; z++)
    {
      const float z_f = (float) z;
      float z_part = 0.5f - fabs(0.5f - z_f * Nz_rec );
      z_part = (z_part * z_part) * dz_sq_rec;

      for (size_t y = 0; y < Y_Size; y++)
      {
        const float y_f = (float) y;
        float y_part = 0.5f - fabs(0.5f - y_f * Ny_rec);
        y_part = (y_part * y_part) * dy_sq_rec;

        const float yz_part = z_part + y_part;
        for (size_t x = 0; x < X_Size; x++)
        {
          const float x_f = (float) x;
          float x_part = 0.5f - fabs(0.5f - x_f * Nx_rec);
          x_part = (x_part * x_part) * dx_sq_rec;

          float  k = c_ref_dt_pi * sqrtf(x_part + yz_part);

          // kappa element
          kappa[(z*Y_Size + y) * X_Size + x ] = (k == 0.0f) ? 1.0f : sin(k)/k;
        }//x
      }//y
    }// z
  }// parallel
}// end of Generate_kappa
//------------------------------------------------------------------------------

/**
 * Generate kappa, absorb_nabla1, absorb_nabla2 for absorbing media.
 *
 */
void TKSpaceFirstOrder3DSolver::Generate_kappa_absorb_nabla1_absorb_nabla2()
{
  #pragma omp parallel
  {
    const float dx_sq_rec = 1.0f / (Parameters->Get_dx()*Parameters->Get_dx());
    const float dy_sq_rec = 1.0f / (Parameters->Get_dy()*Parameters->Get_dy());
    const float dz_sq_rec = 1.0f / (Parameters->Get_dz()*Parameters->Get_dz());

    const float c_ref_dt_2 = Parameters->Get_c_ref() * Parameters->Get_dt() * 0.5f;
    const float pi_2       = float(M_PI) * 2.0f;

    const size_t Nx = Parameters->GetFullDimensionSizes().X;
    const size_t Ny = Parameters->GetFullDimensionSizes().Y;
    const size_t Nz = Parameters->GetFullDimensionSizes().Z;

    const float Nx_rec   = 1.0f / (float) Nx;
    const float Ny_rec   = 1.0f / (float) Ny;
    const float Nz_rec   = 1.0f / (float) Nz;

    const size_t X_Size  = Parameters->GetReducedDimensionSizes().X;
    const size_t Y_Size  = Parameters->GetReducedDimensionSizes().Y;
    const size_t Z_Size  = Parameters->GetReducedDimensionSizes().Z;

    float * kappa           = Get_kappa().GetRawData();
    float * absorb_nabla1   = Get_absorb_nabla1().GetRawData();
    float * absorb_nabla2   = Get_absorb_nabla2().GetRawData();
    const float alpha_power = Parameters->Get_alpha_power();

    #pragma omp for schedule (static)
    for (size_t z = 0; z < Z_Size; z++)
    {
      const float z_f = (float) z;
      float z_part = 0.5f - fabs(0.5f - z_f * Nz_rec );
      z_part = (z_part * z_part) * dz_sq_rec;

      for (size_t y = 0; y < Y_Size; y++)
      {
        const float y_f = (float) y;
        float y_part = 0.5f - fabs(0.5f - y_f * Ny_rec);
        y_part = (y_part * y_part) * dy_sq_rec;

        const float yz_part = z_part + y_part;

        size_t i = (z*Y_Size + y) * X_Size;

        for (size_t x = 0; x < X_Size; x++)
        {
          const float x_f = (float) x;

          float x_part = 0.5f - fabs(0.5f - x_f * Nx_rec);
          x_part = (x_part * x_part) * dx_sq_rec;


          float  k         = pi_2 * sqrt(x_part + yz_part);
          float  c_ref_k   = c_ref_dt_2 * k;

          absorb_nabla1[i] = pow(k, alpha_power - 2);
          absorb_nabla2[i] = pow(k, alpha_power - 1);

          kappa[i]         =  (c_ref_k == 0.0f) ? 1.0f : sin(c_ref_k)/c_ref_k;

          if (absorb_nabla1[i] ==  std::numeric_limits<float>::infinity()) absorb_nabla1[i] = 0.0f;
          if (absorb_nabla2[i] ==  std::numeric_limits<float>::infinity()) absorb_nabla2[i] = 0.0f;

          i++;
        }//x
      }//y
    }// z
  }// parallel
}// end of Generate_kappa_absorb_nabla1_absorb_nabla2
//------------------------------------------------------------------------------

/**
 * Generate absorb_tau and absorb_eta in for heterogenous media.
 */
void TKSpaceFirstOrder3DSolver::Generate_absorb_tau_absorb_eta_matrix()
{
  // test for scalars
  if ((Parameters->Get_alpha_coeff_scallar_flag()) && (Parameters->Get_c0_scalar_flag()))
  {
    const float alpha_power = Parameters->Get_alpha_power();
    const float tan_pi_y_2  = tan(float(M_PI_2)* alpha_power);
    const float alpha_db_neper_coeff = (100.0f * pow(1.0e-6f / (2.0f * (float) M_PI), alpha_power)) /
                                       (20.0f * (float) M_LOG10E);

    const float alpha_coeff_2 = 2.0f * Parameters->Get_alpha_coeff_scallar() * alpha_db_neper_coeff;

    Parameters->Get_absorb_tau_scalar() =  (-alpha_coeff_2) * pow(Parameters->Get_c0_scalar(),alpha_power-1);
    Parameters->Get_absorb_eta_scalar() =    alpha_coeff_2  * pow(Parameters->Get_c0_scalar(),alpha_power) * tan_pi_y_2;
  }
  else
  {
    #pragma omp parallel
    {
      const size_t Z_Size  = Parameters->GetFullDimensionSizes().Z;
      const size_t Y_Size  = Parameters->GetFullDimensionSizes().Y;
      const size_t X_Size  = Parameters->GetFullDimensionSizes().X;

      float * absorb_tau = Get_absorb_tau().GetRawData();
      float * absorb_eta = Get_absorb_eta().GetRawData();

      float * alpha_coeff;
      size_t  alpha_shift;

      if (Parameters->Get_alpha_coeff_scallar_flag())
      {
        alpha_coeff = &(Parameters->Get_alpha_coeff_scallar());
        alpha_shift = 0;
      }
      else
      {
        alpha_coeff = Get_Temp_1_RS3D().GetRawData();
        alpha_shift = 1;
      }

      float * c0;
      size_t  c0_shift;
      if (Parameters->Get_c0_scalar_flag())
      {
        c0 = &(Parameters->Get_c0_scalar());
        c0_shift = 0;
      }
      else
      {
        c0 = Get_c2().GetRawData();
        c0_shift = 1;
      }

      const float alpha_power = Parameters->Get_alpha_power();
      const float tan_pi_y_2  = tan(float(M_PI_2)* alpha_power);

      //alpha = 100*alpha.*(1e-6/(2*pi)).^y./
      //                  (20*log10(exp(1)));
      const float alpha_db_neper_coeff = (100.0f * pow(1.0e-6f / (2.0f * (float) M_PI), alpha_power)) /
                                         (20.0f * (float) M_LOG10E);

      #pragma omp for schedule (static)
      for (size_t z = 0; z < Z_Size; z++)
      {
        for (size_t y = 0; y < Y_Size; y++)
        {
          size_t i = (z*Y_Size + y) * X_Size;
          for (size_t x = 0; x < X_Size; x++)
          {
            const float alpha_coeff_2 = 2.0f * alpha_coeff[i * alpha_shift] * alpha_db_neper_coeff;
            absorb_tau[i] = (-alpha_coeff_2) * pow(c0[i * c0_shift],alpha_power-1);
            absorb_eta[i] =   alpha_coeff_2  * pow(c0[i * c0_shift],alpha_power) * tan_pi_y_2;
            i++;
          }//x
        }//y
      }// z
    }// parallel
  } // absorb_tau and aborb_eta = matrics
}// end of Generate_absorb_tau_absorb_eta_matrix
//------------------------------------------------------------------------------

/**
 * Prepare dt./ rho0  for non-uniform grid.
 *
 */
void TKSpaceFirstOrder3DSolver::Calculate_dt_rho0_non_uniform()
{
  #pragma omp parallel
  {
    float * dt_rho0_sgx   = Get_dt_rho0_sgx().GetRawData();
    float * dt_rho0_sgy   = Get_dt_rho0_sgy().GetRawData();
    float * dt_rho0_sgz   = Get_dt_rho0_sgz().GetRawData();

    const float dt = Parameters->Get_dt();

    const float * duxdxn_sgx = Get_dxudxn_sgx().GetRawData();
    const float * duydyn_sgy = Get_dyudyn_sgy().GetRawData();
    const float * duzdzn_sgz = Get_dzudzn_sgz().GetRawData();

    const size_t Z_Size = Get_dt_rho0_sgx().GetDimensionSizes().Z;
    const size_t Y_Size = Get_dt_rho0_sgx().GetDimensionSizes().Y;
    const size_t X_Size = Get_dt_rho0_sgx().GetDimensionSizes().X;

    const size_t SliceSize = (X_Size * Y_Size );

    #pragma omp for schedule (static)
    for (size_t z = 0; z < Z_Size; z++)
    {
      register size_t i = z * SliceSize;
      for (size_t y = 0; y < Y_Size; y++)
      {
        for (size_t x = 0; x < X_Size; x++)
        {
          dt_rho0_sgx[i] = (dt * duxdxn_sgx[x]) / dt_rho0_sgx[i];
          i++;
        } // x
      } // y
    } // z

    #pragma omp for schedule (static)
    for (size_t z = 0; z < Z_Size; z++)
    {
      register size_t i = z * SliceSize;
      for (size_t y = 0; y < Y_Size; y++)
      {
        const float duydyn_el = duydyn_sgy[y];
        for (size_t x = 0; x < X_Size; x++)
        {
          dt_rho0_sgy[i] = (dt * duydyn_el) / dt_rho0_sgy[i];
          i++;
        } // x
      } // y
    } // z

    #pragma omp for schedule (static)
    for (size_t z = 0; z < Z_Size; z++)
    {
      register size_t i = z* SliceSize;
      const float duzdzn_el = duzdzn_sgz[z];
      for (size_t y = 0; y < Y_Size; y++)
      {
        for (size_t x = 0; x < X_Size; x++)
        {
          dt_rho0_sgz[i] = (dt * duzdzn_el) / dt_rho0_sgz[i];
          i++;
        } // x
      } // y
    } // z
  } // parallel
}// end of Calculate_dt_rho0_non_uniform
//------------------------------------------------------------------------------

/**
 * Calculate p0 source when necessary.
 *
 */
void TKSpaceFirstOrder3DSolver::Calculate_p0_source()
{
  Get_p().CopyData(Get_p0_source_input());

  const float * p0 = Get_p0_source_input().GetRawData();

  float * c2;
  size_t c2_shift;

  if (Parameters->Get_c0_scalar_flag())
  {
    c2 = &Parameters->Get_c0_scalar();
    c2_shift = 0;
  }
  else
  {
    c2 = Get_c2().GetRawData();
    c2_shift = 1;
  }

  float * rhox = Get_rhox().GetRawData();
  float * rhoy = Get_rhoy().GetRawData();
  float * rhoz = Get_rhoz().GetRawData();

  //  add the initial pressure to rho as a mass source
  float tmp;

  #pragma omp parallel for schedule (static) private(tmp)
  for (size_t i = 0; i < Get_rhox().GetTotalElementCount(); i++)
  {
    tmp = p0[i] / (3.0f* c2[i * c2_shift]);
    rhox[i] = tmp;
    rhoy[i] = tmp;
    rhoz[i] = tmp;
  }

  //------------------------------------------------------------------------//
  //--  compute u(t = t1 + dt/2) based on the assumption u(dt/2) = -u(-dt/2) --//
  //--    which forces u(t = t1) = 0 --//
  //------------------------------------------------------------------------//
  Compute_ddx_kappa_fft_p(Get_p(),
                          Get_FFT_X_temp(),Get_FFT_Y_temp(),Get_FFT_Z_temp(),
                          Get_kappa(),
                          Get_ddx_k_shift_pos(),Get_ddy_k_shift_pos(),Get_ddz_k_shift_pos()
                          );

  if (Parameters->Get_rho0_scalar_flag())
  {
    if (Parameters->Get_nonuniform_grid_flag())
    { // non uniform grid
      Get_ux_sgx().Compute_dt_rho_sg_mul_ifft_div_2_scalar_nonuniform_x(Parameters->Get_rho0_sgx_scalar(),
                                                                        Get_dxudxn_sgx(),
                                                                        Get_FFT_X_temp());
      Get_uy_sgy().Compute_dt_rho_sg_mul_ifft_div_2_scalar_nonuniform_y(Parameters->Get_rho0_sgy_scalar(),
                                                                        Get_dyudyn_sgy(),
                                                                        Get_FFT_Y_temp());
      Get_uz_sgz().Compute_dt_rho_sg_mul_ifft_div_2_scalar_nonuniform_z(Parameters->Get_rho0_sgz_scalar(),
                                                                        Get_dzudzn_sgz(),
                                                                        Get_FFT_Z_temp());
    }
    else
    { //uniform grid, heterogeneous
      Get_ux_sgx().Compute_dt_rho_sg_mul_ifft_div_2(Parameters->Get_rho0_sgx_scalar(), Get_FFT_X_temp());
      Get_uy_sgy().Compute_dt_rho_sg_mul_ifft_div_2(Parameters->Get_rho0_sgy_scalar(), Get_FFT_Y_temp());
      Get_uz_sgz().Compute_dt_rho_sg_mul_ifft_div_2(Parameters->Get_rho0_sgz_scalar(), Get_FFT_Z_temp());
    }
  }
  else
  { // homogeneous, non-unifrom grid
    // divide the matrix by 2 and multiply with st./rho0_sg
    Get_ux_sgx().Compute_dt_rho_sg_mul_ifft_div_2(Get_dt_rho0_sgx(), Get_FFT_X_temp());
    Get_uy_sgy().Compute_dt_rho_sg_mul_ifft_div_2(Get_dt_rho0_sgy(), Get_FFT_Y_temp());
    Get_uz_sgz().Compute_dt_rho_sg_mul_ifft_div_2(Get_dt_rho0_sgz(), Get_FFT_Z_temp());
  }
}// end of Calculate_p0_source
//------------------------------------------------------------------------------

/**
 * Compute c^2.
 *
 */
void TKSpaceFirstOrder3DSolver::Compute_c2()
{
  if (Parameters->Get_c0_scalar_flag())
  { //scalar
    float c = Parameters->Get_c0_scalar();
    Parameters->Get_c0_scalar() = c * c;
  }
  else
  {
    float * c2 =  Get_c2().GetRawData();

    #pragma omp parallel for schedule (static)
    for (size_t i=0; i< Get_c2().GetTotalElementCount(); i++)
    {
      c2[i] = c2[i] * c2[i];
    }
  }// matrix
}// ComputeC2
//------------------------------------------------------------------------------

/**
 *  Compute part of the new velocity term - gradient of p
 *  represented by:
 *  bsxfun(\@times, ddx_k_shift_pos, kappa .* p_k)
 *
 * @param [in]  X_Matrix - 3D pressure matrix
 * @param [out] FFT_X - matrix to store input for iFFT (p) /dx
 * @param [out] FFT_Y - matrix to store input for iFFT (p) /dy
 * @param [out] FFT_Z - matrix to store input for iFFT (p) /dz
 *
 * @param [in]  kappa - Real matrix of kappa
 *
 * @param [in]  ddx - precomputed value of ddx_k_shift_pos
 * @param [in]  ddy - precomputed value of ddy_k_shift_pos
 * @param [in]  ddz - precomputed value of ddz_k_shift_pos
 */
void TKSpaceFirstOrder3DSolver::Compute_ddx_kappa_fft_p(TRealMatrix&        X_Matrix,
                                                        TFFTWComplexMatrix& FFT_X,
                                                        TFFTWComplexMatrix& FFT_Y,
                                                        TFFTWComplexMatrix& FFT_Z,
                                                        const TRealMatrix&  kappa,
                                                        const TComplexMatrix& ddx,
                                                        const TComplexMatrix& ddy,
                                                        const TComplexMatrix& ddz)
{
  // Compute FFT of X
  FFT_X.Compute_FFT_3D_R2C(X_Matrix);

  #pragma omp parallel
  {
    float * p_k_x_data = FFT_X.GetRawData();
    float * p_k_y_data = FFT_Y.GetRawData();
    float * p_k_z_data = FFT_Z.GetRawData();

    const float * kappa_data = kappa.GetRawData();
    const float * ddx_data   = ddx.GetRawData();
    const float * ddy_data   = ddy.GetRawData();
    const float * ddz_data   = ddz.GetRawData();

    const size_t Z_Size = FFT_X.GetDimensionSizes().Z;
    const size_t Y_Size = FFT_X.GetDimensionSizes().Y;
    const size_t X_Size = FFT_X.GetDimensionSizes().X;

    const size_t SliceSize = (X_Size * Y_Size ) << 1;

    #pragma omp for schedule (static)
    for (size_t z = 0; z < Z_Size; z++)
    {
      register size_t i = z * SliceSize;
      const size_t z2 = z<<1;
      for (size_t y = 0; y < Y_Size; y++)
      {
        const size_t y2 = y<<1;
        for (size_t x = 0; x < X_Size;  x++)
        {
          // kappa ./ p_k
          const float kappa_el  = kappa_data[i>>1];
          const float p_k_el_re = p_k_x_data[i]   * kappa_el;
          const float p_k_el_im = p_k_x_data[i+1] * kappa_el;
          const size_t x2 = x<<1;

          //bsxfun(ddx...)
          p_k_x_data[i]   = p_k_el_re * ddx_data[x2]   - p_k_el_im * ddx_data[x2+1];
          p_k_x_data[i+1] = p_k_el_re * ddx_data[x2+1] + p_k_el_im * ddx_data[x2];

          //bsxfun(ddy...)
          p_k_y_data[i]   = p_k_el_re * ddy_data[y2]   - p_k_el_im * ddy_data[y2+1];
          p_k_y_data[i+1] = p_k_el_re * ddy_data[y2+1] + p_k_el_im * ddy_data[y2];

          //bsxfun(ddz...)
          p_k_z_data[i]   = p_k_el_re * ddz_data[z2]   - p_k_el_im * ddz_data[z2+1];
          p_k_z_data[i+1] = p_k_el_re * ddz_data[z2+1] + p_k_el_im * ddz_data[z2];

          i +=2;
        } // x
      } // y
    } // z
  }// parallel
}// end of KSpaceFirstOrder3DSolver
//------------------------------------------------------------------------------


/**
 * Compute new values for duxdx, duydy, duzdz.
 *
 */
void  TKSpaceFirstOrder3DSolver::Compute_duxyz()
{
  Get_FFT_X_temp().Compute_FFT_3D_R2C(Get_ux_sgx());
  Get_FFT_Y_temp().Compute_FFT_3D_R2C(Get_uy_sgy());
  Get_FFT_Z_temp().Compute_FFT_3D_R2C(Get_uz_sgz());

  #pragma omp parallel
  {
    float * Temp_FFT_X_Data  = Get_FFT_X_temp().GetRawData();
    float * Temp_FFT_Y_Data  = Get_FFT_Y_temp().GetRawData();
    float * Temp_FFT_Z_Data  = Get_FFT_Z_temp().GetRawData();

    const float * kappa   = Get_kappa().GetRawData();

    const size_t FFT_Z_dim = Get_FFT_X_temp().GetDimensionSizes().Z;
    const size_t FFT_Y_dim = Get_FFT_X_temp().GetDimensionSizes().Y;
    const size_t FFT_X_dim = Get_FFT_X_temp().GetDimensionSizes().X;

    const size_t SliceSize = (FFT_X_dim * FFT_Y_dim) << 1;
    const float  Divider = 1.0f / Get_ux_sgx().GetTotalElementCount();

    const TFloatComplex * ddx = (TFloatComplex *) Get_ddx_k_shift_neg().GetRawData();
    const TFloatComplex * ddy = (TFloatComplex *) Get_ddy_k_shift_neg().GetRawData();
    const TFloatComplex * ddz = (TFloatComplex *) Get_ddz_k_shift_neg().GetRawData();


    #pragma omp for schedule (static)
    for (size_t z = 0; z < FFT_Z_dim; z++)
    {
      register size_t i = z * SliceSize;

      const float ddz_neg_re = ddz[z].real;
      const float ddz_neg_im = ddz[z].imag;
      for (size_t y = 0; y < FFT_Y_dim; y++)
      {
        const float ddy_neg_re = ddy[y].real;
        const float ddy_neg_im = ddy[y].imag;
        for (size_t x = 0; x < FFT_X_dim; x++)
        {
          float FFT_el_x_re = Temp_FFT_X_Data[i];
          float FFT_el_x_im = Temp_FFT_X_Data[i+1];

          float FFT_el_y_re = Temp_FFT_Y_Data[i];
          float FFT_el_y_im = Temp_FFT_Y_Data[i+1];

          float FFT_el_z_re = Temp_FFT_Z_Data[i];
          float FFT_el_z_im = Temp_FFT_Z_Data[i+1];

          const float kappa_el = kappa[i >> 1];

          FFT_el_x_re   *= kappa_el;
          FFT_el_x_im   *= kappa_el;

          FFT_el_y_re   *= kappa_el;
          FFT_el_y_im   *= kappa_el;

          FFT_el_z_re   *= kappa_el;
          FFT_el_z_im   *= kappa_el;


          Temp_FFT_X_Data[i]     = ((FFT_el_x_re * ddx[x].real) -
                                    (FFT_el_x_im * ddx[x].imag)
                                   ) * Divider;
          Temp_FFT_X_Data[i + 1] = ((FFT_el_x_im * ddx[x].real) +
                                    (FFT_el_x_re * ddx[x].imag)
                                   )* Divider;

          Temp_FFT_Y_Data[i]     = ((FFT_el_y_re * ddy_neg_re) -
                                    (FFT_el_y_im * ddy_neg_im)
                                   ) * Divider;
          Temp_FFT_Y_Data[i + 1] = ((FFT_el_y_im * ddy_neg_re) +
                                    (FFT_el_y_re * ddy_neg_im)
                                   )* Divider;

          Temp_FFT_Z_Data[i]     = ((FFT_el_z_re * ddz_neg_re) -
                                    (FFT_el_z_im * ddz_neg_im)
                                   ) * Divider;
          Temp_FFT_Z_Data[i + 1] = ((FFT_el_z_im * ddz_neg_re) +
                                    (FFT_el_z_re * ddz_neg_im)
                                   )* Divider;

          i+=2;
        } // x
      } // y
    } // z
  } // parallel;

  Get_FFT_X_temp().Compute_FFT_3D_C2R(Get_duxdx());
  Get_FFT_Y_temp().Compute_FFT_3D_C2R(Get_duydy());
  Get_FFT_Z_temp().Compute_FFT_3D_C2R(Get_duzdz());

 //-------------------------------------------------------------------------//
 //--------------------- Non linear grid -----------------------------------//
 //-------------------------------------------------------------------------//
  if (Parameters->Get_nonuniform_grid_flag() != 0)
  {
    #pragma omp parallel
    {
      float * duxdx = Get_duxdx().GetRawData();
      float * duydy = Get_duydy().GetRawData();
      float * duzdz = Get_duzdz().GetRawData();

      const float * duxdxn = Get_dxudxn().GetRawData();
      const float * duydyn = Get_dyudyn().GetRawData();
      const float * duzdzn = Get_dzudzn().GetRawData();

      const size_t Z_Size = Get_duxdx().GetDimensionSizes().Z;
      const size_t Y_Size = Get_duxdx().GetDimensionSizes().Y;
      const size_t X_Size = Get_duxdx().GetDimensionSizes().X;

      const size_t SliceSize = (X_Size * Y_Size );

      #pragma omp for schedule (static)
      for (size_t z = 0; z < Z_Size; z++)
      {
        register size_t i = z* SliceSize;
        for (size_t y = 0; y < Y_Size; y++)
        {
          for (size_t x = 0; x < X_Size; x++)
          {
            duxdx[i] *= duxdxn[x];
            i++;
          } // x
        } // y
      } // z

      #pragma omp for schedule (static)
      for (size_t z = 0; z < Z_Size; z++)
      {
        register size_t i = z * SliceSize;
        for (size_t y = 0; y < Y_Size; y++)
        {
          const float dyudyn_el = duydyn[y];
          for (size_t x = 0; x < X_Size; x++)
          {
            duydy[i] *=  dyudyn_el;
            i++;
          } // x
        } // y
      } // z

      #pragma omp for schedule (static)
      for (size_t z = 0; z < Z_Size; z++)
      {
        const float duzdzn_el = duzdzn[z];
        register size_t i = z * SliceSize;
        for (size_t y = 0; y < Y_Size; y++)
        {
          for (size_t x = 0; x < X_Size; x++)
          {
            duzdz[i] *=  duzdzn_el;
            i++;
          } // x
        } // y
      } // z
    } // parallel
 }// nonlinear
}// end of Compute_duxyz
//------------------------------------------------------------------------------


/**
 * Calculate new values of rhox, rhoy and rhoz for non-linear case.
 *
 */
void TKSpaceFirstOrder3DSolver::Compute_rhoxyz_nonlinear()
{
  const size_t Z_Size = Get_rhox().GetDimensionSizes().Z;
  const size_t Y_Size = Get_rhox().GetDimensionSizes().Y;
  const size_t X_Size = Get_rhox().GetDimensionSizes().X;

  const float dt2   = 2.0f * Parameters->Get_dt();
  const float dt_el = Parameters->Get_dt();
  const size_t SliceSize = Y_Size * X_Size;

  #pragma omp parallel
  {
    float * rhox_data  = Get_rhox().GetRawData();
    float * rhoy_data  = Get_rhoy().GetRawData();
    float * rhoz_data  = Get_rhoz().GetRawData();

    const float * pml_x_data = Get_pml_x().GetRawData();
    const float * duxdx_data = Get_duxdx().GetRawData();
    const float * duydy_data = Get_duydy().GetRawData();
    const float * duzdz_data = Get_duzdz().GetRawData();

    // rho0 is a scalar
    if (Parameters->Get_rho0_scalar())
    {
      const float dt_rho0 = Parameters->Get_rho0_scalar() * dt_el;

      #pragma omp for schedule (static)
      for (size_t z = 0; z < Z_Size; z++)
      {
        register size_t i = z * SliceSize;
        for (size_t y = 0; y < Y_Size; y++)
        {
          for (size_t x = 0; x < X_Size; x++)
          {
            const float pml_x   = pml_x_data[x];
                  float du      = duxdx_data[i];

            rhox_data[i] = pml_x * (
                                ((pml_x * rhox_data[i]) - (dt_rho0 * du))/
                                (1.0f + (dt2 * du))
                                );
            i++;
          } // x
        }// y
      }// z

      #pragma omp for schedule (static)
      for (size_t z = 0; z < Z_Size; z++)
      {
        register size_t i = z * SliceSize;
        for (size_t y = 0; y < Y_Size; y++)
        {
          const float pml_y = Get_pml_y()[y];
          for (size_t x = 0; x < X_Size; x++)
          {
            float du = duydy_data[i];

            rhoy_data[i] = pml_y * (
                                ((pml_y * rhoy_data[i]) - (dt_rho0 * du))/
                                (1.0f + (dt2 * du))
                                );
            i++;
          } // x
        }// y
      }// z


      #pragma omp for schedule (static)
      for (size_t z = 0; z < Z_Size; z++)
      {
        register size_t i = z * SliceSize;
        const float pml_z = Get_pml_z()[z];
        for (size_t y = 0; y < Y_Size; y++)
        {
          for (size_t x = 0; x < X_Size; x++)
          {
            float du      = duzdz_data[i];

            rhoz_data[i] = pml_z * (
                                ((pml_z * rhoz_data[i]) - (dt_rho0 * du))/
                                (1.0f + (dt2 * du))
                                );
            i++;
          } // x
        }// y
      }// z
    }
    else
    { //----------------------------------------------------------------//
      // rho0 is a matrix
      const float * rho0_data  = Get_rho0().GetRawData();

      #pragma omp for schedule (static)
      for (size_t z = 0; z < Z_Size; z++)
      {
        register size_t i = z * SliceSize;
        for (size_t y = 0; y < Y_Size; y++)
        {
          for (size_t x = 0; x < X_Size; x++)
          {
            const float pml_x   = pml_x_data[x];
            const float dt_rho0 = dt_el * rho0_data[i];
                  float du      = duxdx_data[i];

            rhox_data[i] = pml_x * (
                                ((pml_x * rhox_data[i]) - (dt_rho0 * du))/
                                (1.0f + (dt2 * du))
                                );
            i++;
          } // x
        }// y
      }// z

      #pragma omp for schedule (static)
      for (size_t z = 0; z < Z_Size; z++)
      {
        register size_t i = z * SliceSize;
        for (size_t y = 0; y < Y_Size; y++)
        {
          const float pml_y = Get_pml_y()[y];
          for (size_t x = 0; x < X_Size; x++)
          {
            const float dt_rho0 = dt_el * rho0_data[i];
                  float du = duydy_data[i];

            rhoy_data[i] = pml_y * (
                                ((pml_y * rhoy_data[i]) - (dt_rho0 * du))/
                                (1.0f + (dt2 * du))
                                );
            i++;
          } // x
        }// y
      }// z

      #pragma omp for schedule (static)
      for (size_t z = 0; z < Z_Size; z++)
      {
        register size_t i = z * SliceSize;
        const float pml_z = Get_pml_z()[z];
        for (size_t y = 0; y < Y_Size; y++)
        {
          for (size_t x = 0; x < X_Size; x++)
          {
            const float dt_rho0 = dt_el * rho0_data[i];
                  float du      = duzdz_data[i];

            rhoz_data[i] = pml_z * (
                                ((pml_z * rhoz_data[i]) - (dt_rho0 * du))/
                                (1.0f + (dt2 * du))
                                );
            i++;
          } // x
        }// y
      }// z
    } // end rho is matrix
  }// parallel
}// end of Compute_rhoxyz_nonlinear
//------------------------------------------------------------------------------


/**
 * Calculate new values of rhox, rhoy and rhoz for linear case.
 *
 */
void TKSpaceFirstOrder3DSolver::Compute_rhoxyz_linear()
{
  const size_t Z_Size = Get_rhox().GetDimensionSizes().Z;
  const size_t Y_Size = Get_rhox().GetDimensionSizes().Y;
  const size_t X_Size = Get_rhox().GetDimensionSizes().X;

  const float dt_el = Parameters->Get_dt();
  const size_t SliceSize =  Y_Size * X_Size;

  #pragma omp parallel
  {
    float * rhox_data  = Get_rhox().GetRawData();
    float * rhoy_data  = Get_rhoy().GetRawData();
    float * rhoz_data  = Get_rhoz().GetRawData();

    const float * pml_x_data = Get_pml_x().GetRawData();
    const float * duxdx_data = Get_duxdx().GetRawData();
    const float * duydy_data = Get_duydy().GetRawData();
    const float * duzdz_data = Get_duzdz().GetRawData();


    if (Parameters->Get_rho0_scalar())
    { // rho0 is a scalar
      const float dt_rho0 = Parameters->Get_rho0_scalar() * dt_el;

      #pragma omp for schedule (static)
      for (size_t z = 0; z < Z_Size; z++)
      {
        register size_t i = z * SliceSize;
        for (size_t y = 0; y < Y_Size; y++)
        {
          for (size_t x = 0; x < X_Size; x++)
          {
            const float pml_x   = pml_x_data[x];

            rhox_data[i] = pml_x * (
                                   ((pml_x * rhox_data[i]) - (dt_rho0 * duxdx_data[i]))
                                   );
            i++;
          } // x
        }// y
      }// z

      #pragma omp for schedule (static)
      for (size_t z = 0; z < Z_Size; z++)
      {
        register size_t i = z * SliceSize;
        for (size_t y = 0; y < Y_Size; y++)
        {
          const float pml_y = Get_pml_y()[y];
          for (size_t x = 0; x < X_Size; x++)
          {
            rhoy_data[i] = pml_y * (
                                   ((pml_y * rhoy_data[i]) - (dt_rho0 * duydy_data[i]))
                                   );
            i++;
          } // x
        }// y
      }// z

      #pragma omp for schedule (static)
      for (size_t z = 0; z < Z_Size; z++)
      {
        register size_t i = z * SliceSize;
        const float pml_z = Get_pml_z()[z];

        for (size_t y = 0; y < Y_Size; y++)
        {
          for (size_t x = 0; x < X_Size; x++)
          {
            rhoz_data[i] = pml_z * (
                                   ((pml_z * rhoz_data[i]) - (dt_rho0 * duzdz_data[i]))
                                   );
            i++;
          } // x
        }// y
      }// z

    }
    else
    { //-----------------------------------------------------//
      // rho0 is a matrix
      const float * rho0_data  = Get_rho0().GetRawData();

      #pragma omp for schedule (static)
      for (size_t z = 0; z < Z_Size; z++)
      {
        register size_t i = z * SliceSize;
        for (size_t y = 0; y < Y_Size; y++)
        {
          for (size_t x = 0; x < X_Size; x++)
          {
            const float pml_x   = pml_x_data[x];
            const float dt_rho0 = dt_el * rho0_data[i];

            rhox_data[i] = pml_x * (
                                   ((pml_x * rhox_data[i]) - (dt_rho0 * duxdx_data[i]))
                                   );

            i++;
          } // x
        }// y
      }// z

      #pragma omp for schedule (static)
      for (size_t z = 0; z < Z_Size; z++)
      {
        register size_t i = z * SliceSize;
        for (size_t y = 0; y < Y_Size; y++)
        {
          const float pml_y = Get_pml_y()[y];
          for (size_t x = 0; x < X_Size; x++)
          {
            const float dt_rho0 = dt_el * rho0_data[i];

            rhoy_data[i] = pml_y * (
                                   ((pml_y * rhoy_data[i]) - (dt_rho0 * duydy_data[i]))
                                   );
            i++;

          } // x
        }// y
      }// z


      #pragma omp for schedule (static)
      for (size_t z = 0; z < Z_Size; z++)
      {
        register size_t i = z * SliceSize;
        const float pml_z = Get_pml_z()[z];

        for (size_t y = 0; y < Y_Size; y++)
        {
          for (size_t x = 0; x < X_Size; x++)
          {
            const float dt_rho0 = dt_el * rho0_data[i];

            rhoz_data[i] = pml_z * (
                                   ((pml_z * rhoz_data[i]) - (dt_rho0 * duzdz_data[i]))
                                   );
            i++;
          } // x
        }// y
      }// z

   } // end rho is a matrix
  }// parallel
}// end of Compute_rhoxyz_linear
//------------------------------------------------------------------------------


/**
 * Calculate three temporary sums in the new pressure formula \n
 * non-linear absorbing case, SSE2 version.
 * @param [out] RHO_Temp  - rhox_sgx + rhoy_sgy + rhoz_sgz
 * @param [out] BonA_Temp - BonA + rho ^2 / 2 rho0  + (rhox_sgx + rhoy_sgy + rhoz_sgz)
 * @param [out] Sum_du    - rho0* (duxdx + duydy + duzdz)
 */
void TKSpaceFirstOrder3DSolver::Calculate_SumRho_BonA_SumDu_SSE2(TRealMatrix & RHO_Temp,
                                                                  TRealMatrix & BonA_Temp,
                                                                  TRealMatrix & Sum_du)
{
  // step of 4
  const size_t TotalElementCount_4 = (RHO_Temp.GetTotalElementCount() >> 2) << 2;

  const float * rhox_data = Get_rhox().GetRawData();
  const float * rhoy_data = Get_rhoy().GetRawData();
  const float * rhoz_data = Get_rhoz().GetRawData();

  const float * dux_data = Get_duxdx().GetRawData();
  const float * duy_data = Get_duydy().GetRawData();
  const float * duz_data = Get_duzdz().GetRawData();


  // set BonA to be either scalar or a matrix
        float * BonA;
        size_t  BonA_shift;
  const bool    BonA_flag = Parameters->Get_BonA_scalar_flag();

  if (BonA_flag)
  {
    BonA = &Parameters->Get_BonA_scalar();
    BonA_shift = 0;
  }
  else
  {
    BonA = Get_BonA().GetRawData();
    BonA_shift = 1;
  }


  // set rho0A to be either scalar or a matrix
        float * rho0_data;
        size_t  rho0_shift;
  const bool    rho0_flag = Parameters->Get_rho0_scalar_flag();

  if (rho0_flag)
  {
    rho0_data = &Parameters->Get_rho0_scalar();
    rho0_shift = 0;
  }
  else
  {
    rho0_data = Get_rho0().GetRawData();
    rho0_shift = 1;
  }

  // compute loop
  #pragma  omp parallel
  {
    float * RHO_Temp_Data  = RHO_Temp.GetRawData();
    float * BonA_Temp_Data = BonA_Temp.GetRawData();
    float * SumDU_Temp_Data= Sum_du.GetRawData();


    const __m128 Two_SSE   = _mm_set1_ps(2.0f);
          __m128 BonA_SSE  = _mm_set1_ps(Parameters->Get_BonA_scalar());
          __m128 rho0_SSE  = _mm_set1_ps(Parameters->Get_rho0_scalar());


   #pragma omp for schedule (static) nowait
   for (size_t i = 0; i < TotalElementCount_4; i+=4)
   {
      if (!BonA_flag) BonA_SSE = _mm_load_ps(&BonA[i]);

      __m128 xmm1 = _mm_load_ps(&rhox_data[i]);
      __m128 xmm2 = _mm_load_ps(&rhoy_data[i]);
      __m128 xmm3 = _mm_load_ps(&rhoz_data[i]);

      if (!rho0_flag)  rho0_SSE = _mm_load_ps(&rho0_data[i]);

      __m128 rho_xyz_sq_SSE;
      __m128 rho_xyz_el_SSE;

      //  register const float rho_xyz_el = rhox_data[i] + rhoy_data[i] + rhoz_data[i];
      rho_xyz_el_SSE = _mm_add_ps(xmm1, xmm2);
      rho_xyz_el_SSE = _mm_add_ps(xmm3, rho_xyz_el_SSE);

      // RHO_Temp_Data[i]  = rho_xyz_el;
      _mm_stream_ps(&RHO_Temp_Data[i], rho_xyz_el_SSE);

      //  BonA_Temp_Data[i] =  ((BonA * (rho_xyz_el * rho_xyz_el)) / (2.0f * rho0_data[i])) + rho_xyz_el;
      rho_xyz_sq_SSE = _mm_mul_ps(rho_xyz_el_SSE, rho_xyz_el_SSE);// (rho_xyz_el * rho_xyz_el)

      xmm1           = _mm_mul_ps(rho_xyz_sq_SSE, BonA_SSE);      //((BonA * (rho_xyz_el * rho_xyz_el))
      xmm2           = _mm_mul_ps(Two_SSE, rho0_SSE);             // (2.0f * rho0_data[i])
      xmm3           = _mm_div_ps(xmm1, xmm2);                    // (BonA * (rho_xyz_el * rho_xyz_el)) /  (2.0f * rho0_data[i])

      xmm1           = _mm_add_ps(xmm3, rho_xyz_el_SSE);          // + rho_xyz_el

      _mm_stream_ps(&BonA_Temp_Data[i], xmm1);   //bypass cache

      xmm1       = _mm_load_ps(&dux_data[i]); //dudx
      xmm2       = _mm_load_ps(&duy_data[i]); //dudu
      xmm3       = _mm_load_ps(&duz_data[i]); //dudz

      __m128 xmm_acc = _mm_add_ps(xmm1, xmm2);
      xmm_acc    = _mm_add_ps(xmm_acc, xmm3);
      xmm_acc    = _mm_mul_ps(xmm_acc, rho0_SSE);

      _mm_stream_ps(&SumDU_Temp_Data[i],xmm_acc);

    // BonA_Temp_Data[i] =  ((BonA * (rho_xyz_el * rho_xyz_el)) / (2.0f * rho0_data[i])) + rho_xyz_el;
    }

    // non SSE code, in OpenMP only the last thread does this
    #ifdef _OPENMP
      if (omp_get_thread_num() == omp_get_num_threads() -1)
    #endif
    {
      for (size_t i = TotalElementCount_4; i < RHO_Temp.GetTotalElementCount() ; i++)
      {
        register const float rho_xyz_el = rhox_data[i] + rhoy_data[i] + rhoz_data[i];

        RHO_Temp_Data[i]   = rho_xyz_el;
        BonA_Temp_Data[i]  = ((BonA[i * BonA_shift] * (rho_xyz_el * rho_xyz_el)) / (2.0f * rho0_data[i* rho0_shift])) + rho_xyz_el;
        SumDU_Temp_Data[i] = rho0_data[i * rho0_shift] * (dux_data[i] + duy_data[i] + duz_data[i]);
      }
    }

  }// parallel
 } // end of Calculate_SumRho_BonA_SumDu_SSE2
//------------------------------------------------------------------------------



 /**
  * Calculate two temporary sums in the new pressure formula, linear absorbing case.
  * @param [out] Sum_rhoxyz  - rhox_sgx + rhoy_sgy + rhoz_sgz
  * @param [out] Sum_rho0_du - rho0* (duxdx + duydy + duzdz);
  */
void TKSpaceFirstOrder3DSolver::Calculate_SumRho_SumRhoDu(TRealMatrix & Sum_rhoxyz,
                                                           TRealMatrix & Sum_rho0_du)
 {
  const size_t TotalElementCount = Parameters->GetFullDimensionSizes().GetElementCount();

  #pragma  omp parallel
  {
    const float * rhox_data = Get_rhox().GetRawData();
    const float * rhoy_data = Get_rhoy().GetRawData();
    const float * rhoz_data = Get_rhoz().GetRawData();

    const float * dux_data = Get_duxdx().GetRawData();
    const float * duy_data = Get_duydy().GetRawData();
    const float * duz_data = Get_duzdz().GetRawData();

    const float * rho0_data = NULL;

    const float rho0_data_el = Parameters->Get_rho0_scalar();
    if (!Parameters->Get_rho0_scalar_flag())
    {
      rho0_data = Get_rho0().GetRawData();
    }


    float * Sum_rhoxyz_Data  = Sum_rhoxyz.GetRawData();
    float * Sum_rho0_du_Data = Sum_rho0_du.GetRawData();

    #pragma omp for schedule (static)
    for (size_t i = 0; i < TotalElementCount; i++)
    {
      Sum_rhoxyz_Data[i] = rhox_data[i] + rhoy_data[i] + rhoz_data[i];
    }


    if (Parameters->Get_rho0_scalar_flag())
    { // scalar
      #pragma omp for schedule (static) nowait
      for (size_t i = 0; i < TotalElementCount; i++)
      {
        Sum_rho0_du_Data[i] = rho0_data_el * (dux_data[i] + duy_data[i] + duz_data[i]);
      }
    }
    else
    { // matrix
      #pragma omp for schedule (static) nowait
      for (size_t i = 0; i < TotalElementCount; i++)
      {
        Sum_rho0_du_Data[i] = rho0_data[i] * (dux_data[i] + duy_data[i] + duz_data[i]);
      }
    }
  } // parallel
}// end of Calculate_SumRho_SumRhoDu
//------------------------------------------------------------------------------


 /**
  * Compute absorbing term with abosrb_nabla1 and absorb_nabla2, SSE2 version. \n
  * Calculate absorb_nabla1 .* fft_1 \n
  * Calculate absorb_nabla2 .* fft_2 \n
  *
  * @param [in,out] FFT_1
  * @param [in,out] FFT_2
  */
void TKSpaceFirstOrder3DSolver::Compute_Absorb_nabla1_2_SSE2(TFFTWComplexMatrix& FFT_1,
                                                              TFFTWComplexMatrix& FFT_2)
{
  const float * nabla1 = Get_absorb_nabla1().GetRawData();
  const float * nabla2 = Get_absorb_nabla2().GetRawData();

  const size_t TotalElementCount     = FFT_1.GetTotalElementCount();
  const size_t TotalElementCount_SSE = (FFT_1.GetTotalElementCount() >> 1) << 1;

  #pragma omp parallel
  {
    float * FFT_1_data  = FFT_1.GetRawData();
    float * FFT_2_data  = FFT_2.GetRawData();

    #pragma omp for schedule (static) nowait
    for (size_t i = 0; i < TotalElementCount_SSE; i+=2)
    {
       __m128 xmm_nabla1 = _mm_set_ps(nabla1[i+1], nabla1[i+1], nabla1[i], nabla1[i]);
       __m128 xmm_FFT_1  = _mm_load_ps(&FFT_1_data[2*i]);

              xmm_FFT_1  = _mm_mul_ps(xmm_nabla1, xmm_FFT_1);
                           _mm_store_ps(&FFT_1_data[2*i], xmm_FFT_1);
    }

    #pragma omp for schedule (static)
    for (size_t i = 0; i < TotalElementCount; i+=2)
    {
      __m128 xmm_nabla2 = _mm_set_ps(nabla2[i+1], nabla2[i+1], nabla2[i], nabla2[i]);
      __m128 xmm_FFT_2  = _mm_load_ps(&FFT_2_data[2*i]);

             xmm_FFT_2  = _mm_mul_ps(xmm_nabla2, xmm_FFT_2);
                          _mm_store_ps(&FFT_2_data[2*i], xmm_FFT_2);
      }

    //-- non SSE code --//
    #ifdef _OPENMP
      if (omp_get_thread_num() == omp_get_num_threads() -1)
    #endif
    {
      for (size_t i = TotalElementCount_SSE; i < TotalElementCount ; i++)
      {
        FFT_1_data[(i<<1)]   *= nabla1[i];
        FFT_1_data[(i<<1)+1] *= nabla1[i];

        FFT_2_data[(i<<1)]   *=  nabla2[i];
        FFT_2_data[(i<<1)+1] *=  nabla2[i];
      }
    }

  }// parallel
 } // end of Compute_Absorb_nabla1_2_SSE2
//------------------------------------------------------------------------------


 /**
  * Sum sub-terms to calculate new pressure, non-linear case.
  * @param [in] Absorb_tau_temp -
  * @param [in] Absorb_eta_temp - BonA + rho ^2 / 2 rho0  + (rhox_sgx + rhoy_sgy + rhoz_sgz)
  * @param [in] BonA_temp       - rho0* (duxdx + duydy + duzdz)
  */
void TKSpaceFirstOrder3DSolver::Sum_Subterms_nonlinear(TRealMatrix& Absorb_tau_temp,
                                                       TRealMatrix& Absorb_eta_temp,
                                                       TRealMatrix& BonA_temp)
{
  float * tau_data;
  float * eta_data;
  float * c2_data;

  size_t  c2_shift;
  size_t  tau_eta_shift;

  const float * Absorb_tau_data = Absorb_tau_temp.GetRawData();
  const float * Absorb_eta_data = Absorb_eta_temp.GetRawData();

  const size_t TotalElementCount = Get_p().GetTotalElementCount();
  const float  Divider = 1.0f / (float) TotalElementCount;

  // c2 scalar
  if (Parameters->Get_c0_scalar_flag())
  {
    c2_data = &Parameters->Get_c0_scalar();
    c2_shift = 0;
  }
  else
  {
    c2_data  = Get_c2().GetRawData();
    c2_shift = 1;
  }

  // eta and tau scalars
  if (Parameters->Get_c0_scalar_flag() && Parameters->Get_alpha_coeff_scallar_flag())
  {
    tau_data = &Parameters->Get_absorb_tau_scalar();
    eta_data = &Parameters->Get_absorb_eta_scalar();
    tau_eta_shift = 0;
  }
  else
  {
    tau_data = Get_absorb_tau().GetRawData();
    eta_data = Get_absorb_eta().GetRawData();
    tau_eta_shift = 1;
  }

  #pragma omp parallel
  {
    const float * BonA_data = BonA_temp.GetRawData();
    float * p_data  = Get_p().GetRawData();

    #pragma omp for schedule (static)
    for (size_t i = 0; i < TotalElementCount; i++)
    {
      p_data[i] = (c2_data[i * c2_shift]) *(
                   BonA_data[i] +
                   (Divider * ((Absorb_tau_data[i] * tau_data[i * tau_eta_shift]) -
                               (Absorb_eta_data[i] * eta_data[i * tau_eta_shift])
                              ))
                );
    }

  }// parallel
}// end of Sum_Subterms_nonlinear
//------------------------------------------------------------------------------


 /**
  * Sum sub-terms to calculate new pressure, linear case.
  * @param [in] Absorb_tau_temp - sub-term with absorb_tau
  * @param [in] Absorb_eta_temp - sub-term with absorb_eta
  * @param [in] Sum_rhoxyz      - rhox_sgx + rhoy_sgy + rhoz_sgz
  */
void TKSpaceFirstOrder3DSolver::Sum_Subterms_linear(TRealMatrix& Absorb_tau_temp,
                                                    TRealMatrix& Absorb_eta_temp,
                                                    TRealMatrix& Sum_rhoxyz)
{
  const float *  tau_data = NULL;
  const float *  eta_data = NULL;
  const float *  c2_data  = NULL;

  size_t c2_shift      = 0;
  size_t tau_eta_shift = 0;

  const float * Absorb_tau_data = Absorb_tau_temp.GetRawData();
  const float * Absorb_eta_data = Absorb_eta_temp.GetRawData();

  const size_t TotalElementCount = Parameters->GetFullDimensionSizes().GetElementCount();
  const float Divider = 1.0f / (float) TotalElementCount;

  // c2 scalar
  if (Parameters->Get_c0_scalar_flag())
  {
    c2_data = &Parameters->Get_c0_scalar();
    c2_shift = 0;
  }
  else
  {
    c2_data = Get_c2().GetRawData();
    c2_shift = 1;
  }

  // eta and tau scalars
  if (Parameters->Get_c0_scalar_flag() && Parameters->Get_alpha_coeff_scallar_flag())
  {
    tau_data = &Parameters->Get_absorb_tau_scalar();
    eta_data = &Parameters->Get_absorb_eta_scalar();
    tau_eta_shift = 0;
  }
  else
  {
    tau_data = Get_absorb_tau().GetRawData();
    eta_data = Get_absorb_eta().GetRawData();
    tau_eta_shift = 1;
  }

  #pragma omp parallel
  {
    const float * Sum_rhoxyz_Data = Sum_rhoxyz.GetRawData();
          float * p_data  = Get_p().GetRawData();

    #pragma omp for schedule (static)
    for (size_t i = 0; i < TotalElementCount; i++)
    {
      p_data[i] = (c2_data[i * c2_shift]) *
                    (Sum_rhoxyz_Data[i] +
                      (Divider * ((Absorb_tau_data[i] * tau_data[i * tau_eta_shift]) -
                                 (Absorb_eta_data[i] * eta_data[i * tau_eta_shift])
                      ))
                );
    }
  }// parallel
}// end of Sum_Subterms_linear
//------------------------------------------------------------------------------



/**
 * Sum sub-terms for new p, non-linear lossless case.
 *
 */
 void TKSpaceFirstOrder3DSolver::Sum_new_p_nonlinear_lossless()
 {
  #pragma omp parallel
  {
    const size_t TotalElementCount = Parameters->GetFullDimensionSizes().GetElementCount();
    float * p_data = Get_p().GetRawData();

    const float * rhox_data = Get_rhox().GetRawData();
    const float * rhoy_data = Get_rhoy().GetRawData();
    const float * rhoz_data = Get_rhoz().GetRawData();

    float * c2_data;
    size_t  c2_shift;

    if (Parameters->Get_c0_scalar_flag())
    {
      c2_data = &Parameters->Get_c0_scalar();
      c2_shift = 0;
    }
    else
    {
      c2_data = Get_c2().GetRawData();
      c2_shift = 1;
    }

    float * BonA_data;
    size_t  BonA_shift;

    if (Parameters->Get_BonA_scalar_flag())
    {
      BonA_data = &Parameters->Get_BonA_scalar();
      BonA_shift = 0;
    }
    else
    {
      BonA_data = Get_BonA().GetRawData();
      BonA_shift = 1;
    }

    float * rho0_data;
    size_t rho0_shift;

    if (Parameters->Get_rho0_scalar_flag())
    {
      rho0_data = &Parameters->Get_rho0_scalar();
      rho0_shift = 0;
    }
    else
    {
      rho0_data = Get_rho0().GetRawData();
      rho0_shift = 1;
    }


    #pragma omp for schedule (static)
    for (size_t i = 0; i < TotalElementCount; i++)
    {
      const float sum_rho = rhox_data[i] + rhoy_data[i] + rhoz_data[i];

      p_data[i] = c2_data[i * c2_shift] *
                    (sum_rho +
                              (BonA_data[i * BonA_shift] * (sum_rho * sum_rho) /
                              (2.0f * rho0_data[i * rho0_shift]))
                    );
    }
  }// parallel
 }// end of Sum_new_p_nonlinear_lossless
//------------------------------------------------------------------------------


 /**
  * Sum sub-terms for new p, linear lossless case.
  *
  */
 void TKSpaceFirstOrder3DSolver::Sum_new_p_linear_lossless()
{
  #pragma omp parallel
  {
    const float * rhox = Get_rhox().GetRawData();
    const float * rhoy = Get_rhoy().GetRawData();
    const float * rhoz = Get_rhoz().GetRawData();
          float * p_data = Get_p().GetRawData();
    const size_t  TotalElementCount = Parameters->GetFullDimensionSizes().GetElementCount();

    if (Parameters->Get_c0_scalar_flag())
    {
      const float c2_element = Parameters->Get_c0_scalar();

      #pragma omp for schedule (static)
      for (size_t i = 0; i < TotalElementCount; i++)
      {
        p_data[i] = c2_element * (rhox[i] + rhoy[i] + rhoz[i]);
      }
    }
    else
    {
      const float * c2_data = Get_c2().GetRawData();

      #pragma omp for schedule (static)
      for (size_t i = 0; i < TotalElementCount; i++)
      {
        p_data[i] = (c2_data[i]) * (rhox[i] + rhoy[i] + rhoz[i]);
      }
    }
  }// parallel
}// end of Sum_new_p_linear_lossless()
//------------------------------------------------------------------------------


/**
 * Compute new p for non-linear case.
 */
 void TKSpaceFirstOrder3DSolver::Compute_new_p_nonlinear()
{
  // rhox + rhoy + rhoz
  if (Parameters->Get_absorbing_flag())
  { // absorbing case

    TRealMatrix & Sum_rhoxyz = Get_Temp_1_RS3D();
    TRealMatrix & BonA_rho_rhoxyz = Get_Temp_2_RS3D();
    TRealMatrix & Sum_du = Get_Temp_3_RS3D();

    Calculate_SumRho_BonA_SumDu_SSE2(Sum_rhoxyz, BonA_rho_rhoxyz, Sum_du);

    // ifftn ( absorb_nabla1 * fftn (rho0 * (duxdx+duydy+duzdz))
    Get_FFT_X_temp().Compute_FFT_3D_R2C(Sum_du);
    Get_FFT_Y_temp().Compute_FFT_3D_R2C(Sum_rhoxyz);

    Compute_Absorb_nabla1_2_SSE2(Get_FFT_X_temp(), Get_FFT_Y_temp());

    TRealMatrix& Absorb_tau_temp = Sum_du;
    TRealMatrix& Absorb_eta_temp = Sum_rhoxyz;

    Get_FFT_X_temp().Compute_FFT_3D_C2R(Absorb_tau_temp);
    Get_FFT_Y_temp().Compute_FFT_3D_C2R(Absorb_eta_temp);

    Sum_Subterms_nonlinear(Absorb_tau_temp, Absorb_eta_temp, BonA_rho_rhoxyz);
  }
  else
  {
    Sum_new_p_nonlinear_lossless();
  }
}// end of Compute_new_p_nonlinear
//------------------------------------------------------------------------------


/**
 * Compute new p for linear case.
 */
 void TKSpaceFirstOrder3DSolver::Compute_new_p_linear()
 {
  // rhox + rhoy + rhoz
  if (Parameters->Get_absorbing_flag())
  { // absorbing case

    TRealMatrix& Sum_rhoxyz = Get_Temp_1_RS3D();
    TRealMatrix& Sum_rho0_du = Get_Temp_2_RS3D();

    Calculate_SumRho_SumRhoDu(Sum_rhoxyz, Sum_rho0_du);

    // ifftn ( absorb_nabla1 * fftn (rho0 * (duxdx+duydy+duzdz))

    Get_FFT_X_temp().Compute_FFT_3D_R2C(Sum_rho0_du);
    Get_FFT_Y_temp().Compute_FFT_3D_R2C(Sum_rhoxyz);

    Compute_Absorb_nabla1_2_SSE2(Get_FFT_X_temp(), Get_FFT_Y_temp());

    TRealMatrix& Absorb_tau_temp = Get_Temp_2_RS3D(); // override Sum_rho0_dx
    TRealMatrix& Absorb_eta_temp = Get_Temp_3_RS3D();

    Get_FFT_X_temp().Compute_FFT_3D_C2R(Absorb_tau_temp);
    Get_FFT_Y_temp().Compute_FFT_3D_C2R(Absorb_eta_temp);

    Sum_Subterms_linear(Absorb_tau_temp, Absorb_eta_temp, Sum_rhoxyz);
  }
  else
  {
    // lossless case
    Sum_new_p_linear_lossless();
  }
 }// end of Compute_new_p_linear
//------------------------------------------------------------------------------




 /*
  * Compute new values of ux_sgx, uy_sgy, uz_sgz.
  */
 void TKSpaceFirstOrder3DSolver::Compute_uxyz()
 {

  Compute_ddx_kappa_fft_p(Get_p(),
                          Get_FFT_X_temp(),Get_FFT_Y_temp(),Get_FFT_Z_temp(),
                          Get_kappa(),
                          Get_ddx_k_shift_pos(),Get_ddy_k_shift_pos(),Get_ddz_k_shift_pos());

   Get_FFT_X_temp().Compute_FFT_3D_C2R(Get_Temp_1_RS3D());
   Get_FFT_Y_temp().Compute_FFT_3D_C2R(Get_Temp_2_RS3D());
   Get_FFT_Z_temp().Compute_FFT_3D_C2R(Get_Temp_3_RS3D());

  #pragma omp parallel
  {
    if (Parameters->Get_rho0_scalar_flag())
    { // scalars
      if (Parameters->Get_nonuniform_grid_flag())
      {
        Get_ux_sgx().Compute_ux_sgx_normalize_scalar_nonuniform(Get_Temp_1_RS3D(),
                                                                Parameters->Get_rho0_sgx_scalar(),
                                                                Get_dxudxn_sgx(), Get_pml_x_sgx());
        Get_uy_sgy().Compute_uy_sgy_normalize_scalar_nonuniform(Get_Temp_2_RS3D(),
                                                                Parameters->Get_rho0_sgy_scalar(),
                                                                Get_dyudyn_sgy(), Get_pml_y_sgy());
        Get_uz_sgz().Compute_uz_sgz_normalize_scalar_nonuniform(Get_Temp_3_RS3D(),
                                                                Parameters->Get_rho0_sgz_scalar(),
                                                                Get_dzudzn_sgz(), Get_pml_z_sgz());
       }
      else
      {
        Get_ux_sgx().Compute_ux_sgx_normalize_scalar_uniform(Get_Temp_1_RS3D(),
                                                             Parameters->Get_rho0_sgx_scalar(),
                                                             Get_pml_x_sgx());
        Get_uy_sgy().Compute_uy_sgy_normalize_scalar_uniform(Get_Temp_2_RS3D(),
                                                             Parameters->Get_rho0_sgy_scalar(),
                                                             Get_pml_y_sgy());
        Get_uz_sgz().Compute_uz_sgz_normalize_scalar_uniform(Get_Temp_3_RS3D(),
                                                             Parameters->Get_rho0_sgz_scalar(),
                                                             Get_pml_z_sgz());
      }
    }
    else
    {// matrices
      Get_ux_sgx().Compute_ux_sgx_normalize(Get_Temp_1_RS3D(),
                                            Get_dt_rho0_sgx(),
                                            Get_pml_x_sgx());
      Get_uy_sgy().Compute_uy_sgy_normalize(Get_Temp_2_RS3D(),
                                            Get_dt_rho0_sgy(),
                                            Get_pml_y_sgy());
      Get_uz_sgz().Compute_uz_sgz_normalize(Get_Temp_3_RS3D(),
                                            Get_dt_rho0_sgz(),
                                            Get_pml_z_sgz());
    }
  } // parallel
}// end of Compute_uxyz()
//------------------------------------------------------------------------------

/**
 * Add u source to the particle velocity.
 */
void TKSpaceFirstOrder3DSolver::Add_u_source()
{
  const size_t t_index = Parameters->Get_t_index();

  if (Parameters->Get_ux_source_flag() > t_index)
  {
    Get_ux_sgx().Add_u_source(Get_ux_source_input(),
                              Get_u_source_index(),
                              t_index,
                              Parameters->Get_u_source_mode(),
                              Parameters->Get_u_source_many());
  }

  if (Parameters->Get_uy_source_flag() > t_index)
  {
    Get_uy_sgy().Add_u_source(Get_uy_source_input(),
                              Get_u_source_index(),
                              t_index,
                              Parameters->Get_u_source_mode(),
                              Parameters->Get_u_source_many());
  }

  if (Parameters->Get_uz_source_flag() > t_index)
  {
    Get_uz_sgz().Add_u_source(Get_uz_source_input(),
                              Get_u_source_index(),
                              t_index,
                              Parameters->Get_u_source_mode(),
                              Parameters->Get_u_source_many());
  }
}// end of Add_u_source
//------------------------------------------------------------------------------



 /**
  * Add in pressure source.
  *
  */
void TKSpaceFirstOrder3DSolver::Add_p_source()
{
  const size_t t_index = Parameters->Get_t_index();

  if (Parameters->Get_p_source_flag() > t_index)
  {
    float * rhox = Get_rhox().GetRawData();
    float * rhoy = Get_rhoy().GetRawData();
    float * rhoz = Get_rhoz().GetRawData();

    float  * p_source_input = Get_p_source_input().GetRawData();
    const size_t * p_source_index = Get_p_source_index().GetRawData();

    const bool   Is_p_source_many  = (Parameters->Get_p_source_many() != 0);
    const size_t Index2D           = (Is_p_source_many) ? t_index * Get_p_source_index().GetTotalElementCount() : t_index;
    const size_t p_source_size     = Get_p_source_index().GetTotalElementCount();

    // replacement
    if (Parameters->Get_p_source_mode() == 0)
    {
      #pragma omp parallel for if (p_source_size > 16384)
      for (size_t i = 0; i < p_source_size; i++)
      {
        const size_t SignalIndex = (Is_p_source_many) ? Index2D + i : Index2D;

        rhox[p_source_index[i]] = p_source_input[SignalIndex];
        rhoy[p_source_index[i]] = p_source_input[SignalIndex];
        rhoz[p_source_index[i]] = p_source_input[SignalIndex];

      }
    }
    // Addition
    else
    {
      #pragma omp parallel for if (p_source_size > 16384)
      for (size_t i = 0; i < p_source_size; i++)
      {
        const size_t SignalIndex = (Is_p_source_many) ? Index2D + i : Index2D;

        rhox[p_source_index[i]] += p_source_input[SignalIndex];
        rhoy[p_source_index[i]] += p_source_input[SignalIndex];
        rhoz[p_source_index[i]] += p_source_input[SignalIndex];
      }
    }// type of replacement
  }// if do at all
}// end of Add_p_source
//------------------------------------------------------------------------------


/**
 * Calculated shifted velocities.
 * \n
 * ux_shifted = real(ifft(bsxfun(\@times, x_shift_neg, fft(ux_sgx, [], 1)), [], 1)); \n
 * uy_shifted = real(ifft(bsxfun(\@times, y_shift_neg, fft(uy_sgy, [], 2)), [], 2)); \n
 * uz_shifted = real(ifft(bsxfun(\@times, z_shift_neg, fft(uz_sgz, [], 3)), [], 3)); \n
 */
void TKSpaceFirstOrder3DSolver::Calculate_shifted_velocity()
{
  const TFloatComplex * x_shift_neg_r  = (TFloatComplex *) Get_x_shift_neg_r().GetRawData();
  const TFloatComplex * y_shift_neg_r  = (TFloatComplex *) Get_y_shift_neg_r().GetRawData();
  const TFloatComplex * z_shift_neg_r  = (TFloatComplex *) Get_z_shift_neg_r().GetRawData();

        TFloatComplex * FFT_shift_temp = (TFloatComplex *) Get_FFT_shift_temp().GetRawData();


  // sizes of frequency spaces
  TDimensionSizes XShiftDims = Parameters->GetFullDimensionSizes();
                  XShiftDims.X = XShiftDims.X/2 + 1;

  TDimensionSizes YShiftDims = Parameters->GetFullDimensionSizes();
                  YShiftDims.Y = YShiftDims.Y/2 + 1;

  TDimensionSizes ZShiftDims = Parameters->GetFullDimensionSizes();
                  ZShiftDims.Z = ZShiftDims.Z/2 + 1;

  // normalization constants for FFTs
  const float DividerX = 1.0f / (float) Parameters->GetFullDimensionSizes().X;
  const float DividerY = 1.0f / (float) Parameters->GetFullDimensionSizes().Y;
  const float DividerZ = 1.0f / (float) Parameters->GetFullDimensionSizes().Z;

  //-------------------------- ux_shifted --------------------------------------
  Get_FFT_shift_temp().Compute_FFT_1DX_R2C(Get_ux_sgx());

  #pragma omp parallel for schedule (static)
  for (size_t z = 0; z < XShiftDims.Z; z++)
  {
    register size_t i = z *  XShiftDims.Y * XShiftDims.X;
    for (size_t y = 0; y < XShiftDims.Y; y++)
    {
      for (size_t x = 0; x < XShiftDims.X; x++)
      {
        TFloatComplex Temp;

        Temp.real = ((FFT_shift_temp[i].real * x_shift_neg_r[x].real) -
                     (FFT_shift_temp[i].imag * x_shift_neg_r[x].imag)
                    ) * DividerX;


        Temp.imag = ((FFT_shift_temp[i].imag * x_shift_neg_r[x].real) +
                     (FFT_shift_temp[i].real * x_shift_neg_r[x].imag)
                    ) * DividerX;

        FFT_shift_temp[i] = Temp;

        i++;
      } // x
    } // y
  }//z*/
  Get_FFT_shift_temp().Compute_FFT_1DX_C2R(Get_ux_shifted());


  //-------------------------- uy_shifted --------------------------------------
  Get_FFT_shift_temp().Compute_FFT_1DY_R2C(Get_uy_sgy());

  #pragma omp parallel for schedule (static)
  for (size_t z = 0; z < YShiftDims.Z; z++)
  {
    register size_t i = z *  YShiftDims.Y * YShiftDims.X;
    for (size_t y = 0; y < YShiftDims.Y; y++)
    {
      for (size_t x = 0; x < YShiftDims.X; x++)
      {
        TFloatComplex Temp;

        Temp.real = ((FFT_shift_temp[i].real * y_shift_neg_r[y].real) -
                     (FFT_shift_temp[i].imag * y_shift_neg_r[y].imag)) *
                      DividerY;


        Temp.imag = ((FFT_shift_temp[i].imag * y_shift_neg_r[y].real) +
                     (FFT_shift_temp[i].real * y_shift_neg_r[y].imag)
                    ) * DividerY;

        FFT_shift_temp[i] = Temp;

        i++;
      } // x
    } // y
  }//z
  Get_FFT_shift_temp().Compute_FFT_1DY_C2R(Get_uy_shifted());


  //-------------------------- uz_shifted --------------------------------------
  Get_FFT_shift_temp().Compute_FFT_1DZ_R2C(Get_uz_sgz());
  #pragma omp parallel for schedule (static)
  for (size_t z = 0; z < ZShiftDims.Z; z++)
  {
    register size_t i = z *  ZShiftDims.Y * ZShiftDims.X;
    for (size_t y = 0; y < ZShiftDims.Y; y++)
    {
      for (size_t x = 0; x < ZShiftDims.X; x++)
      {
        TFloatComplex Temp;

        Temp.real = ((FFT_shift_temp[i].real * z_shift_neg_r[z].real) -
                     (FFT_shift_temp[i].imag * z_shift_neg_r[z].imag)) *
                      DividerZ;


        Temp.imag = ((FFT_shift_temp[i].imag * z_shift_neg_r[z].real) +
                     (FFT_shift_temp[i].real * z_shift_neg_r[z].imag)
                    ) * DividerZ;

        FFT_shift_temp[i] = Temp;

        i++;
      } // x
    } // y
  }//z
  Get_FFT_shift_temp().Compute_FFT_1DZ_C2R(Get_uz_shifted());
}// end of Calculate_shifted_velocity()
//------------------------------------------------------------------------------


/**
 * Compute the main time loop of KSpaceFirstOrder3D.
 *
 */
void TKSpaceFirstOrder3DSolver::Compute_MainLoop()
{

  ActPercent = 0;
  // set ActPercent to correspond the t_index after recovery
  if (Parameters->Get_t_index() > 0)
  {
    ActPercent = (100 * Parameters->Get_t_index()) / Parameters->Get_Nt();
  }

  PrintOtputHeader();

  IterationTime.Start();

  while (Parameters->Get_t_index() < Parameters->Get_Nt() && (!IsTimeToCheckpoint()))
  {
    const size_t t_index = Parameters->Get_t_index();

    Compute_uxyz();
    // add in the velocity u source term
    Add_u_source();

    // add in the transducer source term (t = t1) to ux
    if (Parameters->Get_transducer_source_flag() > t_index)
    {
     Get_ux_sgx().AddTransducerSource(Get_u_source_index(), Get_delay_mask(), Get_transducer_source_input());
    }

    Compute_duxyz();

    if (Parameters->Get_nonlinear_flag())
    {
      Compute_rhoxyz_nonlinear();
    }
    else
    {
      Compute_rhoxyz_linear();
    }


     // add in the source pressure term
     Add_p_source();

    if (Parameters->Get_nonlinear_flag())
    {
      Compute_new_p_nonlinear();
    }
    else
    {
      Compute_new_p_linear();
    }

    // calculate initial pressure
    if ((t_index == 0) && (Parameters->Get_p0_source_flag() == 1)) Calculate_p0_source();

    StoreSensorData();
    PrintStatisitcs();
    Parameters->Increment_t_index();

  }// time loop
}// end of Compute_Main_Loop()
//------------------------------------------------------------------------------


/**
 * Print progress statistics.
 *
 */
void TKSpaceFirstOrder3DSolver::PrintStatisitcs()
{
  const size_t Nt =  Parameters->Get_Nt();
  const size_t t_index = Parameters->Get_t_index();


  if (t_index > ((ActPercent * Nt) / 100) )
  {
    ActPercent += Parameters->GetVerboseInterval();

    IterationTime.Stop();

    const double ElTime = IterationTime.GetElapsedTime();
    const double ElTimeWithLegs = IterationTime.GetElapsedTime() + SimulationTime.GetCumulatedElapsedTimeOverPreviousLegs();
    const double ToGo   = ((ElTimeWithLegs / (double) (t_index + 1)) *  double(Nt)) - ElTimeWithLegs;

    struct tm *current;
    time_t now;
    time(&now);
    now += ToGo;
    current = localtime(&now);

    fprintf(stdout, "%5li%c      %9.3fs      %9.3fs      %02i/%02i/%02i %02i:%02i:%02i\n",
            (100 * t_index) / Nt ,'%',
            ElTime, ToGo,
            current->tm_mday, current->tm_mon + 1, current->tm_year - 100,
            current->tm_hour, current->tm_min, current->tm_sec
            );

    fflush(stdout);
  }
}// end of KSpaceFirstOrder3DSolver
//------------------------------------------------------------------------------

/**
 * Print the header of the progress statistics.
 */
void TKSpaceFirstOrder3DSolver::PrintOtputHeader()
{
  fprintf(stdout,"-------------------------------------------------------------\n");
  fprintf(stdout,"....................... Simulation ..........................\n");
  fprintf(stdout,"Progress...ElapsedTime........TimeToGo......TimeOfTermination\n");
}// end of PrintOtputHeader
//------------------------------------------------------------------------------

/**
 * Is time to checkpoint
 * @return true if it is necessary to stop to checkpoint?
 */
bool TKSpaceFirstOrder3DSolver::IsTimeToCheckpoint()
{
  if (!Parameters->IsCheckpointEnabled()) return false;

  TotalTime.Stop();

  return (TotalTime.GetElapsedTime() > (float) Parameters->GetCheckpointInterval());

}// end of IsTimeToCheckpoint
//------------------------------------------------------------------------------


/**
 * Post processing the quantities, closing the output streams and storing the
 * sensor mask.
 *
 */
void TKSpaceFirstOrder3DSolver::PostPorcessing()
{
  if (Parameters->IsStore_p_final())
  {
    Get_p().WriteDataToHDF5File(Parameters->HDF5_OutputFile, p_final_Name, Parameters->GetCompressionLevel());
  }// p_final

  if (Parameters->IsStore_u_final())
  {
    Get_ux_sgx().WriteDataToHDF5File(Parameters->HDF5_OutputFile, ux_final_Name, Parameters->GetCompressionLevel());
    Get_uy_sgy().WriteDataToHDF5File(Parameters->HDF5_OutputFile, uy_final_Name, Parameters->GetCompressionLevel());
    Get_uz_sgz().WriteDataToHDF5File(Parameters->HDF5_OutputFile, uz_final_Name, Parameters->GetCompressionLevel());
  }// u_final

  // Apply post-processing and close
  OutputStreamContainer.PostProcessStreams();
  OutputStreamContainer.CloseStreams();


  // store sensor mask if wanted
  if (Parameters->IsCopySensorMask())
  {
    if (Parameters->Get_sensor_mask_type() == TParameters::smt_index)
    {
      Get_sensor_mask_index().RecomputeIndicesToMatlab();
      Get_sensor_mask_index().WriteDataToHDF5File(Parameters->HDF5_OutputFile,sensor_mask_index_Name,
                                                  Parameters->GetCompressionLevel());
    }
    if (Parameters->Get_sensor_mask_type() == TParameters::smt_corners)
    {
      Get_sensor_mask_corners().RecomputeIndicesToMatlab();
      Get_sensor_mask_corners().WriteDataToHDF5File(Parameters->HDF5_OutputFile,sensor_mask_corners_Name,
                                                    Parameters->GetCompressionLevel());
    }
  }
}// end of PostPorcessing
//------------------------------------------------------------------------------


/**
 * Store sensor data.
 *
 */
void TKSpaceFirstOrder3DSolver::StoreSensorData()
{
  // Unless the time for sampling has come, exit
  if (Parameters->Get_t_index() >= Parameters->GetStartTimeIndex())
  {
    if (Parameters->IsStore_u_non_staggered_raw())
    {
      Calculate_shifted_velocity();
    }
    OutputStreamContainer.SampleStreams();
  }
}// end of StoreData
//------------------------------------------------------------------------------

/**
 * Save checkpoint data into the checkpoint file, flush aggregated outputs into
 * the output file.
 */
void TKSpaceFirstOrder3DSolver::SaveCheckpointData()
{
  // export FFTW wisdom
  TFFTWComplexMatrix::ExportWisdom();


  // Create Checkpoint file
  THDF5_File & HDF5_CheckpointFile = Parameters->HDF5_CheckpointFile;
  // if it happens and the file is opened (from the recovery, close it)
  if (HDF5_CheckpointFile.IsOpened()) HDF5_CheckpointFile.Close();
  // Create the new file (overwrite the old one)
  HDF5_CheckpointFile.Create(Parameters->GetCheckpointFileName().c_str());


  //--------------------- Store Matrices ------------------------------//

  // Store all necessary matrices in Checkpoint file
  MatrixContainer.StoreDataIntoCheckpointHDF5File(HDF5_CheckpointFile);

  // Write t_index
  HDF5_CheckpointFile.WriteScalarValue(HDF5_CheckpointFile.GetRootGroup(),
                                       t_index_Name,
                                       Parameters->Get_t_index());
  // store basic dimension sizes (Nx, Ny, Nz) - Nt is not necessary
  HDF5_CheckpointFile.WriteScalarValue(HDF5_CheckpointFile.GetRootGroup(),
                                       Nx_Name,
                                       Parameters->GetFullDimensionSizes().X);
  HDF5_CheckpointFile.WriteScalarValue(HDF5_CheckpointFile.GetRootGroup(),
                                       Ny_Name,
                                       Parameters->GetFullDimensionSizes().Y);
  HDF5_CheckpointFile.WriteScalarValue(HDF5_CheckpointFile.GetRootGroup(),
                                       Nz_Name,
                                       Parameters->GetFullDimensionSizes().Z);


  // Write checkpoint file header
  THDF5_FileHeader CheckpointFileHeader = Parameters->HDF5_FileHeader;

  CheckpointFileHeader.SetFileType(THDF5_FileHeader::hdf5_ft_checkpoint);
  CheckpointFileHeader.SetCodeName(GetCodeName());
  CheckpointFileHeader.SetActualCreationTime();

  CheckpointFileHeader.WriteHeaderToCheckpointFile(HDF5_CheckpointFile);

  // Close the checkpoint file
  HDF5_CheckpointFile.Close();

  // checkpoint only if necessary (t_index > start_index) - here we're at  step + 1
  if (Parameters->Get_t_index() > Parameters->GetStartTimeIndex())
  {
    OutputStreamContainer.CheckpointStreams();
  }
  OutputStreamContainer.CloseStreams();
}// end of SaveCheckpointData()
//------------------------------------------------------------------------------


/**
 * Write statistics and the header into the output file.
 */
void TKSpaceFirstOrder3DSolver::WriteOutputDataInfo()
{
  // write t_index into the output file
  Parameters->HDF5_OutputFile.WriteScalarValue(Parameters->HDF5_OutputFile.GetRootGroup(),
                                               t_index_Name,
                                               Parameters->Get_t_index());

  // Write scalars
  Parameters->SaveScalarsToHDF5File(Parameters->HDF5_OutputFile);
  THDF5_FileHeader & HDF5_FileHeader = Parameters->HDF5_FileHeader;

  // Write File header

  HDF5_FileHeader.SetCodeName(GetCodeName());
  HDF5_FileHeader.SetMajorFileVersion();
  HDF5_FileHeader.SetMinorFileVersion();
  HDF5_FileHeader.SetActualCreationTime();
  HDF5_FileHeader.SetFileType(THDF5_FileHeader::hdf5_ft_output);
  HDF5_FileHeader.SetHostName();

  HDF5_FileHeader.SetMemoryConsumption(ShowMemoryUsageInMB());

  // Stop total timer here
  TotalTime.Stop();
  HDF5_FileHeader.SetExecutionTimes(GetCumulatedTotalTime(),
                                    GetCumulatedDataLoadTime(),
                                    GetCumulatedPreProcessingTime(),
                                    GetCumulatedSimulationTime(),
                                    GetCumulatedPostProcessingTime());

  HDF5_FileHeader.SetNumberOfCores();
  HDF5_FileHeader.WriteHeaderToOutputFile(Parameters->HDF5_OutputFile);
}// end of WriteOutputDataInfo
//------------------------------------------------------------------------------


/**
 *
 * Restore cumulated elapsed time form Output file (header stored in TParameters)
 * Open the header, read this and store into TTimeMeasure classes.
 */
void TKSpaceFirstOrder3DSolver::RestoreCumulatedElapsedFromOutputFile()
{
  double ElapsedTotalTime, ElapsedDataLoadTime, ElapsedPreProcessingTime,
         ElapsedSimulationTime, ElapsedPostProcessingTime;

  // Get execution times stored in the output file header
  Parameters->HDF5_FileHeader.GetExecutionTimes(ElapsedTotalTime,
                                                ElapsedDataLoadTime,
                                                ElapsedPreProcessingTime,
                                                ElapsedSimulationTime,
                                                ElapsedPostProcessingTime);

  TotalTime.SetCumulatedElapsedTimeOverPreviousLegs(ElapsedTotalTime);
  DataLoadTime.SetCumulatedElapsedTimeOverPreviousLegs(ElapsedDataLoadTime);
  PreProcessingTime.SetCumulatedElapsedTimeOverPreviousLegs(ElapsedPreProcessingTime);
  SimulationTime.SetCumulatedElapsedTimeOverPreviousLegs(ElapsedSimulationTime);
  PostProcessingTime.SetCumulatedElapsedTimeOverPreviousLegs(ElapsedPostProcessingTime);

}// end of RestoreCumulatedElapsedFromOutputFile
//------------------------------------------------------------------------------


/**
 * Check the output file has the correct format and version.
 * @throw ios::failure if an error happens
 */
void TKSpaceFirstOrder3DSolver::CheckOutputFile()
{
  // The header has already been read
  THDF5_FileHeader & OutputFileHeader = Parameters->HDF5_FileHeader;
  THDF5_File       & OutputFile       = Parameters->HDF5_OutputFile;

  // test file type
  if (OutputFileHeader.GetFileType() != THDF5_FileHeader::hdf5_ft_output)
  {
    char ErrorMessage[256] = "";
    sprintf(ErrorMessage,
            KSpaceFirstOrder3DSolver_ERR_FMT_IncorrectOutputFileFormat,
            Parameters->GetOutputFileName().c_str());
    throw ios::failure(ErrorMessage);
  }

  // test file major version
  if (!OutputFileHeader.CheckMajorFileVersion())
  {
    char ErrorMessage[256] = "";
    sprintf(ErrorMessage,
            Parameters_ERR_FMT_IncorrectMajorHDF5FileVersion,
            Parameters->GetCheckpointFileName().c_str(),
            OutputFileHeader.GetCurrentHDF5_MajorVersion().c_str());
    throw ios::failure(ErrorMessage);
  }

  // test file minor version
  if (!OutputFileHeader.CheckMinorFileVersion())
  {
    char ErrorMessage[256] = "";
    sprintf(ErrorMessage,
            Parameters_ERR_FMT_IncorrectMinorHDF5FileVersion,
            Parameters->GetCheckpointFileName().c_str(),
            OutputFileHeader.GetCurrentHDF5_MinorVersion().c_str());
    throw ios::failure(ErrorMessage);
  }


  // Check dimension sizes
  TDimensionSizes OutputDimSizes;
  OutputFile.ReadScalarValue(OutputFile.GetRootGroup(),
                             Nx_Name,
                             OutputDimSizes.X);

  OutputFile.ReadScalarValue(OutputFile.GetRootGroup(),
                             Ny_Name,
                             OutputDimSizes.Y);

  OutputFile.ReadScalarValue(OutputFile.GetRootGroup(),
                             Nz_Name,
                             OutputDimSizes.Z);

 if (Parameters->GetFullDimensionSizes() != OutputDimSizes)
 {
    char ErrorMessage[256] = "";
    sprintf(ErrorMessage,
            KSpaceFirstOrder3DSolver_ERR_FMT_OutputDimensionsDoNotMatch,
            OutputDimSizes.X,
            OutputDimSizes.Y,
            OutputDimSizes.Z,
            Parameters->GetFullDimensionSizes().X,
            Parameters->GetFullDimensionSizes().Y,
            Parameters->GetFullDimensionSizes().Z);

   throw ios::failure(ErrorMessage);
 }
}// end of CheckOutputFile
//------------------------------------------------------------------------------


/**
 * Check the file type and the version of the checkpoint file.
 * @throw ios::failure if an error happens
 *
 */
void TKSpaceFirstOrder3DSolver::CheckCheckpointFile()
{
  // read the header and check the file version
  THDF5_FileHeader CheckpointFileHeader;
  THDF5_File &     HDF5_CheckpointFile = Parameters->HDF5_CheckpointFile;

  CheckpointFileHeader.ReadHeaderFromCheckpointFile(HDF5_CheckpointFile);

  // test file type
  if (CheckpointFileHeader.GetFileType() != THDF5_FileHeader::hdf5_ft_checkpoint)
  {
    char ErrorMessage[256] = "";
    sprintf(ErrorMessage,
            KSpaceFirstOrder3DSolver_ERR_FMT_IncorrectCheckpointFileFormat,
            Parameters->GetCheckpointFileName().c_str());
    throw ios::failure(ErrorMessage);
  }

  // test file major version
  if (!CheckpointFileHeader.CheckMajorFileVersion())
  {
    char ErrorMessage[256] = "";
    sprintf(ErrorMessage,
            Parameters_ERR_FMT_IncorrectMajorHDF5FileVersion,
            Parameters->GetCheckpointFileName().c_str(),
            CheckpointFileHeader.GetCurrentHDF5_MajorVersion().c_str());
    throw ios::failure(ErrorMessage);
  }

  // test file minor version
  if (!CheckpointFileHeader.CheckMinorFileVersion())
  {
    char ErrorMessage[256] = "";
    sprintf(ErrorMessage,
            Parameters_ERR_FMT_IncorrectMinorHDF5FileVersion,
            Parameters->GetCheckpointFileName().c_str(),
            CheckpointFileHeader.GetCurrentHDF5_MinorVersion().c_str());
    throw ios::failure(ErrorMessage);
  }


  // Check dimension sizes
  TDimensionSizes CheckpointDimSizes;
  HDF5_CheckpointFile.ReadScalarValue(HDF5_CheckpointFile.GetRootGroup(),
                                      Nx_Name,
                                      CheckpointDimSizes.X);

  HDF5_CheckpointFile.ReadScalarValue(HDF5_CheckpointFile.GetRootGroup(),
                                      Ny_Name,
                                      CheckpointDimSizes.Y);

  HDF5_CheckpointFile.ReadScalarValue(HDF5_CheckpointFile.GetRootGroup(),
                                      Nz_Name,
                                      CheckpointDimSizes.Z);

 if (Parameters->GetFullDimensionSizes() != CheckpointDimSizes)
 {
    char ErrorMessage[256] = "";
    sprintf(ErrorMessage,
            KSpaceFirstOrder3DSolver_ERR_FMT_CheckpointDimensionsDoNotMatch,
            CheckpointDimSizes.X,
            CheckpointDimSizes.Y,
            CheckpointDimSizes.Z,
            Parameters->GetFullDimensionSizes().X,
            Parameters->GetFullDimensionSizes().Y,
            Parameters->GetFullDimensionSizes().Z);

   throw ios::failure(ErrorMessage);
 }
}// end of CheckCheckpointFile
//------------------------------------------------------------------------------

//----------------------------------------------------------------------------//
//                            Private methods                                 //
//----------------------------------------------------------------------------//
