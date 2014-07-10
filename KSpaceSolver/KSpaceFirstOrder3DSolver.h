/**
 * @file        KSpaceFirstOrder3DSolver.h
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au
 * @brief       The header file containing the main class of the project
 *              responsible for the entire simulation.
 * @version     kspaceFirstOrder3D 2.14
 * @date        12 July     2012, 10:27   (created)\n
 *              07 July     2014, 15:40   (revised)
 *
 * @section License
 * This file is part of the C++ extension of the k-Wave Toolbox (http://www.k-wave.org).\n
 * Copyright (C) 2012 Jiri Jaros and Bradley Treeby.
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

#ifndef TKSPACE3DSOLVER_H
#define	TKSPACE3DSOLVER_H

#include <fftw3.h>

#include <Parameters/Parameters.h>

#include <MatrixClasses/MatrixContainer.h>
#include <MatrixClasses/RealMatrix.h>
#include <MatrixClasses/ComplexMatrix.h>
#include <MatrixClasses/LongMatrix.h>
#include <MatrixClasses/OutputHDF5Stream.h>
#include <MatrixClasses/UXYZ_SGXYZMatrix.h>
#include <MatrixClasses/FFTWComplexMatrix.h>

#include <Utils/TimeMeasure.h>

using namespace std;



/**
 * \class TKSpaceFirstOrder3DSolver
 * \brief Class responsible for running the k-space first order 3D method.
 *
 */
class TKSpaceFirstOrder3DSolver
{
  public:
    /// Constructor
    TKSpaceFirstOrder3DSolver();

    /// Destructor
    virtual ~TKSpaceFirstOrder3DSolver();


    /// Memory allocation
    virtual void AllocateMemory();
    /// Memory deallocation
    virtual void FreeMemory();

    /// Load simulation data from the input file
    virtual void LoadInputData();

    /// Compute the 3D kspace first order simulation
    virtual void Compute();

    /// Print parameters of the simulation
    virtual void PrintParametersOfSimulation(FILE * file);

    /// Get memory usage in MB
    virtual size_t ShowMemoryUsageInMB();

    /// Get code name
    string GetCodeName() {return "kspaceFirstOrder3D-OMP v1.1"; };

    /// Print the code name and license
    void   PrintFullNameCodeAndLicense(FILE * file);

    /// Get total simulation time
    double GetTotalTime()          const { return TotalTime.GetElapsedTime();         };

    /// Get pre-processing time
    double GetPreProcessingTime()  const { return PreProcessingTime.GetElapsedTime(); };

    /// Get data load time
    double GetDataLoadTime()       const { return DataLoadTime.GetElapsedTime();      };

    /// Get simulation time (time loop)
    double GetSimulationTime()     const { return SimulationTime.GetElapsedTime();    };

    /// Get post-processing time
    double GetPostProcessingTime() const { return PostProcessingTime.GetElapsedTime();};

    /// Get total simulation time cumulated over all legs
    double GetCumulatedTotalTime()          const { return TotalTime.GetComulatedElaspedTimeOverAllLegs();          };
    /// Get pre-processing time cumulated over all legs
    double GetCumulatedPreProcessingTime()  const { return PreProcessingTime.GetComulatedElaspedTimeOverAllLegs();  };
    /// Get data load time cumulated over all legs
    double GetCumulatedDataLoadTime()       const { return DataLoadTime.GetComulatedElaspedTimeOverAllLegs();       };
    /// Get simulation time (time loop) cumulated over all legs
    double GetCumulatedSimulationTime()     const { return SimulationTime.GetComulatedElaspedTimeOverAllLegs();     };
    /// Get post-processing time cumulated over all legs
    double GetCumulatedPostProcessingTime() const { return PostProcessingTime.GetComulatedElaspedTimeOverAllLegs(); };

  protected:

    /// Copy constructor not allowed for public
    TKSpaceFirstOrder3DSolver(const TKSpaceFirstOrder3DSolver& src);
    /// operator = not allowed for public
    TKSpaceFirstOrder3DSolver& operator = (const TKSpaceFirstOrder3DSolver& src);

    /// Initialize FFT plans
    void InitializeFFTWPlans();

    /// Compute pre-processing phase
    void PreProcessingPhase();

    /// Compute the main time loop of the kspaceFirstOrder3D
    void Compute_MainLoop();

    /// Post processing, and closing the output streams
    void PostPorcessing();

    /// Store sensor data
    void StoreSensorData();
    /// Store intensity data
    void StoreIntensityData();

    /// Save checkpoint data
    void SaveCheckpointData();

    /// Write statistics and header into the output file
    void WriteOutputDataInfo();


    /// compute new values of for ux_sgx, uy_sgy, uz_sgz
    void Compute_uxyz();

    /// Compute new values of for duxdx, duydy, dzdz
    void Compute_duxyz();


    /// Compute new values of rhox, rhoy, rhoz for non-linear case
    void Compute_rhoxyz_nonlinear();
    /// Compute new values of rhox, rhoy, rhoz for linear case
    void Compute_rhoxyz_linear();

    /// Add u source to the particle velocity.
    void Add_u_source();

    /// Add in pressure source.
    void Add_p_source();


    /// Generate kappa matrix for non-absorbing media
    void Generate_kappa();

    /// Generate kappa matrix, absorb_nabla1, absorb_nabla2 for absorbing media
    void Generate_kappa_absorb_nabla1_absorb_nabla2();

    /// Generate absorb_tau, absorb_eta for heterogenous media
    void Generate_absorb_tau_absorb_eta_matrix();

    /// Calculate dt ./ rho0 for non-uniform grids
    void Caclucalte_dt_rho0_non_uniform();

    /// Calculate p0_source
    void Calculate_p0_source();

    /// Calculate c^2
    void Compute_c2();

    /// Compute part of the new velocity - gradient in p
    void Compute_ddx_kappa_fft_p(TRealMatrix& X_Matrix,
                                 TFFTWComplexMatrix& FFT_X,
                                 TFFTWComplexMatrix& FFT_Y,
                                 TFFTWComplexMatrix& FFT_Z,
                                 TRealMatrix& kappa,
                                 TComplexMatrix& ddx,
                                 TComplexMatrix& ddy,
                                 TComplexMatrix& ddz);

    /// Calculate new p, non-linear case
    void Compute_new_p_nonlinear();
    /// Calculate new p linear-case, absorbing
    void Compute_new_p_linear();


    /// Calculate three temporary sums in the new pressure formula, non-linear absorbing case, SSE2 version
    void Calculate_SumRho_BonA_SumDu_SSE2(TRealMatrix& RHO_Temp,
                                          TRealMatrix& BonA_Temp,
                                          TRealMatrix& Sum_du);
    /// Calculate two temporary sums in the new pressure formula, linear absorbing case
    void Calculate_SumRho_SumRhoDu(TRealMatrix& Sum_rhoxyz,
                                   TRealMatrix& Sum_rho0_du);

    /// Compute absorbing term with abosrb_nabla1 and absorb_nabla2, SSE2 version
    void Compute_Absorb_nabla1_2_SSE2(TFFTWComplexMatrix& FFT_1,
                                      TFFTWComplexMatrix& FFT_2);

    /// Sum sub-terms to calculate new pressure, non-linear case
    void Sum_Subterms_nonlinear(TRealMatrix& Absorb_tau_temp,
                                TRealMatrix& Absorb_eta_temp,
                                TRealMatrix& BonA_temp);

    /// Sum sub-terms to calculate new pressure, linear case
    void Sum_Subterms_linear(TRealMatrix& Absorb_tau_temp,
                             TRealMatrix& Absorb_eta_temp,
                             TRealMatrix& Sum_rhoxyz);

    /// Sum sub-terms for new p, linear lossless case
    void Sum_new_p_nonlinear_lossless();
    ///  Sum sub-terms for new p, linear lossless case
    void Sum_new_p_linear_lossless();


    /// Print progress statistics
    void PrintStatisitcs();

    /// Print the header of the progress statistics
    void PrintOtputHeader();

    /// Is time to checkpoint (save actual state on disk)
    bool IsTimeToCheckpoint();

    /// Was the loop interrupted to checkpoint
    bool IsCheckpointInterruption() const
    {
      return (Parameters->Get_t_index() != Parameters->Get_Nt());
    };

    /// Reads the header of the output file and sets the cumulative elapsed time from the first log
    void RestoreCumulatedElapsedFromOutputFile(THDF5_File& HDF5_OutputFile);


     //------------------------- Get matrices --------------------------------//
    /// Get the kappa matrix from the container
    TRealMatrix& Get_kappa()
    {
      return MatrixContainer.GetRealMatrix(kappa);
    };
    /// Get the c^2 matrix from the container
    TRealMatrix& Get_c2()
    {
      return MatrixContainer.GetRealMatrix(c2);
    };

    /// Get the p matrix from the container
    TRealMatrix& Get_p()
    {
      return MatrixContainer.GetRealMatrix(p);
    };

    /// Get the ux_sgx matrix from the container
    Tuxyz_sgxyzMatrix& Get_ux_sgx()
    {
      return MatrixContainer.GetUxyz_sgxyzMatrix(ux_sgx);
    };
    /// Get the uy_sgy matrix from the container
    Tuxyz_sgxyzMatrix& Get_uy_sgy()
    {
      return MatrixContainer.GetUxyz_sgxyzMatrix(uy_sgy);
    };
    /// Get the uz_sgz matrix from the container
    Tuxyz_sgxyzMatrix& Get_uz_sgz()
    {
      return MatrixContainer.GetUxyz_sgxyzMatrix(uz_sgz);
    };

    /// Get the duxdx matrix from the container
    TRealMatrix& Get_duxdx()
    {
      return MatrixContainer.GetRealMatrix(duxdx);
    };
    /// Get the duydy matrix from the container
    TRealMatrix& Get_duydy()
    {
      return MatrixContainer.GetRealMatrix(duydy);
    };
    /// Get the duzdz matrix from the container
    TRealMatrix& Get_duzdz()
    {
      return MatrixContainer.GetRealMatrix(duzdz);
    };

    /// Get the dt.*rho0_sgx matrix from the container
    TRealMatrix& Get_dt_rho0_sgx()
    {
      return MatrixContainer.GetRealMatrix(dt_rho0_sgx);
    };
    /// Get the dt.*rho0_sgy matrix from the container
    TRealMatrix& Get_dt_rho0_sgy()
    {
      return MatrixContainer.GetRealMatrix(dt_rho0_sgy);
    };
    /// Get the dt.*rho0_sgz matrix from the container
    TRealMatrix& Get_dt_rho0_sgz()
    {
      return MatrixContainer.GetRealMatrix(dt_rho0_sgz);
    };

    /// Get the rhox matrix from the container
    TRealMatrix& Get_rhox()
    {
      return MatrixContainer.GetRealMatrix(rhox);
    };
    /// Get the rhoy matrix from the container
    TRealMatrix& Get_rhoy()
    {
      return MatrixContainer.GetRealMatrix(rhoy);
    };
    /// Get the rhoz matrix from the container
    TRealMatrix& Get_rhoz()
    {
      return MatrixContainer.GetRealMatrix(rhoz);
    };
    /// Get the rho0 matrix from the container
    TRealMatrix& Get_rho0()
    {
      return MatrixContainer.GetRealMatrix(rho0);
    };

    /// Get the ddx_k_shift_pos matrix from the container
    TComplexMatrix& Get_ddx_k_shift_pos()
    {
      return MatrixContainer.GetComplexMatrix(ddx_k_shift_pos);
    };
    /// Get the ddy_k_shift_pos matrix from the container
    TComplexMatrix& Get_ddy_k_shift_pos()
    {
      return MatrixContainer.GetComplexMatrix(ddy_k_shift_pos);
    };
    /// Get the ddz_k_shift_pos matrix from the container
    TComplexMatrix& Get_ddz_k_shift_pos()
    {
      return MatrixContainer.GetComplexMatrix(ddz_k_shift_pos);
    };
    /// Get the ddx_k_shift_neg matrix from the container
    TComplexMatrix& Get_ddx_k_shift_neg()
    {
      return MatrixContainer.GetComplexMatrix(ddx_k_shift_neg);
    };
    /// Get the ddy_k_shift_neg matrix from the container
    TComplexMatrix& Get_ddy_k_shift_neg()
    {
      return MatrixContainer.GetComplexMatrix(ddy_k_shift_neg);
    };
    /// Get the ddz_k_shift_neg matrix from the container
    TComplexMatrix& Get_ddz_k_shift_neg()
    {
      return MatrixContainer.GetComplexMatrix(ddz_k_shift_neg);
    };


    /// Get the pml_x_sgx matrix from the container
    TRealMatrix& Get_pml_x_sgx()
    {
      return MatrixContainer.GetRealMatrix(pml_x_sgx);
    };
    /// Get the pml_y_sgy matrix from the container
    TRealMatrix& Get_pml_y_sgy()
    {
      return MatrixContainer.GetRealMatrix(pml_y_sgy);
    };
    /// Get the pml_z_sgz matrix from the container
    TRealMatrix& Get_pml_z_sgz()
    {
      return MatrixContainer.GetRealMatrix(pml_z_sgz);
    };

    /// Get the pml_x matrix from the container
    TRealMatrix& Get_pml_x()
    {
      return MatrixContainer.GetRealMatrix(pml_x);
    };
    /// Get the pml_y matrix from the container
    TRealMatrix& Get_pml_y()
    {
      return MatrixContainer.GetRealMatrix(pml_y);
    };
    /// Get the pml_z matrix from the container
    TRealMatrix& Get_pml_z()
    {
      return MatrixContainer.GetRealMatrix(pml_z);
    };


    /// Get the dxudxn matrix from the container
    TRealMatrix& Get_dxudxn()
    {
      return MatrixContainer.GetRealMatrix(dxudxn);
    };
    /// Get the dyudyn matrix from the container
    TRealMatrix& Get_dyudyn()
    {
      return MatrixContainer.GetRealMatrix(dyudyn);
    };
    /// Get the dzudzn matrix from the container
    TRealMatrix& Get_dzudzn()
    {
      return MatrixContainer.GetRealMatrix(dzudzn);
    };

    /// Get the dxudxn_sgx matrix from the container
    TRealMatrix& Get_dxudxn_sgx()
    {
      return MatrixContainer.GetRealMatrix(dxudxn_sgx);
    };
    /// Get the dyudyn_sgy matrix from the container
    TRealMatrix& Get_dyudyn_sgy()
    {
      return MatrixContainer.GetRealMatrix(dyudyn_sgy);
    };
    /// Get the dzudzn_sgz matrix from the container
    TRealMatrix& Get_dzudzn_sgz()
    {
      return MatrixContainer.GetRealMatrix(dzudzn_sgz);
    };


    /// Get the BonA matrix from the container
    TRealMatrix& Get_BonA()
    {
      return MatrixContainer.GetRealMatrix(BonA);
    };
    /// Get the absorb_tau matrix from the container
    TRealMatrix& Get_absorb_tau()
    {
      return MatrixContainer.GetRealMatrix(absorb_tau);
    };
    /// Get the absorb_eta matrix from the container
    TRealMatrix& Get_absorb_eta()
    {
      return MatrixContainer.GetRealMatrix(absorb_eta);
    };

    /// Get the absorb_nabla1 matrix from the container
    TRealMatrix& Get_absorb_nabla1()
    {
      return MatrixContainer.GetRealMatrix(absorb_nabla1);
    };
    /// Get the absorb_nabla2 matrix from the container
    TRealMatrix& Get_absorb_nabla2()
    {
      return MatrixContainer.GetRealMatrix(absorb_nabla2);
    };


    //-- Index matrices --//

    /// Get the sensor_mask_index matrix from the container
    TLongMatrix& Get_sensor_mask_index()
    {
      return MatrixContainer.GetLongMatrix(sensor_mask_index);
    };
    /// Get the sensor_mask_corners matrix from the container
    TLongMatrix& Get_sensor_mask_corners()
    {
      return MatrixContainer.GetLongMatrix(sensor_mask_corners);
    };

    /// Get the u_source_index matrix from the container
    TLongMatrix& Get_u_source_index()
    {
      return MatrixContainer.GetLongMatrix(u_source_index);
    };
    /// Get the p_source_index matrix from the container
    TLongMatrix& Get_p_source_index()
    {
      return MatrixContainer.GetLongMatrix(p_source_index);
    };
    /// Get the delay_mask matrix from the container
    TLongMatrix& Get_delay_mask()
    {
      return MatrixContainer.GetLongMatrix(delay_mask);
    }


    //-- sources  --//

    /// Get the transducer_source_input matrix from the container
    TRealMatrix& Get_transducer_source_input()
    {
      return MatrixContainer.GetRealMatrix(transducer_source_input);
    };

    /// Get the p_source_input matrix from the container
    TRealMatrix& Get_p_source_input()
    {
      return MatrixContainer.GetRealMatrix(p_source_input);
    };

    /// Get the p0_source_input from the container
    TRealMatrix& Get_p0_source_input()
    {
      return MatrixContainer.GetRealMatrix(p0_source_input);
    };


    /// Get the ux_source_input matrix from the container
    TRealMatrix& Get_ux_source_input()
    {
      return MatrixContainer.GetRealMatrix(ux_source_input);
    };
    /// Get the uy_source_input matrix from the container
    TRealMatrix& Get_uy_source_input()
    {
      return MatrixContainer.GetRealMatrix(uy_source_input);
    };
    /// Get the uz_source_input matrix from the container
    TRealMatrix& Get_uz_source_input()
    {
      return MatrixContainer.GetRealMatrix(uz_source_input);
    };


    //--Temporary matrices --//

    /// Get the Temp_1_RS3D matrix from the container
    TRealMatrix& Get_Temp_1_RS3D()
    {
      return MatrixContainer.GetRealMatrix(Temp_1_RS3D);
    };
    /// Get the Temp_2_RS3D matrix from the container
    TRealMatrix& Get_Temp_2_RS3D()
    {
      return MatrixContainer.GetRealMatrix(Temp_2_RS3D);
    };
    /// Get the Temp_3_RS3D matrix from the container
    TRealMatrix& Get_Temp_3_RS3D()
    {
      return MatrixContainer.GetRealMatrix(Temp_3_RS3D);
    };


    /// Get the FFT_X_temp from the container
    TFFTWComplexMatrix& Get_FFT_X_temp()
    {
      return MatrixContainer.GetFFTWComplexMatrix(FFT_X_temp);
    };
    /// Get the FFT_Y_temp from the container
    TFFTWComplexMatrix& Get_FFT_Y_temp()
    {
      return MatrixContainer.GetFFTWComplexMatrix(FFT_Y_temp);
    };
    /// Get the FFT_Z_temp from the container
    TFFTWComplexMatrix& Get_FFT_Z_temp()
    {return MatrixContainer.GetFFTWComplexMatrix(FFT_Z_temp);
    };

  private:

    /// Matrix container with all the matrix classes
    TMatrixContainer MatrixContainer;
    /// Output stream container
    TOutputStreamContainer OutputStreamContainer;

    /// Percentage of the simulation done
    int                ActPercent;

    /// Global parameters of the simulation
    TParameters *      Parameters;

    /// Total time of the simulation
    TTimeMesssure       TotalTime;
    /// Pre-processing time of the simulation
    TTimeMesssure       PreProcessingTime;
    /// Data load time of the simulation
    TTimeMesssure       DataLoadTime;
    /// Simulation time of the simulation
    TTimeMesssure       SimulationTime;
    /// Post-processing time of the simulation
    TTimeMesssure       PostProcessingTime;
    /// Iteration time of the simulation
    TTimeMesssure       IterationTime;


};// end of  TKSpace3DSolver
//------------------------------------------------------------------------------

#endif	/* TKSPACE3DSOLVER_H */

