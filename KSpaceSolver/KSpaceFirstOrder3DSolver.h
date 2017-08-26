/**
 * @file        KSpaceFirstOrder3DSolver.h
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The header file containing the main class of the project
 *              responsible for the entire simulation.
 *
 * @version     kspaceFirstOrder3D 2.16
 *
 * @date        12 July      2012, 10:27 (created)\n
 *              26 August    2017, 17:14 (revised)
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

#ifndef TKSPACE3DSOLVER_H
#define	TKSPACE3DSOLVER_H

#include <fftw3.h>

#include <Parameters/Parameters.h>

#include <MatrixClasses/MatrixContainer.h>
#include <MatrixClasses/RealMatrix.h>
#include <MatrixClasses/ComplexMatrix.h>
#include <MatrixClasses/IndexMatrix.h>
#include <OutputStreams/BaseOutputStream.h>
#include <MatrixClasses/VelocityMatrix.h>
#include <MatrixClasses/FftwComplexMatrix.h>

#include <Utils/TimeMeasure.h>

using namespace std;


/**
 * @class TKSpaceFirstOrder3DSolver
 * @brief Class responsible for running the k-space first order 3D method.
 * @details Class responsible for running the k-space first order 3D method. This
 *          class maintain the whole k-wave (implements the time loop).
 *
 */
class TKSpaceFirstOrder3DSolver
{
  public:
    /// Constructor.
    TKSpaceFirstOrder3DSolver();
    /// Destructor.
    virtual ~TKSpaceFirstOrder3DSolver();


    /// Memory allocation.
    virtual void AllocateMemory();
    /// Memory deallocation.
    virtual void FreeMemory();

    /// Load simulation data from the input file.
    virtual void LoadInputData();

    /// Compute the 3D kspace first order simulation.
    virtual void Compute();

    /// Print parameters of the simulation.
    virtual void PrintParametersOfSimulation(FILE * file);

    /// Get memory usage in MB.
    virtual size_t ShowMemoryUsageInMB();

    /// Get code name.
    string GetCodeName() {return "kspaceFirstOrder3D-OMP v1.1"; };

    /// Print the code name and license.
    void   PrintFullNameCodeAndLicense(FILE * file);

    /// Set processor affinity.
    void   SetProcessorAffinity();

    /// Get total simulation time.
    double GetTotalTime()          const { return TotalTime.getElapsedTime();         };

    /// Get pre-processing time.
    double GetPreProcessingTime()  const { return PreProcessingTime.getElapsedTime(); };

    /// Get data load time.
    double GetDataLoadTime()       const { return DataLoadTime.getElapsedTime();      };

    /// Get simulation time (time loop).
    double GetSimulationTime()     const { return SimulationTime.getElapsedTime();    };

    /// Get post-processing time.
    double GetPostProcessingTime() const { return PostProcessingTime.getElapsedTime();};

    /// Get total simulation time cumulated over all legs.
    double GetCumulatedTotalTime()          const { return TotalTime.getElapsedTimeOverAllLegs();          };
    /// Get pre-processing time cumulated over all legs.
    double GetCumulatedPreProcessingTime()  const { return PreProcessingTime.getElapsedTimeOverAllLegs();  };
    /// Get data load time cumulated over all legs.
    double GetCumulatedDataLoadTime()       const { return DataLoadTime.getElapsedTimeOverAllLegs();       };
    /// Get simulation time (time loop) cumulated over all legs.
    double GetCumulatedSimulationTime()     const { return SimulationTime.getElapsedTimeOverAllLegs();     };
    /// Get post-processing time cumulated over all legs.
    double GetCumulatedPostProcessingTime() const { return PostProcessingTime.getElapsedTimeOverAllLegs(); };

  protected:

    /// Copy constructor not allowed for public.
    TKSpaceFirstOrder3DSolver(const TKSpaceFirstOrder3DSolver& src);
    /// operator = not allowed for public.
    TKSpaceFirstOrder3DSolver& operator = (const TKSpaceFirstOrder3DSolver& src);

    /// Initialize FFT plans.
    void InitializeFFTWPlans();

    /// Compute pre-processing phase.
    void PreProcessingPhase();

    /// Compute the main time loop of the kspaceFirstOrder3D.
    void Compute_MainLoop();

    /// Post processing, and closing the output streams.
    void PostPorcessing();

    /// Store sensor data.
    void StoreSensorData();

    /// Save checkpoint data.
    void SaveCheckpointData();

    /// Write statistics and header into the output file.
    void WriteOutputDataInfo();


    /// compute new values of for ux_sgx, uy_sgy, uz_sgz.
    void Compute_uxyz();

    /// Compute new values of for duxdx, duydy, dzdz.
    void Compute_duxyz();


    /// Compute new values of rhox, rhoy, rhoz for non-linear case.
    void Compute_rhoxyz_nonlinear();
    /// Compute new values of rhox, rhoy, rhoz for linear case.
    void Compute_rhoxyz_linear();

    /// Add u source to the particle velocity.
    void Add_u_source();

    /// Add in pressure source.
    void Add_p_source();


    /// Generate kappa matrix for non-absorbing media.
    void Generate_kappa();

    /// Generate kappa matrix, absorb_nabla1, absorb_nabla2 for absorbing media.
    void Generate_kappa_absorb_nabla1_absorb_nabla2();

    /// Generate absorb_tau, absorb_eta for heterogenous media.
    void Generate_absorb_tau_absorb_eta_matrix();

    /// Calculate dt ./ rho0 for non-uniform grids.
    void Calculate_dt_rho0_non_uniform();

    /// Calculate p0_source.
    void Calculate_p0_source();

    /// Calculate c^2.
    void Compute_c2();

    /// Compute part of the new velocity - gradient of p.
    void Compute_ddx_kappa_fft_p(RealMatrix&        X_Matrix,
                                 FftwComplexMatrix& FFT_X,
                                 FftwComplexMatrix& FFT_Y,
                                 FftwComplexMatrix& FFT_Z,
                                 const RealMatrix&  kappa,
                                 const ComplexMatrix& ddx,
                                 const ComplexMatrix& ddy,
                                 const ComplexMatrix& ddz);

    /// Calculate new p, non-linear case.
    void Compute_new_p_nonlinear();
    /// Calculate new p linear-case, absorbing.
    void Compute_new_p_linear();


    /// Calculate three temporary sums in the new pressure formula, non-linear absorbing case, SSE2 version.
    void Calculate_SumRho_BonA_SumDu_SSE2(RealMatrix& RHO_Temp,
                                          RealMatrix& BonA_Temp,
                                          RealMatrix& Sum_du);
    /// Calculate two temporary sums in the new pressure formula, linear absorbing case.
    void Calculate_SumRho_SumRhoDu(RealMatrix& Sum_rhoxyz,
                                   RealMatrix& Sum_rho0_du);

    /// Compute absorbing term with abosrb_nabla1 and absorb_nabla2, SSE2 version.
    void Compute_Absorb_nabla1_2_SSE2(FftwComplexMatrix& FFT_1,
                                      FftwComplexMatrix& FFT_2);

    /// Sum sub-terms to calculate new pressure, non-linear case.
    void Sum_Subterms_nonlinear(RealMatrix& Absorb_tau_temp,
                                RealMatrix& Absorb_eta_temp,
                                RealMatrix& BonA_temp);

    /// Sum sub-terms to calculate new pressure, linear case.
    void Sum_Subterms_linear(RealMatrix& Absorb_tau_temp,
                             RealMatrix& Absorb_eta_temp,
                             RealMatrix& Sum_rhoxyz);

    /// Sum sub-terms for new p, linear lossless case.
    void Sum_new_p_nonlinear_lossless();
    ///  Sum sub-terms for new p, linear lossless case.
    void Sum_new_p_linear_lossless();

    /// Calculate ux_shifted, uy_shifted and uz_shifted.
    void Calculate_shifted_velocity();

    /// Print progress statistics.
    void PrintStatisitcs();

    /// Print the header of the progress statistics.
    void PrintOtputHeader();

    /// Is time to checkpoint (save actual state on disk)?
    bool IsTimeToCheckpoint();

    /// Was the loop interrupted to checkpoint?
    bool IsCheckpointInterruption() const
    {
      return (mParameters.getTimeIndex() != mParameters.getNt());
    };

    /// Check the output file has the correct format and version.
    void CheckOutputFile();
    /// Reads the header of the output file and sets the cumulative elapsed time from the first log.
    void RestoreCumulatedElapsedFromOutputFile();

    /// Check the checkpoint file has the correct format and version.
    void CheckCheckpointFile();


     //------------------------- Get matrices --------------------------------//
    /// Get the kappa matrix from the container.
    RealMatrix& Get_kappa()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(kappa);
    };
    /// Get the c^2 matrix from the container.
    RealMatrix& Get_c2()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(c2);
    };

    /// Get the p matrix from the container.
    RealMatrix& Get_p()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(p);
    };

    /// Get the ux_sgx matrix from the container.
    VelocityMatrix& Get_ux_sgx()
    {
      return MatrixContainer.GetMatrix<VelocityMatrix>(ux_sgx);
    };
    /// Get the uy_sgy matrix from the container.
    VelocityMatrix& Get_uy_sgy()
    {
      return MatrixContainer.GetMatrix<VelocityMatrix>(uy_sgy);
    };
    /// Get the uz_sgz matrix from the container.
    VelocityMatrix& Get_uz_sgz()
    {
      return MatrixContainer.GetMatrix<VelocityMatrix>(uz_sgz);
    };

    /// Get the ux_shifted matrix from the container.
    VelocityMatrix& Get_ux_shifted()
    {
      return MatrixContainer.GetMatrix<VelocityMatrix>(ux_shifted);
    };
    /// Get the uy_shifted matrix from the container.
    VelocityMatrix& Get_uy_shifted()
    {
      return MatrixContainer.GetMatrix<VelocityMatrix>(uy_shifted);
    };
    /// Get the uz_shifted matrix from the container.
    VelocityMatrix& Get_uz_shifted()
    {
      return MatrixContainer.GetMatrix<VelocityMatrix>(uz_shifted);
    };

    /// Get the duxdx matrix from the container.
    RealMatrix& Get_duxdx()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(duxdx);
    };
    /// Get the duydy matrix from the container.
    RealMatrix& Get_duydy()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(duydy);
    };
    /// Get the duzdz matrix from the container.
    RealMatrix& Get_duzdz()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(duzdz);
    };

    /// Get the dt.*rho0_sgx matrix from the container.
    RealMatrix& Get_dt_rho0_sgx()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(dt_rho0_sgx);
    };
    /// Get the dt.*rho0_sgy matrix from the container.
    RealMatrix& Get_dt_rho0_sgy()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(dt_rho0_sgy);
    };
    /// Get the dt.*rho0_sgz matrix from the container.
    RealMatrix& Get_dt_rho0_sgz()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(dt_rho0_sgz);
    };

    /// Get the rhox matrix from the container.
    RealMatrix& Get_rhox()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(rhox);
    };
    /// Get the rhoy matrix from the container.
    RealMatrix& Get_rhoy()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(rhoy);
    };
    /// Get the rhoz matrix from the container.
    RealMatrix& Get_rhoz()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(rhoz);
    };
    /// Get the rho0 matrix from the container.
    RealMatrix& Get_rho0()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(rho0);
    };

    /// Get the ddx_k_shift_pos matrix from the container.
    ComplexMatrix& Get_ddx_k_shift_pos()
    {
      return MatrixContainer.GetMatrix<ComplexMatrix>(ddx_k_shift_pos);
    };
    /// Get the ddy_k_shift_pos matrix from the container.
    ComplexMatrix& Get_ddy_k_shift_pos()
    {
      return MatrixContainer.GetMatrix<ComplexMatrix>(ddy_k_shift_pos);
    };
    /// Get the ddz_k_shift_pos matrix from the container.
    ComplexMatrix& Get_ddz_k_shift_pos()
    {
      return MatrixContainer.GetMatrix<ComplexMatrix>(ddz_k_shift_pos);
    };
    /// Get the ddx_k_shift_neg matrix from the container.
    ComplexMatrix& Get_ddx_k_shift_neg()
    {
      return MatrixContainer.GetMatrix<ComplexMatrix>(ddx_k_shift_neg);
    };
    /// Get the ddy_k_shift_neg matrix from the container.
    ComplexMatrix& Get_ddy_k_shift_neg()
    {
      return MatrixContainer.GetMatrix<ComplexMatrix>(ddy_k_shift_neg);
    };
    /// Get the ddz_k_shift_neg matrix from the container.
    ComplexMatrix& Get_ddz_k_shift_neg()
    {
      return MatrixContainer.GetMatrix<ComplexMatrix>(ddz_k_shift_neg);
    };

    /// Get the x_shift_neg_r matrix from the container.
    ComplexMatrix& Get_x_shift_neg_r()
    {
      return MatrixContainer.GetMatrix<ComplexMatrix>(x_shift_neg_r);
    };
    /// Get the y_shift_neg_r from the container.
    ComplexMatrix& Get_y_shift_neg_r()
    {
      return MatrixContainer.GetMatrix<ComplexMatrix>(y_shift_neg_r);
    };
    /// Get the y_shift_neg_r from the container.
    ComplexMatrix& Get_z_shift_neg_r()
    {
      return MatrixContainer.GetMatrix<ComplexMatrix>(z_shift_neg_r);
    };

    /// Get the pml_x_sgx matrix from the container.
    RealMatrix& Get_pml_x_sgx()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(pml_x_sgx);
    };
    /// Get the pml_y_sgy matrix from the container.
    RealMatrix& Get_pml_y_sgy()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(pml_y_sgy);
    };
    /// Get the pml_z_sgz matrix from the container.
    RealMatrix& Get_pml_z_sgz()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(pml_z_sgz);
    };

    /// Get the pml_x matrix from the container.
    RealMatrix& Get_pml_x()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(pml_x);
    };
    /// Get the pml_y matrix from the container.
    RealMatrix& Get_pml_y()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(pml_y);
    };
    /// Get the pml_z matrix from the container.
    RealMatrix& Get_pml_z()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(pml_z);
    };


    /// Get the dxudxn matrix from the container.
    RealMatrix& Get_dxudxn()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(dxudxn);
    };
    /// Get the dyudyn matrix from the container.
    RealMatrix& Get_dyudyn()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(dyudyn);
    };
    /// Get the dzudzn matrix from the container.
    RealMatrix& Get_dzudzn()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(dzudzn);
    };

    /// Get the dxudxn_sgx matrix from the container.
    RealMatrix& Get_dxudxn_sgx()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(dxudxn_sgx);
    };
    /// Get the dyudyn_sgy matrix from the container.
    RealMatrix& Get_dyudyn_sgy()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(dyudyn_sgy);
    };
    /// Get the dzudzn_sgz matrix from the container.
    RealMatrix& Get_dzudzn_sgz()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(dzudzn_sgz);
    };


    /// Get the BonA matrix from the container.
    RealMatrix& Get_BonA()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(BonA);
    };
    /// Get the absorb_tau matrix from the container.
    RealMatrix& Get_absorb_tau()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(absorb_tau);
    };
    /// Get the absorb_eta matrix from the container.
    RealMatrix& Get_absorb_eta()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(absorb_eta);
    };

    /// Get the absorb_nabla1 matrix from the container.
    RealMatrix& Get_absorb_nabla1()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(absorb_nabla1);
    };
    /// Get the absorb_nabla2 matrix from the container.
    RealMatrix& Get_absorb_nabla2()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(absorb_nabla2);
    };


    //-- Index matrices --//

    /// Get the sensor_mask_index matrix from the container.
    IndexMatrix& Get_sensor_mask_index()
    {
      return MatrixContainer.GetMatrix<IndexMatrix>(sensor_mask_index);
    };
    /// Get the sensor_mask_corners matrix from the container.
    IndexMatrix& Get_sensor_mask_corners()
    {
      return MatrixContainer.GetMatrix<IndexMatrix>(sensor_mask_corners);
    };

    /// Get the u_source_index matrix from the container.
    IndexMatrix& Get_u_source_index()
    {
      return MatrixContainer.GetMatrix<IndexMatrix>(u_source_index);
    };
    /// Get the p_source_index matrix from the container.
    IndexMatrix& Get_p_source_index()
    {
      return MatrixContainer.GetMatrix<IndexMatrix>(p_source_index);
    };
    /// Get the delay_mask matrix from the container.
    IndexMatrix& Get_delay_mask()
    {
      return MatrixContainer.GetMatrix<IndexMatrix>(delay_mask);
    }


    //-- sources  --//

    /// Get the transducer_source_input matrix from the container.
    RealMatrix& Get_transducer_source_input()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(transducer_source_input);
    };

    /// Get the p_source_input matrix from the container.
    RealMatrix& Get_p_source_input()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(p_source_input);
    };

    /// Get the p0_source_input from the container.
    RealMatrix& Get_p0_source_input()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(p0_source_input);
    };


    /// Get the ux_source_input matrix from the container.
    RealMatrix& Get_ux_source_input()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(ux_source_input);
    };
    /// Get the uy_source_input matrix from the container.
    RealMatrix& Get_uy_source_input()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(uy_source_input);
    };
    /// Get the uz_source_input matrix from the container.
    RealMatrix& Get_uz_source_input()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(uz_source_input);
    };


    //--Temporary matrices --//

    /// Get the Temp_1_RS3D matrix from the container.
    RealMatrix& Get_Temp_1_RS3D()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(Temp_1_RS3D);
    };
    /// Get the Temp_2_RS3D matrix from the container.
    RealMatrix& Get_Temp_2_RS3D()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(Temp_2_RS3D);
    };
    /// Get the Temp_3_RS3D matrix from the container.
    RealMatrix& Get_Temp_3_RS3D()
    {
      return MatrixContainer.GetMatrix<RealMatrix>(Temp_3_RS3D);
    };


    /// Get the FFT_X_temp from the container.
    FftwComplexMatrix& Get_FFT_X_temp()
    {
      return MatrixContainer.GetMatrix<FftwComplexMatrix>(FFT_X_temp);
    };
    /// Get the FFT_Y_temp from the container.
    FftwComplexMatrix& Get_FFT_Y_temp()
    {
      return MatrixContainer.GetMatrix<FftwComplexMatrix>(FFT_Y_temp);
    };
    /// Get the FFT_Z_temp from the container.
    FftwComplexMatrix& Get_FFT_Z_temp()
    {
      return MatrixContainer.GetMatrix<FftwComplexMatrix>(FFT_Z_temp);
    };
    /// Get the FFT_shift_temp the container.
    FftwComplexMatrix& Get_FFT_shift_temp()
    {
      return MatrixContainer.GetMatrix<FftwComplexMatrix>(FFT_shift_temp);
    };

  private:

    /// Matrix container with all the matrix classes.
    TMatrixContainer MatrixContainer;
    /// Output stream container.
    TOutputStreamContainer OutputStreamContainer;

    /// Percentage of the simulation done.
    size_t             ActPercent;

    /// Global parameters of the simulation.
    Parameters&       mParameters;

    /// Total time of the simulation.
    TimeMeasure       TotalTime;
    /// Pre-processing time of the simulation.
    TimeMeasure       PreProcessingTime;
    /// Data load time of the simulation.
    TimeMeasure       DataLoadTime;
    /// Simulation time of the simulation.
    TimeMeasure       SimulationTime;
    /// Post-processing time of the simulation.
    TimeMeasure       PostProcessingTime;
    /// Iteration time of the simulation.
    TimeMeasure       IterationTime;

};// end of  TKSpaceFirstOrder3DSolver
//------------------------------------------------------------------------------

#endif	/* TKSPACE3DSOLVER_H */

