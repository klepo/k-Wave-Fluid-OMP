/**
 * @file      KSpaceFirstOrderSolver.h
 *
 * @author    Jiri Jaros \n
 *            Faculty of Information Technology\n
 *            Brno University of Technology \n
 *            jarosjir@fit.vutbr.cz
 *
 * @brief     The header file containing the main class of the project responsible for the entire simulation.
 *
 * @version   kspaceFirstOrder 2.17
 *
 * @date      12 July      2012, 10:27 (created)\n
 *            14 March     2019, 12:35 (revised)
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

#ifndef KSPACE_FIRST_ORDER_SOLVER_H
#define KSPACE_FIRST_ORDER_SOLVER_H

#ifdef __linux__
#include <unistd.h>
#endif

#ifdef _WIN64
#include <Windows.h>
#endif

#include <Parameters/Parameters.h>

#include <Containers/MatrixContainer.h>
#include <Containers/OutputStreamContainer.h>

#include <MatrixClasses/RealMatrix.h>
#include <MatrixClasses/ComplexMatrix.h>
#include <MatrixClasses/IndexMatrix.h>
#include <OutputStreams/BaseOutputStream.h>
#include <MatrixClasses/FftwComplexMatrix.h>

#include <Utils/TimeMeasure.h>

/**
 * @class   KSpaceFirstOrderSolver
 * @brief   Class responsible for running the k-space first order method in 2D and 3D media.
 * @details Class responsible for running the k-space first order method in 2D and 3D media. This class maintain
 *          the whole k-wave (implements the time loop).
 *
 */
class KSpaceFirstOrderSolver
{
  public:
    /// Constructor.
    KSpaceFirstOrderSolver();
    /// Copy constructor not allowed for public.
    KSpaceFirstOrderSolver(const KSpaceFirstOrderSolver&) = delete;
    /// Destructor.
    virtual ~KSpaceFirstOrderSolver();
    /// operator= not allowed for public.
    KSpaceFirstOrderSolver& operator=(const KSpaceFirstOrderSolver&) = delete;

    /// Memory allocation.
    virtual void allocateMemory();
    /// Memory deallocation.
    virtual void freeMemory();

    /**
     * @brief Load simulation data.
     *
     * If checkpointing is enabled, this may include reading data from checkpoint and output file.
     */
    virtual void loadInputData();
    /**
     * @brief This method computes k-space First Order 2D/3D simulation.
     *
     * It launches calculation on a given dataset going through FFT initialization, pre-processing, main loop and
     * post-processing phases.
     */
    virtual void compute();

    /**
     * @brief  Get memory usage in MB on the host side.
     * @return Memory consumed on the host side in MB.
     */
    virtual size_t getMemoryUsage() const;
    /**
     * @brief  Get available memory in MB.
     * @return Available memory in MB.
     */
    virtual size_t getAvailableMemory() const;
    /**
     * @brief  Get peak memory usage in MB.
     * @return Peak memory consumed in MB.
     */
    virtual size_t getPeakMemoryUsage() const;
    /**
     * @brief  Get current memory usage in MB.
     * @return Current memory usage on in MB.
     */
    virtual size_t getCurrentMemoryUsage() const;
    /**
     * @brief  Get total RAM memory in MB.
     * @return Total RAM memory in MB.
     */
    virtual size_t getTotalMemory() const;

    /**
     * @brief  Get code name - release code version.
     * @return Release code version.
     */
    std::string getCodeName() const;

    /// Print the code name and license.
    void   printFullCodeNameAndLicense() const;

    /**
     * @brief   Set processor affinity.
     * @warning This may not work on some OS, it should be done by user before launching the code. Moreover, the user
     * may want to change the thread placement, e.g. on NUMA systems. This the routine was disabled in ver 2.16
     */
    void   setProcessorAffinity();

    /**
     * @brief  Get total simulation time.
     * @return Total simulation time in seconds.
     */
    double getTotalTime()          const { return mTotalTime.getElapsedTime();          };
    /**
     * @brief  Get pre-processing time.
     * @return Pre-processing time in seconds.
     */
    double getPreProcessingTime()  const { return mPreProcessingTime.getElapsedTime();  };
    /**
     * @brief  Get data load time.
     * @return Time to load data in seconds.
     */
    double getDataLoadTime()       const { return mDataLoadTime.getElapsedTime();       };
    /**
     * @brief  Get simulation time (time loop).
     * @return Time to execute the simulation in seconds.
     */
    double getSimulationTime()     const { return mSimulationTime.getElapsedTime();     };
    /**
     * @brief  Get post-processing time.
     * @return Time to postprocess simulation data in seconds.
     */
    double getPostProcessingTime() const { return mPostProcessingTime.getElapsedTime(); };

    /**
     * @brief  Get total simulation time accumulated over all legs.
     * @return Total execution time in seconds accumulated over all legs.
     */
    double getCumulatedTotalTime()          const { return mTotalTime.getElapsedTimeOverAllLegs();          };
    /**
     * @brief  Get pre-processing time accumulated over all legs.
     * @return Time to load data in seconds accumulated over all legs.
     */
    double getCumulatedPreProcessingTime()  const { return mPreProcessingTime.getElapsedTimeOverAllLegs();  };
    /**
     * @brief  Get simulation time (time loop) accumulated over all legs.
     * @return Time to execute the simulation in seconds accumulated over all legs.
     */
    double getCumulatedDataLoadTime()       const { return mDataLoadTime.getElapsedTimeOverAllLegs();       };
    /**
     * @brief  Get simulation time (time loop) accumulated over all legs.
     * @return Time to execute the simulation in seconds accumulated over all legs.
     */
    double getCumulatedSimulationTime()     const { return mSimulationTime.getElapsedTimeOverAllLegs();     };
    /**
     * @brief  Get post-processing time accumulated over all legs.
     * @return Time to post-processing simulation data in seconds accumulated over all legs.
     */
    double getCumulatedPostProcessingTime() const { return mPostProcessingTime.getElapsedTimeOverAllLegs(); };

  protected:

    /// Initialize FFTW plans.
    void InitializeFftwPlans();

    /**
     * @brief  Compute pre-processing phase.
     * @tparam simulationDimension - Dimensionality of the simulation.
     *
     * Initialize all indices, pre-compute constants such as c^2, rho0Sgx * dt  and create kappa,
     * absorbEta, absorbTau, absorbNabla1, absorbNabla2  matrices.
     */
    template<Parameters::SimulationDimension simulationDimension>
    void preProcessing();
    /**
     * @brief Compute the main time loop of the kspaceFirstOrder solver.
     * @tparam simulationDimension - Dimensionality of the simulation.
     */
    template<Parameters::SimulationDimension simulationDimension>
    void computeMainLoop();
    /// Post processing, and closing the output streams.
    template<Parameters::SimulationDimension simulationDimension>
    void postProcessing();

    /// Store sensor data.
    void storeSensorData();
    /// Write statistics and header into the output file.
    void writeOutputDataInfo();
    /// Save checkpoint data and flush aggregated outputs into the output file.
    void saveCheckpointData();

    /**
     * @brief Compute average intensities from stored p and u without compression.
     * @tparam simulationDimension - Dimensionality of the simulation.
     */
    template<Parameters::SimulationDimension simulationDimension>
    void computeAverageIntensities();

    /**
     * @brief Compute average intensities from stored p and u compression coefficients.
     * @tparam simulationDimension - Dimensionality of the simulation.
     */
    template<Parameters::SimulationDimension simulationDimension>
    void computeAverageIntensitiesC();

    /**
     * @brief Compute Q term (volume rate of heat deposition) from average intensities.
     * @tparam simulationDimension - Dimensionality of the simulation.
     * @param intensityXAvgStreamIndex - Average intensity x stream index.
     * @param intensityYAvgStreamIndex - Average intensity y stream index.
     * @param intensityZAvgStreamIndex - Average intensity z stream index.
     * @param qTermStreamIdx - Q term stream index.
     */
    template<Parameters::SimulationDimension simulationDimension>
    void computeQTerm(OutputStreamContainer::OutputStreamIdx intensityXAvgStreamIndex,
                      OutputStreamContainer::OutputStreamIdx intensityYAvgStreamIndex,
                      OutputStreamContainer::OutputStreamIdx intensityZAvgStreamIndex,
                      OutputStreamContainer::OutputStreamIdx qTermStreamIdx);

    /**
     * @brief  Compute new values of acoustic velocity in all used dimensions (UxSgx, UySgy, UzSgz).
     * @tparam simulationDimension - Dimensionality of the simulation.
     *
     * <b>Matlab code:</b> \code
     *  p_k = fftn(p);
     *  ux_sgx = bsxfun(@times, pml_x_sgx, ...
     *       bsxfun(@times, pml_x_sgx, ux_sgx) ...
     *       - dt .* rho0_sgx_inv .* real(ifftn( bsxfun(@times, ddx_k_shift_pos, kappa .* fftn(p)) )) ...
     *       );
     *  uy_sgy = bsxfun(@times, pml_y_sgy, ...
     *       bsxfun(@times, pml_y_sgy, uy_sgy) ...
     *       - dt .* rho0_sgy_inv .* real(ifftn( bsxfun(@times, ddy_k_shift_pos, kappa .* fftn(p)) )) ...
     *       );
     *  uz_sgz = bsxfun(@times, pml_z_sgz, ...
     *       bsxfun(@times, pml_z_sgz, uz_sgz) ...
     *       - dt .* rho0_sgz_inv .* real(ifftn( bsxfun(@times, ddz_k_shift_pos, kappa .* fftn(p)) )) ...
     *       );
     \endcode
     */
    template<Parameters::SimulationDimension simulationDimension>
    void computeVelocity();

    /**
     * @brief  Compute new values of acoustic velocity gradients.
     * @tparam simulationDimension - Dimensionality of the simulation.
     *
     * <b>Matlab code:</b> \code
     *  duxdx = real(ifftn( bsxfun(@times, ddx_k_shift_neg, kappa .* fftn(ux_sgx)) ));
     *  duydy = real(ifftn( bsxfun(@times, ddy_k_shift_neg, kappa .* fftn(uy_sgy)) ));
     *  duzdz = real(ifftn( bsxfun(@times, ddz_k_shift_neg, kappa .* fftn(uz_sgz)) ));
     * \endcode
     */
    template<Parameters::SimulationDimension simulationDimension>
    void computeVelocityGradient();

    /**
     * @brief  Compute new values of acoustic density for nonlinear case.
     * @tparam simulationDimension - Dimensionality of the simulation.
     *
     * <b>Matlab code:</b> \code
     *  rho0_plus_rho = 2 .* (rhox + rhoy + rhoz) + rho0;
     *  rhox = bsxfun(@times, pml_x, bsxfun(@times, pml_x, rhox) - dt .* rho0_plus_rho .* duxdx);
     *  rhoy = bsxfun(@times, pml_y, bsxfun(@times, pml_y, rhoy) - dt .* rho0_plus_rho .* duydy);
     *  rhoz = bsxfun(@times, pml_z, bsxfun(@times, pml_z, rhoz) - dt .* rho0_plus_rho .* duzdz);
     * \endcode
     */
    template<Parameters::SimulationDimension simulationDimension>
    void computeDensityNonliner();

    /**
     * @brief  Compute new values of acoustic density for linear case.
     * @tparam simulationDimension - Dimensionality of the simulation.
     *
     * <b>Matlab code:</b> \code
     *  rhox = bsxfun(@times, pml_x, bsxfun(@times, pml_x, rhox) - dt .* rho0 .* duxdx);
     *  rhoy = bsxfun(@times, pml_y, bsxfun(@times, pml_y, rhoy) - dt .* rho0 .* duydy);
     *  rhoz = bsxfun(@times, pml_z, bsxfun(@times, pml_z, rhoz) - dt .* rho0 .* duzdz);
     * \endcode
     *
     */
    template<Parameters::SimulationDimension simulationDimension>
    void computeDensityLinear();

    /**
     * @brief  Compute acoustic pressure for nonlinear case.
     * @tparam simulationDimension - Dimensionality of the simulation.
     *
     * <b>Matlab code:</b> \code
     *  case 'lossless'
     *
     *      % calculate p using a nonlinear adiabatic equation of state
     *      p = c.^2 .* (rhox + rhoy + rhoz + medium.BonA .* (rhox + rhoy + rhoz).^2 ./ (2 .* rho0));
     *
     *  case 'absorbing'
     *      % calculate p using a nonlinear absorbing equation of state
     *      p = c.^2 .* (...
     *          (rhox + rhoy + rhoz) ...
     *          + absorb_tau .* real(ifftn( absorb_nabla1 .* fftn(rho0 .* (duxdx + duydy + duzdz)) ))...
     *          - absorb_eta .* real(ifftn( absorb_nabla2 .* fftn(rhox + rhoy + rhoz) ))...
     *          + medium.BonA .*(rhox + rhoy + rhoz).^2 ./ (2 .* rho0) ...
     *          );
     * \endcode
     */
    template<Parameters::SimulationDimension simulationDimension>
    void computePressureNonlinear();

    /**
     * @brief  Compute acoustic pressure for linear case.
     * @tparam simulationDimension - Dimensionality of the simulation.
     *
     * <b>Matlab code:</b> \code
     *  case 'lossless'
     *
     *      % calculate p using a linear adiabatic equation of state
     *      p = c.^2 .* (rhox + rhoy + rhoz);
     *
     *  case 'absorbing'
     *      % calculate p using a linear absorbing equation of state
     *      p = c.^2 .* ( ...
     *          (rhox + rhoy + rhoz) ...
     *          + absorb_tau .* real(ifftn( absorb_nabla1 .* fftn(rho0 .* (duxdx + duydy + duzdz)) )) ...
     *          - absorb_eta .* real(ifftn( absorb_nabla2 .* fftn(rhox + rhoy + rhoz) )) ...
     *          );
     *\endcode
     */
    template<Parameters::SimulationDimension simulationDimension>
    void computePressureLinear();


    /// Add in all velocity sources.
    void addVelocitySource();
    /**
      * @brief Add in velocity source terms.
      * @param [in] velocityMatrix      - Velocity matrix to add the source in.
      * @param [in] velocitySourceInput - Source input to add.
      * @param [in] velocitySourceIndex - Source geometry index matrix.
      */
    void computeVelocitySourceTerm(RealMatrix&        velocityMatrix,
                                   const RealMatrix&  velocitySourceInput,
                                   const IndexMatrix& velocitySourceIndex);

    /**
     * @brief  Add in pressure source.
     * @tparam simulationDimension - Dimensionality of the simulation.
     */
    template<Parameters::SimulationDimension simulationDimension>
    void addPressureSource();
    /**
     * @brief Calculate initial pressure source.
     * @tparam simulationDimension - Dimensionality of the simulation.
     *
     * <b>Matlab code:</b> \code
     *  % add the initial pressure to rho as a mass source
     *  p = source.p0;
     *  rhox = source.p0 ./ (3 .* c.^2);
     *  rhoy = source.p0 ./ (3 .* c.^2);
     *  rhoz = source.p0 ./ (3 .* c.^2);
     *
     *  % compute u(t = t1 + dt/2) based on the assumption u(dt/2) = -u(-dt/2)
     *  % which forces u(t = t1) = 0
     *  ux_sgx = dt .* rho0_sgx_inv .* real(ifftn( bsxfun(@times, ddx_k_shift_pos, kappa .* fftn(p)) )) / 2;
     *  uy_sgy = dt .* rho0_sgy_inv .* real(ifftn( bsxfun(@times, ddy_k_shift_pos, kappa .* fftn(p)) )) / 2;
     *  uz_sgz = dt .* rho0_sgz_inv .* real(ifftn( bsxfun(@times, ddz_k_shift_pos, kappa .* fftn(p)) )) / 2;
     * \endcode
     */
    template<Parameters::SimulationDimension simulationDimension>
    void addInitialPressureSource();
    /// Add transducer data source to velocity x component.
    void addTransducerSource();


    /// Generate kappa matrix for lossless media.
    void generateKappa();
    /// Generate sourceKappa matrix for additive sources.
    void generateSourceKappa();
    /// Generate kappa matrix, absorbNabla1, absorbNabla2 for absorbing medium.
    void generateKappaAndNablas();
    /// Generate absorbTau, absorbEta for heterogenous medium.
    void generateTauAndEta();

    /**
     * @brief Calculate dt ./ rho0 for nonuniform grids.
     * @tparam simulationDimension - Dimensionality of the simulation.
     */
    template<Parameters::SimulationDimension simulationDimension>
    void generateInitialDenisty();
    /// Calculate square of velocity
    void computeC2();

    /**
     * @brief  Compute velocity for the initial pressure problem, heterogeneous medium, uniform grid.
     * @tparam simulationDimension - Dimensionality of the simulation.
     *
     * <b> Matlab code: </b> \code
     *  ux_sgx = dt ./ rho0_sgx .* ifft(ux_sgx).
     *  uy_sgy = dt ./ rho0_sgy .* ifft(uy_sgy).
     *  uz_sgz = dt ./ rho0_sgz .* ifft(uz_sgz).
     * \endcode
     */
    template<Parameters::SimulationDimension simulationDimension>
    void computeInitialVelocityHeterogeneous();
    /**
     * @brief  Compute velocity for the initial pressure problem, homogeneous medium, uniform grid.
     * @tparam simulationDimension - Dimensionality of the simulation.
     *
     * <b> Matlab code: </b> \code
     *  ux_sgx = dt ./ rho0_sgx .* ifft(ux_sgx).
     *  uy_sgy = dt ./ rho0_sgy .* ifft(uy_sgy).
     *  uz_sgz = dt ./ rho0_sgz .* ifft(uz_sgz).
     * \endcode
     *
     */
     template<Parameters::SimulationDimension simulationDimension>
     void computeInitialVelocityHomogeneousUniform();
    /**
     * @brief  Compute acoustic velocity for initial pressure problem, homogenous medium, nonuniform grid.
     * @tparam simulationDimension - Dimensionality of the simulation.
     *
     * <b> Matlab code: </b> \code
     *  ux_sgx = dt ./ rho0_sgx .* dxudxn_sgx .* ifft(ux_sgx)
     *  uy_sgy = dt ./ rho0_sgy .* dyudxn_sgy .* ifft(uy_sgy)
     *  uz_sgz = dt ./ rho0_sgz .* dzudzn_sgz .* ifft(uz_sgz)
     * \endcode
     */
    template<Parameters::SimulationDimension simulationDimension>
    void computeInitialVelocityHomogeneousNonuniform();

    /**
     * @brief  Compute acoustic velocity for heterogeneous medium and a uniform grid.
     * @tparam simulationDimension - Dimensionality of the simulation.
     *
     * <b> Matlab code: </b> \code
     *  ux_sgx = bsxfun(@times, pml_x_sgx, bsxfun(@times, pml_x_sgx, ux_sgx) - dt .* rho0_sgx_inv .* real(ifftX)
     *  uy_sgy = bsxfun(@times, pml_y_sgy, bsxfun(@times, pml_y_sgy, uy_sgy) - dt .* rho0_sgy_inv .* real(ifftY)
     *  uz_sgz = bsxfun(@times, pml_z_sgz, bsxfun(@times, pml_z_sgz, uz_sgz) - dt .* rho0_sgz_inv .* real(ifftZ)
     * \endcode
     */
    template<Parameters::SimulationDimension simulationDimension>
    void computeVelocityHeterogeneous();

    /**
     * @brief  Compute acoustic velocity for homogeneous medium and a uniform grid.
     * @tparam simulationDimension - Dimensionality of the simulation.
     *
     * <b> Matlab code: </b> \code
     *  ux_sgx = bsxfun(@times, pml_x_sgx, bsxfun(@times, pml_x_sgx, ux_sgx) - dt .* rho0_sgx_inv .* real(ifftX)
     *  uy_sgy = bsxfun(@times, pml_y_sgy, bsxfun(@times, pml_y_sgy, uy_sgy) - dt .* rho0_sgy_inv .* real(ifftY)
     *  uz_sgz = bsxfun(@times, pml_z_sgz, bsxfun(@times, pml_z_sgz, uz_sgz) - dt .* rho0_sgz_inv .* real(ifftZ)
     *\endcode
     */
    template<Parameters::SimulationDimension simulationDimension>
    void computeVelocityHomogeneousUniform();
    /**
     * @brief  Compute acoustic velocity for homogenous medium and nonuniform grid.
     * @tparam simulationDimension - Dimensionality of the simulation.
     *
     * <b> Matlab code: </b> \code
     *  ux_sgx = bsxfun(@times, pml_x_sgx, bsxfun(@times, pml_x_sgx, ux_sgx)  ...
     *                  - dt .* rho0_sgx_inv .* dxudxnSgx.* real(ifftX))
     *  uy_sgy = bsxfun(@times, pml_y_sgy, bsxfun(@times, pml_y_sgy, uy_sgy) ...
     *                  - dt .* rho0_sgy_inv .* dyudynSgy.* real(ifftY)
     *  uz_sgz = bsxfun(@times, pml_z_sgz, bsxfun(@times, pml_z_sgz, uz_sgz)
     *                  - dt .* rho0_sgz_inv .* dzudznSgz.* real(ifftZ)
     *\endcode
     */
    template<Parameters::SimulationDimension simulationDimension>
    void computeVelocityHomogeneousNonuniform();

    /**
     * @brief  Compute part of the new velocity term - gradient of pressure.
     * @tparam simulationDimension - Dimensionality of the simulation.*
     *
     * <b>Matlab code:</b> \code
     *  bsxfun(\@times, ddx_k_shift_pos, kappa .* fftn(p))
     * \endcode
     */
    template<Parameters::SimulationDimension simulationDimension>
    void computePressureGradient();
    /**
     * @brief Calculate three temporary sums in the new pressure formula before taking the FFT,
     *        nonlinear absorbing case.
     *
     * @tparam simulationDimension      - Dimensionality of the simulation.
     * @tparam bOnAScalarFlag           - is nonlinearity homogenous?
     * @tparam rho0ScalarFlag           - is density homogeneous?
     * @param [out] densitySum          - rhoX + rhoY + rhoZ
     * @param [out] nonlinearTerm       - BOnA + densitySum ^2 / 2 * rho0
     * @param [out] velocityGradientSum - rho0* (duxdx + duydy + duzdz)
     */
    template<Parameters::SimulationDimension simulationDimension,
             bool bOnAScalarFlag,
             bool rho0ScalarFlag>
    void computePressureTermsNonlinear(RealMatrix& densitySum,
                                       RealMatrix& nonlinearTerm,
                                       RealMatrix& velocityGradientSum);
    /**
     * @brief Calculate two temporary sums in the new pressure formula before taking the FFT,
     *        linear absorbing case.
     *
     * @tparam simulationDimension      - Dimensionality of the simulation.
     * @param [out] densitySum          - rhoxSgx + rhoySgy + rhozSgz
     * @param [out] velocityGradientSum - rho0* (duxdx + duydy + duzdz)
     */
    template<Parameters::SimulationDimension simulationDimension>
    void computePressureTermsLinear(RealMatrix& densitySum,
                                    RealMatrix& velocityGradientSum);

    /**
     * @brief Compute absorbing term with abosrbNabla1 and absorbNabla2.
     *
     * @param [in,out] fftPart1 - fftPart1 = absorbNabla1 .* fftPart1
     * @param [in,out] fftPart2 - fftPart1 = absorbNabla1 .* fftPart2
     */
    void computeAbsorbtionTerm(FftwComplexMatrix& fftPart1,
                               FftwComplexMatrix& fftPart2);

    /**
     * @brief Sum sub-terms to calculate new pressure, after FFTs, nonlinear case.
     *
     * @tparam c0ScalarFlag        - is sound speed homogeneous?
     * @tparam areTauAndEtaScalars - are absorbTau and absorbEeta scalars
     * @param [in] absorbTauTerm    - tau component
     * @param [in] absorbEtaTerm    - eta component  of the pressure term
     * @param [in] nonlinearTerm    - rho0 * (duxdx + duydy + duzdz)
     */
    template<bool c0ScalarFlag, bool areTauAndEtaScalars>
    void sumPressureTermsNonlinear(const RealMatrix& absorbTauTerm,
                                   const RealMatrix& absorbEtaTerm,
                                   const RealMatrix& nonlinearTerm);
    /**
     * @brief Sum sub-terms to calculate new pressure, after FFTs, linear case.
     *
     * @tparam c0ScalarFlag        - is sound speed homogeneous?
     * @tparam areTauAndEtaScalars - are absorbTau and absorbEeta homogeneous?
     * @param [in] absorbTauTerm - tau component
     * @param [in] absorbEtaTerm - eta component  of the pressure term
     * @param [in] densitySum    - Sum of three components of density (rhoXSgx + rhoYSgy + rhoZSgx)
     */
    template<bool c0ScalarFlag, bool areTauAndEtaScalars>
    void sumPressureTermsLinear(const RealMatrix& absorbTauTerm,
                                const RealMatrix& absorbEtaTerm,
                                const RealMatrix& densitySum);

    /**
     * @brief Sum sub-terms for new pressure, linear lossless case.
     *
     * @tparam simulationDimension - Dimensionality of the simulation.
     * @tparam c0ScalarFlag        - is sound speed homogeneous?
     * @tparam nonlinearFlag       - is nonlinearity homogeneous?
     * @tparam rho0ScalarFlag      - is density homogeneous?
     */
    template<Parameters::SimulationDimension simulationDimension,
             bool c0ScalarFlag,
             bool nonlinearFlag,
             bool rho0ScalarFlag>
    void sumPressureTermsNonlinearLossless();

    /**
     * @brief Sum sub-terms for new pressure, linear lossless case.
     * @tparam simulationDimension - Dimensionality of the simulation.
     */
    template<Parameters::SimulationDimension simulationDimension>
    void sumPressureTermsLinearLossless();

    /**
     * @brief Compute shifted velocity for --u_non_staggered flag.
     * @tparam simulationDimension - Dimensionality of the simulation.
     *
     * <b>Matlab code:</b> \code
     *  ux_shifted = real(ifft(bsxfun(\@times, x_shift_neg, fft(ux_sgx, [], 1)), [], 1));
     *  uy_shifted = real(ifft(bsxfun(\@times, y_shift_neg, fft(uy_sgy, [], 2)), [], 2));
     *  uz_shifted = real(ifft(bsxfun(\@times, z_shift_neg, fft(uz_sgz, [], 3)), [], 3));
     * \endcode
     */
    template<Parameters::SimulationDimension simulationDimension>
    void computeShiftedVelocity();

    /// Print progress statistics.
    void printStatistics();

    /**
     * @brief  Was the loop interrupted to checkpoint?
     * @return true if the simulation has been interrupted.
     */
    bool isCheckpointInterruption() const;

    /**
     * @brief Check the output file has the correct format and version.
     * @throw ios::failure - If an error happens.
     */
    void checkOutputFile();
    /**
     * @brief Check the file type and the version of the checkpoint file.
     * @throw ios::failure - If an error happens
     */
    void checkCheckpointFile();

    /// Reads the header of the output file and sets the cumulative elapsed time from the first log.
    void loadElapsedTimeFromOutputFile();

    /**
     * @brief Compute 1D index using 3 spatial coordinates and the size of the matrix.
     * @param [in] z              - z coordinate
     * @param [in] y              - y coordinate
     * @param [in] x              - x coordinate
     * @param [in] dimensionSizes - Size of the matrix.
     * @return
     */
    size_t get1DIndex(const size_t          z,
                      const size_t          y,
                      const size_t          x,
                      const DimensionSizes& dimensionSizes);
    //----------------------------------------------- Get matrices ---------------------------------------------------//
    /**
     * @brief  Get the kappa matrix from the container.
     * @return kappa matrix
     */
    RealMatrix& getKappa()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kKappa);
    };
    /**
     * @brief  Get the sourceKappa matrix from the container.
     * @return kappa matrix
     */
    RealMatrix& getSourceKappa()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kSourceKappa);
    };

    /**
     * @brief  Get the c^2 matrix from the container.
     * @return c^2 matrix.
     */
    RealMatrix& getC2()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kC2);
    };

    /**
     * @brief  Get pressure matrix
     * @return Pressure matrix
     */
    RealMatrix& getP()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kP);
    };

    //--------------------------------------------- Velocity matrices ------------------------------------------------//
    /**
     * @brief  Get velocity matrix on staggered grid in x direction.
     * @return Velocity matrix on staggered grid.
     */
    RealMatrix& getUxSgx()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kUxSgx);
    };
    /**
     * @brief  Get velocity matrix on staggered grid in y direction.
     * @return Velocity matrix on staggered grid.
     */
    RealMatrix& getUySgy()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kUySgy);
    };
    /**
     * @brief  Get velocity matrix on staggered grid in z direction.
     * @return Velocity matrix on staggered grid.
     */
    RealMatrix& getUzSgz()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kUzSgz);
    };

    /**
     * @brief  Get velocity shifted on normal grid in x direction.
     * @return Unstaggeted velocity matrix.
     */
    RealMatrix& getUxShifted()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kUxShifted);
    };
    /**
     * @brief  Get velocity shifted on normal grid in y direction.
     * @return Unstaggered velocity matrix.
     */
    RealMatrix& getUyShifted()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kUyShifted);
    };
    /**
     * @brief  Get velocity shifted on normal grid in z direction.
     * @return Unstaggered velocity matrix.
     */
    RealMatrix& getUzShifted()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kUzShifted);
    };

    //----------------------------------------- Velocity gradient matrices -------------------------------------------//
    /**
     * @brief  Get velocity gradient on in x direction.
     * @return Velocity gradient matrix.
     */
    RealMatrix& getDuxdx()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kDuxdx);
    };
    /**
     * @brief  Get velocity gradient on in y direction.
     * @return Velocity gradient matrix.
     */
    RealMatrix& getDuydy()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kDuydy);
    };
    /**
     * @brief  Get velocity gradient on in z direction.
     * @return Velocity gradient matrix.
     */
    RealMatrix& getDuzdz()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kDuzdz);
    };

    //---------------------------------------------- Density matrices ------------------------------------------------//
    /**
     * @brief  Get dt * rho0Sgx matrix (time step size * ambient velocity on staggered grid in x direction).
     * @return Density matrix
     */
    RealMatrix& getDtRho0Sgx()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kDtRho0Sgx);
    };
    /**
     * @brief  Get dt * rho0Sgy matrix (time step size * ambient velocity on staggered grid in y direction).
     * @return Density matrix
     */
    RealMatrix& getDtRho0Sgy()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kDtRho0Sgy);
    };
    /**
     * @brief  Get dt * rho0Sgz matrix (time step size * ambient velocity on staggered grid in z direction).
     * @return Density matrix
     */
    RealMatrix& getDtRho0Sgz()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kDtRho0Sgz);
    };

    /**
     * @brief  Get density matrix in x direction.
     * @return Density matrix.
     */
    RealMatrix& getRhoX()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kRhoX);
    };
    /**
     * @brief  Get density matrix in y direction.
     * @return Density matrix.
     */
    RealMatrix& getRhoY()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kRhoY);
    };
    /**
     * @brief  Get density matrix in z direction.
     * @return Density matrix.
     */
    RealMatrix& getRhoZ()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kRhoZ);
    };
    /**
     * @brief  Get ambient density matrix.
     * @return Density matrix.
     */
    RealMatrix& getRho0()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kRho0);
    };

    //----------------------------------------------- Shift matrices -------------------------------------------------//
    /**
     * @brief  Get positive Fourier shift in x.
     * @return Shift matrix.
     */
    ComplexMatrix& getDdxKShiftPos()
    {
      return mMatrixContainer.getMatrix<ComplexMatrix>(MatrixContainer::MatrixIdx::kDdxKShiftPosR);
    };
    /**
     * @brief  Get positive Fourier shift in y.
     * @return Shift matrix.
     */
    ComplexMatrix& getDdyKShiftPos()
    {
      return mMatrixContainer.getMatrix<ComplexMatrix>(MatrixContainer::MatrixIdx::kDdyKShiftPos);
    };
    /**
     * @brief  Get positive Fourier shift in z.
     * @return Shift matrix.
     */
    ComplexMatrix& getDdzKShiftPos()
    {
      return mMatrixContainer.getMatrix<ComplexMatrix>(MatrixContainer::MatrixIdx::kDdzKShiftPos);
    };
    /**
     * @brief  Get negative Fourier shift in x.
     * @return Shift matrix.
     */
    ComplexMatrix& getDdxKShiftNeg()
    {
      return mMatrixContainer.getMatrix<ComplexMatrix>(MatrixContainer::MatrixIdx::kDdxKShiftNegR);
    };
    /**
     * @brief  Get negative Fourier shift in y.
     * @return Shift matrix.
     */
    ComplexMatrix& getDdyKShiftNeg()
    {
      return mMatrixContainer.getMatrix<ComplexMatrix>(MatrixContainer::MatrixIdx::kDdyKShiftNeg);
    };
    /**
     * @brief  Get negative Fourier shift in z.
     * @return shift matrix.
     */
    ComplexMatrix& getDdzKShiftNeg()
    {
      return mMatrixContainer.getMatrix<ComplexMatrix>(MatrixContainer::MatrixIdx::kDdzKShiftNeg);
    };

    /**
     * @brief  Get negative shift for non-staggered velocity in x.
     * @return Shift matrix.
     */
    ComplexMatrix& getXShiftNegR()
    {
      return mMatrixContainer.getMatrix<ComplexMatrix>(MatrixContainer::MatrixIdx::kXShiftNegR);
    };
    /**
     * @brief  Get negative shift for non-staggered velocity in y.
     * @return Shift matrix.
     */
    ComplexMatrix& getYShiftNegR()
    {
      return mMatrixContainer.getMatrix<ComplexMatrix>(MatrixContainer::MatrixIdx::kYShiftNegR);
    };
    /**
     * @brief  Get negative shift for non-staggered velocity in z.
     * @return Shift matrix.
     */
    ComplexMatrix& getZShiftNegR()
    {
      return mMatrixContainer.getMatrix<ComplexMatrix>(MatrixContainer::MatrixIdx::kZShiftNegR);
    };

    //------------------------------------------------ PML matrices --------------------------------------------------//
    /**
     * @brief  Get PML on staggered grid x.
     * @return PML matrix.
     */
    RealMatrix& getPmlXSgx()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kPmlXSgx);
    };
    /**
     * @brief  Get PML on staggered grid y.
     * @return PML matrix.
     */
    RealMatrix& getPmlYSgy()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kPmlYSgy);
    };
    /**
     * @brief  Get PML on staggered grid z.
     * @return PML matrix.
     */
    RealMatrix& getPmlZSgz()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kPmlZSgz);
    };

    /**
     * @brief  Get PML in x.
     * @return PML matrix.
     */
    RealMatrix& getPmlX()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kPmlX);
    };
    /**
     * @brief  Get PML in y.
     * @return PML matrix.
     */
    RealMatrix& getPmlY()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kPmlY);
    };
    /**
     * @brief  Get PML in z.
     * @return PML matrix.
     */
    RealMatrix& getPmlZ()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kPmlZ);
    };

    //------------------------------------------- Nonlinear grid matrices --------------------------------------------//
    /**
     * @brief  Non uniform grid acoustic velocity in x.
     * @return Velocity matrix.
     */
    RealMatrix& getDxudxn()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kDxudxn);
    };
    /**
     * @brief  Non uniform grid acoustic velocity in y.
     * @return Velocity matrix.
     */
    RealMatrix& getDyudyn()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kDyudyn);
    };
    /**
     * @brief  Non uniform grid acoustic velocity in z.
     * @return Velocity matrix.
     */
    RealMatrix& getDzudzn()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kDzudzn);
    };
    /**
     * @brief  Non uniform grid acoustic velocity on staggered grid x.
     * @return Velocity matrix.
     */
    RealMatrix& getDxudxnSgx()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kDxudxnSgx);
    };
    /**
     * @brief  Non uniform grid acoustic velocity on staggered grid x.
     * @return Velocity matrix.
     */
    RealMatrix& getDyudynSgy()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kDyudynSgy);
    };
    /**
     * @brief  Non uniform grid acoustic velocity on staggered grid x.
     * @return Velocity matrix.
     */
    RealMatrix& getDzudznSgz()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kDzudznSgz);
    };

    //-------------------------------------- Nonlinear and absorption matrices ---------------------------------------//
    /**
     * @brief  Get B on A (nonlinear coefficient).
     * @return Nonlinear coefficient.
     */
    RealMatrix& getBOnA()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kBOnA);
    };
    /**
     * @brief  Get absorbing coefficient Tau.
     * @return Absorbing coefficient.
     */
    RealMatrix& getAbsorbTau()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kAbsorbTau);
    };
    /**
     * @brief  Get absorbing coefficient Eta.
     * @return Absorbing coefficient.
     */
    RealMatrix& getAbsorbEta()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kAbsorbEta);
    };

    /**
     * @brief  Get absorbing coefficient Nabla1.
     * @return Absorbing coefficient.
     */
    RealMatrix& getAbsorbNabla1()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kAbsorbNabla1);
    };
    /**
     * @brief  Get absorbing coefficient Nabla2.
     * @return Absorbing coefficient.
     */
    RealMatrix& getAbsorbNabla2()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kAbsorbNabla2);
    };

    //----------------------------------------------- Index matrices -------------------------------------------------//
    /**
     * @brief  Get linear sensor mask (spatial geometry of the sensor).
     * @return Sensor mask data.
     */
    IndexMatrix& getSensorMaskIndex()
    {
      return mMatrixContainer.getMatrix<IndexMatrix>(MatrixContainer::MatrixIdx::kSensorMaskIndex);
    };
    /**
     * @brief  Get cuboid corners sensor mask. (Spatial geometry of multiple sensors).
     * @return Sensor mask data.
     */
    IndexMatrix& getSensorMaskCorners()
    {
      return mMatrixContainer.getMatrix<IndexMatrix>(MatrixContainer::MatrixIdx::kSensorMaskCorners);
    };
    /**
     * @brief  Get velocity source geometry data.
     * @return Source geometry indices
     */
    IndexMatrix& getVelocitySourceIndex()
    {
      return mMatrixContainer.getMatrix<IndexMatrix>(MatrixContainer::MatrixIdx::kVelocitySourceIndex);
    };
    /**
     * @brief  Get pressure source geometry data.
     * @return Source geometry indices
     */
    IndexMatrix& getPressureSourceIndex()
    {
      return mMatrixContainer.getMatrix<IndexMatrix>(MatrixContainer::MatrixIdx::kPressureSourceIndex);
    };
    /**
     * @brief  Get delay mask for many types sources
     * @return delay mask.
     */
    IndexMatrix& getDelayMask()
    {
      return mMatrixContainer.getMatrix<IndexMatrix>(MatrixContainer::MatrixIdx::kDelayMask);
    }

    //-------------------------------------------------- Sources  ----------------------------------------------------//
    /**
     * @brief  Get transducer source input data (signal).
     * @return Transducer source input data.
     */
    RealMatrix& getTransducerSourceInput()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kTransducerSourceInput);
    };
    /**
     * @brief  Get pressure source input data (signal).
     * @return Pressure source input data.
     */
    RealMatrix& getPressureSourceInput()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kPressureSourceInput);
    };
    /**
     * @brief  Get initial pressure source input data (whole matrix).
     * @return Initial pressure source input data.
     */
    RealMatrix& getInitialPressureSourceInput()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kInitialPressureSourceInput);
    };

    /**
     * @brief  Get Velocity source input data in x direction.
     * @return Velocity source input data.
     */
    RealMatrix& GetVelocityXSourceInput()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kVelocityXSourceInput);
    };
    /**
     * @brief  Get Velocity source input data in y direction.
     * @return Velocity source input data.
     */
    RealMatrix& GetVelocityYSourceInput()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kVelocityYSourceInput);
    };
    /**
     * @brief  Get Velocity source input data in z direction.
     * @return Velocity source input data.
     */
    RealMatrix& getVelocityZSourceInput()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kVelocityZSourceInput);
    };


    //--------------------------------------------- Temporary matrices -----------------------------------------------//
    /**
     * @brief  Get first real 2D/3D temporary matrix.
     * @return Temporary real 2D/3D matrix.
     */
    RealMatrix& getTemp1RealND()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kTemp1RealND);
    };
    /**
     * @brief  Get second real 2D/3D temporary matrix.
     * @return Temporary real 2D/3D matrix.
     */
    RealMatrix& getTemp2RealND()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kTemp2RealND);
    };
    /**
     * @brief  Get third real 3D temporary matrix.T
     * This matrix is only present for 3D simulations,
     * @return Temporary real 3D matrix.
     */
    RealMatrix& getTemp3RealND()
    {
      return mMatrixContainer.getMatrix<RealMatrix>(MatrixContainer::MatrixIdx::kTemp3RealND);
    };

    /**
     * @brief  Get temporary matrix for 1D fft in x.
     * @return Temporary complex 3D matrix.
     */
    FftwComplexMatrix& getTempFftwX()
    {
      return mMatrixContainer.getMatrix<FftwComplexMatrix>(MatrixContainer::MatrixIdx::kTempFftwX);
    };
    /**
     * @brief  Get temporary matrix for 1D fft in y.
     * @return Temporary complex 3D matrix.
     */
    FftwComplexMatrix& getTempFftwY()
    {
      return mMatrixContainer.getMatrix<FftwComplexMatrix>(MatrixContainer::MatrixIdx::kTempFftwY);
    };
    /**
     * @brief  Get temporary matrix for 1D fft in z.
     * @return Temporary complex 3D matrix.
     */
    FftwComplexMatrix& getTempFftwZ()
    {
      return mMatrixContainer.getMatrix<FftwComplexMatrix>(MatrixContainer::MatrixIdx::kTempFftwZ);
    };
    /**
     * @brief  Get temporary matrix for fft shift.
     * @return Temporary complex 3D matrix.
     */
    FftwComplexMatrix& getTempFftwShift()
    {
      return mMatrixContainer.getMatrix<FftwComplexMatrix>(MatrixContainer::MatrixIdx::kTempFftwShift);
    };

  private:

    /// Matrix container with all the matrix classes.
    MatrixContainer       mMatrixContainer;
    /// Output stream container.
    OutputStreamContainer mOutputStreamContainer;
    /// Global parameters of the simulation.
    Parameters&           mParameters;

    /// Percentage of the simulation done.
    float                 mActPercent;

    /// Total time of the simulation.
    TimeMeasure mTotalTime;
    /// Pre-processing time of the simulation.
    TimeMeasure mPreProcessingTime;
    /// Data load time of the simulation.
    TimeMeasure mDataLoadTime;
    /// Simulation time of the simulation.
    TimeMeasure mSimulationTime;
    /// Post-processing time of the simulation.
    TimeMeasure mPostProcessingTime;
    /// Iteration time of the simulation.
    TimeMeasure mIterationTime;

    size_t mBlockSizeDefault = 0;
    size_t mBlockSizeDefaultC = 0;

};// end of KSpaceFirstOrderSolver
//----------------------------------------------------------------------------------------------------------------------

#endif	/* KSPACE_FIRST_ORDER_SOLVER_H */

