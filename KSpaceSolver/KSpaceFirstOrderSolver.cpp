/**
 * @file      KSpaceFirstOrderSolver.cpp
 *
 * @author    Jiri Jaros \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            jarosjir@fit.vutbr.cz
 *
 * @brief     The implementation file containing the main class of the project responsible for the entire simulation.
 *
 * @version   kspaceFirstOrder 2.17
 *
 * @date      12 July      2012, 10:27 (created) \n
 *            15 March     2019, 12:36 (revised)
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

// Linux build
#ifdef __linux__
  #include <sys/resource.h>
#endif

// Windows build
#ifdef _WIN64
  #define _CRT_SECURE_NO_WARNINGS
  #define _USE_MATH_DEFINES
  #include <Windows.h>
  #include <Psapi.h>
  #pragma comment(lib, "Psapi.lib")
#endif

#ifdef _OPENMP
  #include <omp.h>
#endif

#include <immintrin.h>
#include <cmath>
#include <ctime>
#include <limits>
#include <string.h>
#include <fstream>

#include <KSpaceSolver/KSpaceFirstOrderSolver.h>
#include <Containers/MatrixContainer.h>
#include <Containers/OutputStreamContainer.h>

#include <MatrixClasses/FftwComplexMatrix.h>
#include <Logger/Logger.h>

using std::ios;
/// shortcut for Simulation dimensions
using SD = Parameters::SimulationDimension;

//--------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------- Constants -----------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Public methods ---------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Constructor of the class.
 */
KSpaceFirstOrderSolver::KSpaceFirstOrderSolver():
        mMatrixContainer(), mOutputStreamContainer(),
        mParameters(Parameters::getInstance()),
        mActPercent(0l),
        mTotalTime(), mPreProcessingTime(), mDataLoadTime (), mSimulationTime(),
        mPostProcessingTime(), mIterationTime()
{
  mTotalTime.start();

  //Switch off HDF5 error messages
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);
}// end of KSpaceFirstOrderSolver
//----------------------------------------------------------------------------------------------------------------------

/**
 * Destructor of the class.
 */
KSpaceFirstOrderSolver::~KSpaceFirstOrderSolver()
{
  freeMemory();
}// end of KSpaceFirstOrderSolver
//----------------------------------------------------------------------------------------------------------------------

/**
 * The method allocates the matrix container, creates all matrices and creates all output streams
 * (however not allocating memory).
 */
void KSpaceFirstOrderSolver::allocateMemory()
{
  Logger::log(Logger::LogLevel::kBasic, kOutFmtMemoryAllocation);
  Logger::flush(Logger::LogLevel::kBasic);

  // add matrices into the container and create all matrices
  mMatrixContainer.init();
  mMatrixContainer.createMatrices();

  // add output streams into container
  //@todo Think about moving under LoadInputData routine...
  mOutputStreamContainer.init(mMatrixContainer);

  Logger::log(Logger::LogLevel::kBasic, kOutFmtDone);
}// end of allocateMemory
//----------------------------------------------------------------------------------------------------------------------

/**
 * The method frees all memory allocated by the class.
 */
void KSpaceFirstOrderSolver::freeMemory()
{
  mMatrixContainer.freeMatrices();

  mOutputStreamContainer.freeStreams();
}// end of freeMemory
//----------------------------------------------------------------------------------------------------------------------

/**
 * Load data from the input file provided by the Parameter class and creates the output time series streams.
 */
void KSpaceFirstOrderSolver::loadInputData()
{
  Logger::log(Logger::LogLevel::kBasic, kOutFmtDataLoading);
  Logger::flush(Logger::LogLevel::kBasic);
  // Start timer
  mDataLoadTime.start();

  // get handles
  Hdf5File& inputFile      = mParameters.getInputFile(); // file is opened (in Parameters)
  Hdf5File& outputFile     = mParameters.getOutputFile();
  Hdf5File& checkpointFile = mParameters.getCheckpointFile();

  Logger::log(Logger::LogLevel::kFull, kOutFmtNoDone);

  // Load data from disk
  Logger::log(Logger::LogLevel::kFull, kOutFmtReadingInputFile);
  Logger::flush(Logger::LogLevel::kFull);

  // Load data from disk
  mMatrixContainer.loadDataFromInputFile();

  // close the input file since we don't need it anymore.
  inputFile.close();

  Logger::log(Logger::LogLevel::kFull, kOutFmtDone);

  // The simulation does not use checkpointing or this is the first turn
  bool recoverFromCheckpoint = (mParameters.isCheckpointEnabled() &&
                                Hdf5File::canAccess(mParameters.getCheckpointFileName()));

  if (recoverFromCheckpoint)
  {
    //------------------------------------- Read data from the checkpoint file ---------------------------------------//
    Logger::log(Logger::LogLevel::kFull, kOutFmtReadingCheckpointFile);
    Logger::flush(Logger::LogLevel::kFull);

    // Open checkpoint file
    checkpointFile.open(mParameters.getCheckpointFileName());

    // Check the checkpoint file
    checkCheckpointFile();

    // read the actual value of t_index
    size_t checkpointedTimeIndex;
    checkpointFile.readScalarValue(checkpointFile.getRootGroup(), kTimeIndexName, checkpointedTimeIndex);
    mParameters.setTimeIndex(checkpointedTimeIndex);

    // Read necessary matrices from the checkpoint file
    mMatrixContainer.loadDataFromCheckpointFile();


    Logger::log(Logger::LogLevel::kFull, kOutFmtDone);

    //--------------------------------------- Read data from the output file -----------------------------------------//
    Logger::log(Logger::LogLevel::kFull, kOutFmtReadingOutputFile);
    Logger::flush(Logger::LogLevel::kFull);

    // Reopen output file for RW access
    outputFile.open(mParameters.getOutputFileName(), H5F_ACC_RDWR);
    //Read file header of the output file
    mParameters.getFileHeader().readHeaderFromOutputFile(outputFile);
    // Check the checkpoint file
    checkOutputFile();
    // Restore elapsed time
    loadElapsedTimeFromOutputFile();

    mOutputStreamContainer.reopenStreams();
    Logger::log(Logger::LogLevel::kFull, kOutFmtDone);
    checkpointFile.close();
  }
  else
  {
    if (mParameters.getOnlyPostProcessingFlag() &&
        (mParameters.getStoreIntensityAvgFlag()
         || mParameters.getStoreIntensityAvgCFlag()
         || mParameters.getStoreQTermFlag()
         || mParameters.getStoreQTermCFlag()) &&
        Hdf5File::canAccess(mParameters.getOutputFileName()))
    {
      // Open output file
      // TODO check existing datasets and their sizes (streams)
      outputFile.open(mParameters.getOutputFileName(), H5F_ACC_RDWR);
    }
    else
    {
      //------------------------------------ First round of multi-leg simulation ---------------------------------------//
      // Create the output file
      Logger::log(Logger::LogLevel::kFull, kOutFmtCreatingOutputFile);
      Logger::flush(Logger::LogLevel::kFull);

      outputFile.create(mParameters.getOutputFileName());
      Logger::log(Logger::LogLevel::kFull, kOutFmtDone);
    }
    // Create the steams, link them with the sampled matrices, however DO NOT allocate memory!
    mOutputStreamContainer.createStreams();
  }

 // Stop timer
  mDataLoadTime.stop();
  if (Logger::getLevel() != Logger::LogLevel::kFull)
  {
    Logger::log(Logger::LogLevel::kBasic, kOutFmtDone);
  }
}// end of loadInputData
//----------------------------------------------------------------------------------------------------------------------

/**
 * This method computes k-space First Order 2D/3D simulation.
 */
void KSpaceFirstOrderSolver::compute()
{
  // fft initialisation and preprocessing
  try
  {
    mPreProcessingTime.start();

    // initilaise all FFTW plans
    InitializeFftwPlans();

    // preprocessing phase generating necessary variables
    if (mParameters.isSimulation3D())
      preProcessing<Parameters::SimulationDimension::k3D>();
    else
      preProcessing<Parameters::SimulationDimension::k2D>();

    mPreProcessingTime.stop();
  }
  catch (const std::exception& e)
  {
    Logger::log(Logger::LogLevel::kBasic, kOutFmtFailed);
    Logger::log(Logger::LogLevel::kBasic, kOutFmtLastSeparator);

    Logger::errorAndTerminate(Logger::wordWrapString(e.what(),kErrFmtPathDelimiters, 9));
  }

  // Logger header for simulation
  Logger::log(Logger::LogLevel::kBasic, kOutFmtElapsedTime, mPreProcessingTime.getElapsedTime());
  Logger::log(Logger::LogLevel::kBasic, kOutFmtCompResourcesHeader);
  Logger::log(Logger::LogLevel::kBasic, kOutFmtCurrentMemory,   getMemoryUsage());

  // Main loop
  if (!mParameters.getOnlyPostProcessingFlag())
  {
    try
    {
      mSimulationTime.start();

      if (mParameters.isSimulation3D())
        computeMainLoop<Parameters::SimulationDimension::k3D>();
      else
        computeMainLoop<Parameters::SimulationDimension::k2D>();


      mSimulationTime.stop();

      Logger::log(Logger::LogLevel::kBasic,kOutFmtSimulationEndSeparator);
    }
    catch (const std::exception& e)
    {
      Logger::log(Logger::LogLevel::kBasic, kOutFmtSimulatoinFinalSeparator);
      Logger::errorAndTerminate(Logger::wordWrapString(e.what(),kErrFmtPathDelimiters, 9));
    }

    Logger::log(Logger::LogLevel::kBasic, kOutFmtElapsedTime, mSimulationTime.getElapsedTime());
  }

  // Post processing region
  mPostProcessingTime.start();
  try
  {
    if (!mParameters.getOnlyPostProcessingFlag() && isCheckpointInterruption())
    { // Checkpoint
      Logger::log(Logger::LogLevel::kBasic, kOutFmtCheckpointCompletedTimeSteps, mParameters.getTimeIndex());
      Logger::log(Logger::LogLevel::kBasic, kOutFmtCheckpointHeader);
      Logger::log(Logger::LogLevel::kBasic, kOutFmtCreatingCheckpoint);
      Logger::flush(Logger::LogLevel::kBasic);

      if (Logger::getLevel() == Logger::LogLevel::kFull)
      {
        Logger::log(Logger::LogLevel::kBasic, kOutFmtNoDone);
      }

      saveCheckpointData();

      if (Logger::getLevel() != Logger::LogLevel::kFull)
      {
        Logger::log(Logger::LogLevel::kBasic, kOutFmtDone);
      }
    }
    else
    { // Finish
      Logger::log(Logger::LogLevel::kBasic, kOutFmtSeparator);
      Logger::log(Logger::LogLevel::kBasic, kOutFmtPostProcessing);
      Logger::flush(Logger::LogLevel::kBasic);

      if (mParameters.isSimulation3D())
      {
        postProcessing<SD::k3D>();
      }
      else
      {
        postProcessing<SD::k2D>();
      }

      // if checkpointing is enabled and the checkpoint file was created in the past, delete it
      if (mParameters.isCheckpointEnabled())
      {
        std::remove(mParameters.getCheckpointFileName().c_str());
      }
      Logger::log(Logger::LogLevel::kBasic, kOutFmtDone);
    }
  }
  catch (const std::exception &e)
  {
    Logger::log(Logger::LogLevel::kBasic, kOutFmtFailed);
    Logger::log(Logger::LogLevel::kBasic, kOutFmtLastSeparator);

    Logger::errorAndTerminate(Logger::wordWrapString(e.what(), kErrFmtPathDelimiters,9));
  }
  mPostProcessingTime.stop();

  // Final data written
  try
  {
    writeOutputDataInfo();
    mParameters.getOutputFile().close();

    Logger::log(Logger::LogLevel::kBasic, kOutFmtElapsedTime, mPostProcessingTime.getElapsedTime());
  }
  catch (const std::exception &e)
  {
    Logger::log(Logger::LogLevel::kBasic, kOutFmtLastSeparator);
    Logger::errorAndTerminate(Logger::wordWrapString(e.what(), kErrFmtPathDelimiters, 9));
  }
}// end of compute()
//----------------------------------------------------------------------------------------------------------------------

/**
 * Get peak memory usage.
 */
size_t KSpaceFirstOrderSolver::getMemoryUsage() const
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
}// end of getMemoryUsage
//----------------------------------------------------------------------------------------------------------------------

/**
 * Get available memory.
 */
size_t KSpaceFirstOrderSolver::getAvailableMemory() const
{
  #ifdef __linux__
    std::string token;
    std::ifstream file("/proc/meminfo");
    while (file >> token)
    {
      if (token == "MemAvailable:")
      {
        unsigned long mem;
        if (file >> mem)
        {
          return mem / 1000;
        }
        else
        {
          return 0;
        }
      }
      // ignore rest of the line
      file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    return 0; // nothing found
  #endif
  // This is wrong, probably without cache and other.
  /*#ifdef __linux__
      long pages = sysconf(_SC_AVPHYS_PAGES);
      long page_size = sysconf(_SC_PAGE_SIZE);
      return pages * page_size / 1000000;
  #endif*/
  #ifdef _WIN64
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return size_t(status.ullAvailPhys) / 1000000;
  #endif
}// end of getAvailableMemory
//----------------------------------------------------------------------------------------------------------------------

/**
 * Get peak memory usage.
 */
size_t KSpaceFirstOrderSolver::getPeakMemoryUsage() const
{
  #ifdef __linux__
    // linux file contains this-process info
    FILE* file = fopen("/proc/self/status", "r");
    char buffer[1024] = "";
    int peakRealMem;
    //int peakVirtMem;
    // read the entire file
    while (fscanf(file, " %1023s", buffer) == 1)
    {
      if (strcmp(buffer, "VmHWM:") == 0)
      { // kilobytes
        fscanf(file, " %d", &peakRealMem);
      }
      /*if (strcmp(buffer, "VmPeak:") == 0)
      {
        fscanf(file, " %d", &peakVirtMem);
      }*/
    }
    fclose(file);
    return size_t(peakRealMem) / 1000;
    //return size_t(peakVirtMem) / 1000;
  #endif
  #ifdef _WIN64
    PROCESS_MEMORY_COUNTERS pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
    return size_t(pmc.PeakWorkingSetSize) / 1000000;
  #endif
}// end of getPeakMemoryUsage
//----------------------------------------------------------------------------------------------------------------------

/**
 * Get current memory usage.
 */
size_t KSpaceFirstOrderSolver::getCurrentMemoryUsage() const
{
  #ifdef __linux__
    // linux file contains this-process info
    FILE* file = fopen("/proc/self/status", "r");
    char buffer[1024] = "";
    int currRealMem;
    //int currVirtMem;
    // read the entire file
    while (fscanf(file, " %1023s", buffer) == 1)
    {
      if (strcmp(buffer, "VmRSS:") == 0)
      { // kilobytes
        fscanf(file, " %d", &currRealMem);
      }
      /*if (strcmp(buffer, "VmSize:") == 0)
      {
        fscanf(file, " %d", &currVirtMem);
      }*/
    }
    fclose(file);
    return size_t(currRealMem) / 1000;
    //return size_t(currVirtMem) / 1000;
  #endif
  #ifdef _WIN64
    PROCESS_MEMORY_COUNTERS pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
    return size_t(pmc.WorkingSetSize) / 1000000;
  #endif
}// end of getCurrentMemoryUsage
//----------------------------------------------------------------------------------------------------------------------

/**
 * Get total memory.
 */
size_t KSpaceFirstOrderSolver::getTotalMemory() const
{
  #ifdef __linux__
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size / 1000000;
  #endif
  #ifdef _WIN64
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return size_t(status.ullTotalPhys) / 1000000;
  #endif
}// end of getTotalMemory
//----------------------------------------------------------------------------------------------------------------------

/**
 * Get release code version.
 */
std::string KSpaceFirstOrderSolver::getCodeName() const
{
  return std::string(kOutFmtKWaveVersion);
}// end of getCodeName
//----------------------------------------------------------------------------------------------------------------------

/**
 * Print full code name and the license.
 */
void KSpaceFirstOrderSolver::printFullCodeNameAndLicense() const
{
  Logger::log(Logger::LogLevel::kBasic, kOutFmtBuildNoDataTime, 10, 11, __DATE__, 8, 8, __TIME__);

  if (mParameters.getGitHash() != "")
  {
    Logger::log(Logger::LogLevel::kBasic, kOutFmtVersionGitHash, mParameters.getGitHash().c_str());
  }
  Logger::log(Logger::LogLevel::kBasic, kOutFmtSeparator);


  // OS detection
  #ifdef __linux__
    Logger::log(Logger::LogLevel::kBasic, kOutFmtLinuxBuild);
  #elif __APPLE__
    Logger::log(Logger::LogLevel::kBasic, kOutFmtMacOsBuild);
  #elif _WIN32
    Logger::log(Logger::LogLevel::kBasic, kOutFmtWindowsBuild);
  #endif

  // Compiler detections
  #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    Logger::log(Logger::LogLevel::kBasic, kOutFmtGnuCompiler, __VERSION__);
  #endif
  #ifdef __INTEL_COMPILER
    Logger::log(Logger::LogLevel::kBasic, kOutFmtIntelCompiler, __INTEL_COMPILER);
  #endif
  #ifdef _MSC_VER
  Logger::log(Logger::LogLevel::kBasic, kOutFmtVisualStudioCompiler, _MSC_VER);
  #endif

     // instruction set
  #if (defined (__AVX2__))
    Logger::log(Logger::LogLevel::kBasic, kOutFmtAVX2);
  #elif (defined (__AVX__))
    Logger::log(Logger::LogLevel::kBasic, kOutFmtAVX);
  #elif (defined (__SSE4_2__))
    Logger::log(Logger::LogLevel::kBasic, kOutFmtSSE42);
  #elif (defined (__SSE4_1__))
    Logger::log(Logger::LogLevel::kBasic, kOutFmtSSE41);
  #elif (defined (__SSE3__))
    Logger::log(Logger::LogLevel::kBasic, kOutFmtSSE3);
  #elif (defined (__SSE2__))
    Logger::log(Logger::LogLevel::kBasic, kOutFmtSSE2);
  #endif

  Logger::log(Logger::LogLevel::kBasic, kOutFmtLicense);

}// end of printFullCodeNameAndLicense
//----------------------------------------------------------------------------------------------------------------------

/**
 * Set processor affinity.
 */
void KSpaceFirstOrderSolver::setProcessorAffinity()
{
  // Linux Build
  #ifdef __linux__
    //GNU compiler
    #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
      setenv("OMP_PROC_BIND","TRUE", 1);
    #endif

    #ifdef __INTEL_COMPILER
      setenv("KMP_AFFINITY","none", 1);
    #endif
  #endif

  // Windows build is always compiled by the Intel Compiler
  #ifdef _WIN64
    _putenv_s("KMP_AFFINITY","none");
  #endif
}//end of setProcessorAffinity
//----------------------------------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------ Protected methods -------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Initialize FFTW plans.
 */
void KSpaceFirstOrderSolver::InitializeFftwPlans()
{

  // initialization of FFTW library
  #ifdef _OPENMP
    fftwf_init_threads();
    fftwf_plan_with_nthreads(mParameters.getNumberOfThreads());
  #endif

  // The simulation does not use checkpointing or this is the first turn
  bool recoverFromPrevState = (mParameters.isCheckpointEnabled() &&
                               Hdf5File::canAccess(mParameters.getCheckpointFileName()));


  #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    // import FFTW wisdom if it is here
    if (recoverFromPrevState)
    {


      Logger::log(Logger::LogLevel::kFull, kOutFmtLoadingFftwWisdom);
      Logger::flush(Logger::LogLevel::kFull);
      // import FFTW wisdom
      try
      {
        // try to find the wisdom in the file that has the same name as the checkpoint file (different extension)
        FftwComplexMatrix::importWisdom();
        Logger::log(Logger::LogLevel::kFull, kOutFmtDone);
      }
      catch (const std::runtime_error& e)
      {
        Logger::log(Logger::LogLevel::kFull, kOutFmtFailed);
      }
    }
  #endif

  Logger::log(Logger::LogLevel::kBasic, kOutFmtFftPlans);
  Logger::flush(Logger::LogLevel::kBasic);

  if (!mParameters.getOnlyPostProcessingFlag())
  {
  // create real to complex plans
  getTempFftwX().createR2CFftPlanND(getP());
  getTempFftwY().createR2CFftPlanND(getP());
  if (mParameters.isSimulation3D())
  {
    getTempFftwZ().createR2CFftPlanND(getP());
  }

  // create real to complex plans
  getTempFftwX().createC2RFftPlanND(getP());
  getTempFftwY().createC2RFftPlanND(getP());
  if (mParameters.isSimulation3D())
  {
    getTempFftwZ().createC2RFftPlanND(getP());
  }
  }

  // if necessary, create 1D shift plans.
  // in this case, the matrix has a bit bigger dimensions to be able to store
  // shifted matrices.
  if (mParameters.getStoreVelocityNonStaggeredRawFlag()
      || mParameters.getStoreVelocityNonStaggeredCFlag()
      || mParameters.getStoreIntensityAvgFlag()
      || mParameters.getStoreQTermFlag()
      || mParameters.getStoreIntensityAvgCFlag()
      || mParameters.getStoreQTermCFlag())
  {
    // X shifts
    getTempFftwShift().createR2CFftPlan1DX(getUxShifted());
    getTempFftwShift().createC2RFftPlan1DX(getUxShifted());

    // Y shifts
    getTempFftwShift().createR2CFftPlan1DY(getUyShifted());
    getTempFftwShift().createC2RFftPlan1DY(getUyShifted());

    // Z shifts
    if (mParameters.isSimulation3D())
    {
      getTempFftwShift().createR2CFftPlan1DZ(getUzShifted());
      getTempFftwShift().createC2RFftPlan1DZ(getUzShifted());
    }
  }// end u_non_staggered

  Logger::log(Logger::LogLevel::kBasic, kOutFmtDone);
}// end of InitializeFftwPlans
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute pre-processing phase.
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::preProcessing()
{
  Logger::log(Logger::LogLevel::kBasic,kOutFmtPreProcessing);
  Logger::flush(Logger::LogLevel::kBasic);

  // get the correct sensor mask and recompute indices
  if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kIndex)
  {
    getSensorMaskIndex().recomputeIndicesToCPP();
  }

  if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kCorners)
  {
    getSensorMaskCorners().recomputeIndicesToCPP();
  }

  if (!mParameters.getOnlyPostProcessingFlag())
  {
    if ((mParameters.getTransducerSourceFlag() != 0) ||
        (mParameters.getVelocityXSourceFlag() != 0)  ||
        (mParameters.getVelocityYSourceFlag() != 0)  ||
        (mParameters.getVelocityZSourceFlag() != 0)
       )
    {
      getVelocitySourceIndex().recomputeIndicesToCPP();
    }

    if (mParameters.getTransducerSourceFlag() != 0)
    {
      getDelayMask().recomputeIndicesToCPP();
    }

    if (mParameters.getPressureSourceFlag() != 0)
    {
      getPressureSourceIndex().recomputeIndicesToCPP();
    }


    // compute dt / rho0_sg...
    if (!mParameters.getRho0ScalarFlag())
    { // non-uniform grid cannot be pre-calculated :-(
      // rho is matrix
      if (mParameters.getNonUniformGridFlag())
      {
        generateInitialDenisty<simulationDimension>();
      }
      else
      {
        getDtRho0Sgx().scalarDividedBy(mParameters.getDt());
        getDtRho0Sgy().scalarDividedBy(mParameters.getDt());
        if (simulationDimension == SD::k3D)
        {
          getDtRho0Sgz().scalarDividedBy(mParameters.getDt());
        }
      }
    }

    // generate different matrices
    if (mParameters.getAbsorbingFlag() != 0)
    {
      generateKappaAndNablas();
      generateTauAndEta();
    }
    else
    {
      generateKappa();
    }

    // Generate sourceKappa
    if (((mParameters.getVelocitySourceMode() == Parameters::SourceMode::kAdditive) ||
         (mParameters.getPressureSourceMode() == Parameters::SourceMode::kAdditive)) &&
        (mParameters.getPressureSourceFlag()  ||
         mParameters.getVelocityXSourceFlag() ||
         mParameters.getVelocityYSourceFlag() ||
         mParameters.getVelocityZSourceFlag()))
    {
      generateSourceKappa();
    }

    // calculate c^2. It has to be after kappa gen... because of c modification
    computeC2();
  }

  Logger::log(Logger::LogLevel::kBasic, kOutFmtDone);
}// end of preProcessing
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute the main time loop of KSpaceFirstOrder solver.
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::computeMainLoop()
{
  mActPercent = 0;
  // set ActPercent to correspond the time index after recovery
  if (mParameters.getTimeIndex() > 0)
  {
    mActPercent = (100 * mParameters.getTimeIndex()) / mParameters.getNt();
  }

  // Progress header
  Logger::log(Logger::LogLevel::kBasic,kOutFmtSimulationHeader);

  mIterationTime.start();

  // execute main loop
  while ((mParameters.getTimeIndex() < mParameters.getNt()) &&
         (!mParameters.isTimeToCheckpoint(mTotalTime)))
  {
    const size_t timeIndex = mParameters.getTimeIndex();
    // compute velocity
    computeVelocity<simulationDimension>();
    // add in the velocity source term
    addVelocitySource();

    // add in the transducer source term (t = t1) to ux
    if (mParameters.getTransducerSourceFlag() > timeIndex)
    {
      // transducer source is added only to the x component of the particle velocity
      addTransducerSource();
    }

    // compute gradient of velocity
    computeVelocityGradient<simulationDimension>();

    if (mParameters.getNonLinearFlag())
    {
      computeDensityNonliner<simulationDimension>();
    }
    else
    {
      computeDensityLinear<simulationDimension>();
    }


     // add in the source pressure term
     addPressureSource<simulationDimension>();

    if (mParameters.getNonLinearFlag())
    {
      computePressureNonlinear<simulationDimension>();
    }
    else
    {
      computePressureLinear<simulationDimension>();
    }

    // calculate initial pressure
    if ((timeIndex == 0) && (mParameters.getInitialPressureSourceFlag() == 1))
    {
      addInitialPressureSource<simulationDimension>();
    }

    storeSensorData();
    printStatistics();
    mParameters.incrementTimeIndex();
  }// time loop
}// end of computeMainLoop
//----------------------------------------------------------------------------------------------------------------------

/**
 * Post processing the quantities, closing the output streams and storing the sensor mask.
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::postProcessing()
{
  if (mParameters.getStorePressureFinalAllFlag())
  {
    getP().writeData(mParameters.getOutputFile(), kPressureFinalName, mParameters.getCompressionLevel());
  }// p_final

  if (mParameters.getStoreVelocityFinalAllFlag())
  {
    getUxSgx().writeData(mParameters.getOutputFile(), kUxFinalName, mParameters.getCompressionLevel());
    getUySgy().writeData(mParameters.getOutputFile(), kUyFinalName, mParameters.getCompressionLevel());
    if (mParameters.isSimulation3D())
    {
      getUzSgz().writeData(mParameters.getOutputFile(), kUzFinalName, mParameters.getCompressionLevel());
    }
  }// u_final

  if (mParameters.getStoreIntensityAvgFlag() ||
      mParameters.getStoreIntensityAvgCFlag() ||
      mParameters.getStoreQTermFlag() ||
      mParameters.getStoreQTermCFlag())
  {
    Logger::log(Logger::LogLevel::kBasic, kOutFmtNoDone);
  }

  // Compute average intensity from stored p and u without compression
  if (mParameters.getStoreQTermFlag() || mParameters.getStoreIntensityAvgFlag())
  {
    Logger::log(Logger::LogLevel::kBasic, kOutFmtComputingAverageIntensity);
    computeAverageIntensities<simulationDimension>();
    Logger::log(Logger::LogLevel::kBasic, kOutFmtDone);
  }

  // Compute average intensity from stored p and u compression coefficients
  if (mParameters.getOnlyPostProcessingFlag() && (mParameters.getStoreQTermCFlag() || mParameters.getStoreIntensityAvgCFlag()))
  {
    Logger::log(Logger::LogLevel::kBasic, kOutFmtComputingAverageIntensityC);
    computeAverageIntensitiesC<simulationDimension>();
    Logger::log(Logger::LogLevel::kBasic, kOutFmtDone);
  }

  // Apply post-processing
  mOutputStreamContainer.postProcessStreams();

  // Compute and store Q term (volume rate of heat deposition) from average intensity
  if (mParameters.getStoreQTermFlag())
  {
    Logger::log(Logger::LogLevel::kBasic, kOutFmtComputingQTerm);
    computeQTerm<simulationDimension>(OutputStreamContainer::OutputStreamIdx::kIntensityXAvg,
                                      OutputStreamContainer::OutputStreamIdx::kIntensityYAvg,
                                      OutputStreamContainer::OutputStreamIdx::kIntensityZAvg,
                                      OutputStreamContainer::OutputStreamIdx::kQTerm);
    Logger::log(Logger::LogLevel::kBasic, kOutFmtDone);
  }

  // Compute and store Q term (volume rate of heat deposition) from average intensity computed using compression
  if (mParameters.getStoreQTermCFlag())
  {
    Logger::log(Logger::LogLevel::kBasic, kOutFmtComputingQTermC);
    computeQTerm<simulationDimension>(OutputStreamContainer::OutputStreamIdx::kIntensityXAvgC,
                                      OutputStreamContainer::OutputStreamIdx::kIntensityYAvgC,
                                      OutputStreamContainer::OutputStreamIdx::kIntensityZAvgC,
                                      OutputStreamContainer::OutputStreamIdx::kQTermC);
    Logger::log(Logger::LogLevel::kBasic, kOutFmtDone);
  }

  if (mParameters.getStoreIntensityAvgFlag() ||
      mParameters.getStoreIntensityAvgCFlag() ||
      mParameters.getStoreQTermFlag() ||
      mParameters.getStoreQTermCFlag())
  {
    Logger::log(Logger::LogLevel::kBasic, kOutFmtEmpty);
  }

  // Apply post-processing 2
  mOutputStreamContainer.postProcessStreams2();

  // Close
  mOutputStreamContainer.closeStreams();

  // store sensor mask if wanted
  if (mParameters.getCopySensorMaskFlag())
  {
    if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kIndex
        && !mParameters.getOutputFile().datasetExists(mParameters.getOutputFile().getRootGroup(), kSensorMaskIndexName))
    {
      getSensorMaskIndex().recomputeIndicesToMatlab();
      getSensorMaskIndex().writeData(mParameters.getOutputFile(), kSensorMaskIndexName,
                                     mParameters.getCompressionLevel());
    }
    if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kCorners
        && !mParameters.getOutputFile().datasetExists(mParameters.getOutputFile().getRootGroup(), kSensorMaskCornersName))
    {
      getSensorMaskCorners().recomputeIndicesToMatlab();
      getSensorMaskCorners().writeData(mParameters.getOutputFile(), kSensorMaskCornersName,
                                       mParameters.getCompressionLevel());
    }
  }
}// end of postProcessing
//----------------------------------------------------------------------------------------------------------------------

/**
 * Store sensor data.
 */
void KSpaceFirstOrderSolver::storeSensorData()
{
  // Unless the time for sampling has come, exit
  if (mParameters.getTimeIndex() >= mParameters.getSamplingStartTimeIndex())
  {
    if (mParameters.getStoreVelocityNonStaggeredRawFlag()
        || mParameters.getStoreVelocityNonStaggeredCFlag()
        || mParameters.getStoreIntensityAvgFlag()
        || mParameters.getStoreQTermFlag()
        || mParameters.getStoreIntensityAvgCFlag()
        || mParameters.getStoreQTermCFlag())
    {
      if (mParameters.isSimulation3D())
      {
        computeShiftedVelocity<SD::k3D>();
      }
      else
      {
        computeShiftedVelocity<SD::k2D>();
      }
    }
    mOutputStreamContainer.sampleStreams();
  }
}// end of storeSensorData
//----------------------------------------------------------------------------------------------------------------------

/**
 * Write statistics and the header into the output file.
 */
void KSpaceFirstOrderSolver::writeOutputDataInfo()
{
  // write timeIndex into the output file
  mParameters.getOutputFile().writeScalarValue(mParameters.getOutputFile().getRootGroup(),
                                               kTimeIndexName,
                                               mParameters.getTimeIndex());

  // Write scalars
  mParameters.saveScalarsToOutputFile();
  Hdf5FileHeader& fileHeader = mParameters.getFileHeader();

  // Write File header
  fileHeader.setCodeName(getCodeName());
  fileHeader.setMajorFileVersion();
  fileHeader.setMinorFileVersion();
  fileHeader.setActualCreationTime();
  fileHeader.setFileType(Hdf5FileHeader::FileType::kOutput);
  fileHeader.setHostName();

  fileHeader.setMemoryConsumption(getMemoryUsage());

  // Stop total timer here
  mTotalTime.stop();
  fileHeader.setExecutionTimes(getCumulatedTotalTime(),
                               getCumulatedDataLoadTime(),
                               getCumulatedPreProcessingTime(),
                               getCumulatedSimulationTime(),
                               getCumulatedPostProcessingTime());

  fileHeader.setNumberOfCores();
  fileHeader.writeHeaderToOutputFile(mParameters.getOutputFile());

  // Write some other info
  hid_t rootGroup = mParameters.getOutputFile().getRootGroup();
  mParameters.getOutputFile().writeStringAttribute(rootGroup,  "/", "mStoreQTermFlag", std::to_string(mParameters.getStoreQTermFlag()));
  mParameters.getOutputFile().writeStringAttribute(rootGroup,  "/", "mStoreQTermCFlag", std::to_string(mParameters.getStoreQTermCFlag()));
  mParameters.getOutputFile().writeStringAttribute(rootGroup,  "/", "mStoreIntensityAvgFlag", std::to_string(mParameters.getStoreIntensityAvgFlag()));
  mParameters.getOutputFile().writeStringAttribute(rootGroup,  "/", "mStoreIntensityAvgCFlag", std::to_string(mParameters.getStoreIntensityAvgCFlag()));
  mParameters.getOutputFile().writeStringAttribute(rootGroup,  "/", "mNoCompressionOverlapFlag", std::to_string(mParameters.getNoCompressionOverlapFlag()));
  mParameters.getOutputFile().writeStringAttribute(rootGroup,  "/", "mPeriod", std::to_string((&CompressHelper::getInstance())->getPeriod()));
  mParameters.getOutputFile().writeStringAttribute(rootGroup,  "/", "mMOS", std::to_string((&CompressHelper::getInstance())->getMos()));
  mParameters.getOutputFile().writeStringAttribute(rootGroup,  "/", "mHarmonics", std::to_string((&CompressHelper::getInstance())->getHarmonics()));
  mParameters.getOutputFile().writeStringAttribute(rootGroup,  "/", "mBlockSizeDefault", std::to_string(mBlockSizeDefault));
  mParameters.getOutputFile().writeStringAttribute(rootGroup,  "/", "mBlockSizeDefaultC", std::to_string(mBlockSizeDefaultC));
  mParameters.getOutputFile().writeStringAttribute(rootGroup,  "/", "mSamplingStartTimeStep", std::to_string(mParameters.getSamplingStartTimeIndex()));
  mParameters.getOutputFile().writeStringAttribute(rootGroup,  "/", "output_file_size", std::to_string(mParameters.getOutputFile().getFileSize()));


}// end of writeOutputDataInfo
//----------------------------------------------------------------------------------------------------------------------

/**
 * Save checkpoint data into the checkpoint file, flush aggregated outputs into the output file.
 */
void KSpaceFirstOrderSolver::saveCheckpointData()
{
  #if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
     Logger::log(Logger::LogLevel::kFull, kOutFmtStoringFftwWisdom);
     Logger::flush(Logger::LogLevel::kFull);
    // export FFTW wisdom
     try
     {
       FftwComplexMatrix::exportWisdom();
       Logger::log(Logger::LogLevel::kFull, kOutFmtDone);
     }
     catch (const std::runtime_error& e)
     {
       Logger::log(Logger::LogLevel::kFull, kOutFmtFailed);
     }
  #endif

  // Create Checkpoint file
  Hdf5File& checkpointFile = mParameters.getCheckpointFile();
  // if it happens and the file is opened (from the recovery, close it)
  if (checkpointFile.isOpen()) checkpointFile.close();

  Logger::log(Logger::LogLevel::kFull, kOutFmtStoringCheckpointData);
  Logger::flush(Logger::LogLevel::kFull);

  // Create the new file (overwrite the old one)
  checkpointFile.create(mParameters.getCheckpointFileName());


  //------------------------------------------------ Store Matrices --------------------------------------------------//
  // Store all necessary matrices in Checkpoint file
  mMatrixContainer.storeDataIntoCheckpointFile();

  // Write t_index
  checkpointFile.writeScalarValue(checkpointFile.getRootGroup(), kTimeIndexName, mParameters.getTimeIndex());

  // store basic dimension sizes (Nx, Ny, Nz) - Nt is not necessary
  checkpointFile.writeScalarValue(checkpointFile.getRootGroup(), kNxName, mParameters.getFullDimensionSizes().nx);
  checkpointFile.writeScalarValue(checkpointFile.getRootGroup(), kNyName, mParameters.getFullDimensionSizes().ny);
  checkpointFile.writeScalarValue(checkpointFile.getRootGroup(), kNzName, mParameters.getFullDimensionSizes().nz);


  // Write checkpoint file header
  Hdf5FileHeader fileHeader = mParameters.getFileHeader();

  fileHeader.setFileType(Hdf5FileHeader::FileType::kCheckpoint);
  fileHeader.setCodeName(getCodeName());
  fileHeader.setActualCreationTime();

  fileHeader.writeHeaderToCheckpointFile(checkpointFile);

  Logger::log(Logger::LogLevel::kFull, kOutFmtDone);

  // checkpoint output streams only if necessary (t_index > start_index) - here we're at  step + 1
  if (mParameters.getTimeIndex() > mParameters.getSamplingStartTimeIndex())
  {
    Logger::log(Logger::LogLevel::kFull,kOutFmtStoringSensorData);
    Logger::flush(Logger::LogLevel::kFull);
    mOutputStreamContainer.checkpointStreams();
    Logger::log(Logger::LogLevel::kFull, kOutFmtDone);
  }
  mOutputStreamContainer.closeStreams();
  // Close the checkpoint file
  checkpointFile.close();
}// end of saveCheckpointData
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute average intensities from stored p and u without compression.
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::computeAverageIntensities()
{
  float* intensityXAvgData = mOutputStreamContainer[OutputStreamContainer::OutputStreamIdx::kIntensityXAvg].getCurrentStoreBuffer();
  float* intensityYAvgData = mOutputStreamContainer[OutputStreamContainer::OutputStreamIdx::kIntensityYAvg].getCurrentStoreBuffer();
  float* intensityZAvgData = (simulationDimension == SD::k3D) ? mOutputStreamContainer[OutputStreamContainer::OutputStreamIdx::kIntensityZAvg].getCurrentStoreBuffer() : nullptr;
  size_t steps = 0;
  if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kIndex)
  {
    steps = mParameters.getOutputFile().getDatasetDimensionSizes(mParameters.getOutputFile().getRootGroup(), kPName).ny;
  }
  else if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kCorners)
  {
    steps = mParameters.getOutputFile().getDatasetDimensionSizes(mParameters.getOutputFile().getRootGroup(), kPName + "/1").nt;
  }

  // Compute shifts
  const float pi2 = static_cast<float>(M_PI) * 2.0f;
  const size_t stepsComplex = steps / 2 + 1;
  FloatComplex* kx = reinterpret_cast<FloatComplex*>(_mm_malloc(stepsComplex * sizeof(FloatComplex), kDataAlignment));
  //#pragma omp simd - Intel exp bug
  for (size_t i = 0; i < stepsComplex; i++)
  {
    const ssize_t shift = ssize_t((i + (steps / 2)) % steps - (steps / 2));
    kx[i] = exp(FloatComplex(0.0f, 1.0f) * (pi2 * 0.5f) * (float(shift) / float(steps)));
  }

  size_t sensorMaskSize = 0;
  size_t numberOfDatasets = 1;

  if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kIndex)
  {
    sensorMaskSize = getSensorMaskIndex().capacity();
  }
  else if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kCorners)
  {
    sensorMaskSize = getSensorMaskCorners().getSizeOfAllCuboids();
    numberOfDatasets = getSensorMaskCorners().getDimensionSizes().ny;
  }

  hid_t datasetP = 0;
  hid_t datasetUx = 0;
  hid_t datasetUy = 0;
  hid_t datasetUz = 0;

  size_t maxBlockSize = mParameters.getBlockSize();
  DimensionSizes fullDims = mParameters.getFullDimensionSizes();

  // Compute max block size for dataset reading
  if (maxBlockSize == 0)
  {
    size_t memory = getAvailableMemory() * 1000000;
    //std::cout << std::endl;
    //std::cout << "getMemoryUsage:        " << getMemoryUsage() << std::endl;
    //std::cout << "getTotalMemory:        " << getTotalMemory() << std::endl;
    //std::cout << "getCurrentMemoryUsage: " << getCurrentMemoryUsage() << std::endl;
    //std::cout << "getPeakMemoryUsage:    " << getPeakMemoryUsage() << std::endl;
    //std::cout << "getAvailableMemory:    " << getAvailableMemory() << std::endl;
    // dataP, dataUx, dataUy, dataUz, fftwTimeShiftMatrix -> 5 matrices x 4 (size of float)
    maxBlockSize = size_t(float(memory) / 20 * 0.98f);
  }

  size_t sliceSize = fullDims.nx * fullDims.ny;
  size_t fullSize = fullDims.nx * fullDims.ny * fullDims.nz;
  mBlockSizeDefault = maxBlockSize / steps;
  // Minimal size is sliceSize
  mBlockSizeDefault = mBlockSizeDefault < sliceSize ? sliceSize : mBlockSizeDefault;
  // Maximal size is fullSize
  mBlockSizeDefault = mBlockSizeDefault > fullSize ? fullSize : mBlockSizeDefault;
  DimensionSizes shiftDims(mBlockSizeDefault, steps / 2 + 1, 1);
  mBlockSizeDefault *= steps;
  // Temporary matrix for time shift
  FftwComplexMatrix* fftwTimeShiftMatrix = new FftwComplexMatrix(shiftDims);

  //std::cout << "maxBlockSize:          " << maxBlockSize <<  std::endl;
  //std::cout << "mBlockSizeDefault:     " << mBlockSizeDefault << std::endl;
  //std::cout << std::endl;

  // For every dataset
  // TODO test with more cuboids
  for (size_t indexOfDataset = 1; indexOfDataset <= numberOfDatasets; indexOfDataset++)
  {
    size_t datasetSize = 0;
    DimensionSizes datasetDimensionSizes;
    size_t sliceSize = 0;
    // Block size for given dataset
    size_t blockSize = fftwTimeShiftMatrix->getDimensionSizes().nx;
    DimensionSizes datasetBlockSizes;

    // Open datasets
    if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kIndex)
    {
      datasetP = mParameters.getOutputFile().openDataset(mParameters.getOutputFile().getRootGroup(), kPName);
      datasetUx = mParameters.getOutputFile().openDataset(mParameters.getOutputFile().getRootGroup(), kUxNonStaggeredName);
      datasetUy = mParameters.getOutputFile().openDataset(mParameters.getOutputFile().getRootGroup(), kUyNonStaggeredName);
      if (simulationDimension == SD::k3D)
      {
        datasetUz = mParameters.getOutputFile().openDataset(mParameters.getOutputFile().getRootGroup(), kUzNonStaggeredName);
      }
      datasetSize = mParameters.getOutputFile().getDatasetSize(mParameters.getOutputFile().getRootGroup(), kPName) / steps;
      // Maximum size is datasetSize
      if (blockSize > datasetSize)
      {
        blockSize = datasetSize;
      }
      datasetBlockSizes = DimensionSizes(blockSize, steps, 1);
    }
    else if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kCorners)
    {
      datasetP = mParameters.getOutputFile().openDataset(mParameters.getOutputFile().getRootGroup(), kPName + "/" + std::to_string(indexOfDataset));
      datasetUx = mParameters.getOutputFile().openDataset(mParameters.getOutputFile().getRootGroup(), kUxNonStaggeredName + "/" + std::to_string(indexOfDataset));
      datasetUy = mParameters.getOutputFile().openDataset(mParameters.getOutputFile().getRootGroup(), kUyNonStaggeredName + "/" + std::to_string(indexOfDataset));
      if (simulationDimension == SD::k3D)
      {
        datasetUz = mParameters.getOutputFile().openDataset(mParameters.getOutputFile().getRootGroup(), kUzNonStaggeredName + "/" + std::to_string(indexOfDataset));
      }
      datasetSize = mParameters.getOutputFile().getDatasetSize(mParameters.getOutputFile().getRootGroup(), kPName + "/" + std::to_string(indexOfDataset)) / steps;
      datasetDimensionSizes = mParameters.getOutputFile().getDatasetDimensionSizes(mParameters.getOutputFile().getRootGroup(), kPName + "/" + std::to_string(indexOfDataset));
      sliceSize = datasetDimensionSizes.nx * datasetDimensionSizes.ny;
      // Maximum size is datasetSize
      if (blockSize > datasetSize)
      {
        blockSize = datasetSize;
      }
      size_t zCount = blockSize / sliceSize;
      blockSize = sliceSize * zCount;
      datasetBlockSizes = DimensionSizes(datasetDimensionSizes.nx, datasetDimensionSizes.ny, zCount, steps);
    }

    Logger::log(Logger::LogLevel::kBasic, kOutFmtNoDone);
    Logger::log(Logger::LogLevel::kBasic, kOutFmtBlockSizePostProcessing, (Logger::formatMessage(kOutFmt2DDomainSizeFormat, blockSize, steps)).c_str(), (blockSize * steps * 4) / 1000000);
    Logger::log(Logger::LogLevel::kBasic, kOutFmtEmpty);

    RealMatrix dataP(DimensionSizes(blockSize, steps, 1));
    RealMatrix dataUx(DimensionSizes(blockSize, steps, 1));
    RealMatrix dataUy(DimensionSizes(blockSize, steps, 1));
    RealMatrix dataUz(DimensionSizes(blockSize, steps, 1));
    fftwTimeShiftMatrix->createR2CFftPlan1DY(dataP);
    fftwTimeShiftMatrix->createC2RFftPlan1DY(dataP);

    for (size_t i = 0; i < datasetSize; i += blockSize)
    {
      // Read block by sensor mask type
      if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kIndex)
      {
        if (i + blockSize > datasetSize)
        {
          blockSize = datasetSize - i;
          datasetBlockSizes = DimensionSizes(blockSize, steps, 1);
          dataP.resize(DimensionSizes(blockSize, steps, 1));
          dataUx.resize(DimensionSizes(blockSize, steps, 1));
          dataUy.resize(DimensionSizes(blockSize, steps, 1));
          dataUz.resize(DimensionSizes(blockSize, steps, 1));
          fftwTimeShiftMatrix->createR2CFftPlan1DY(dataP);
          fftwTimeShiftMatrix->createC2RFftPlan1DY(dataP);
        }

        mParameters.getOutputFile().readHyperSlab(datasetP, DimensionSizes(i, 0, 0), datasetBlockSizes, dataP.getData());
        mParameters.getOutputFile().readHyperSlab(datasetUx, DimensionSizes(i, 0, 0), datasetBlockSizes, dataUx.getData());
        mParameters.getOutputFile().readHyperSlab(datasetUy, DimensionSizes(i, 0, 0), datasetBlockSizes, dataUy.getData());
        if (simulationDimension == SD::k3D)
        {
          mParameters.getOutputFile().readHyperSlab(datasetUz, DimensionSizes(i, 0, 0), datasetBlockSizes, dataUz.getData());
        }
      }
      else if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kCorners)
      {
        if (i + blockSize > datasetSize)
        {
          blockSize = datasetSize - i;
          datasetBlockSizes = DimensionSizes(datasetDimensionSizes.nx, datasetDimensionSizes.ny, blockSize / sliceSize, steps);
          dataP.resize(DimensionSizes(blockSize, steps, 1));
          dataUx.resize(DimensionSizes(blockSize, steps, 1));
          dataUy.resize(DimensionSizes(blockSize, steps, 1));
          dataUz.resize(DimensionSizes(blockSize, steps, 1));
          fftwTimeShiftMatrix->createR2CFftPlan1DY(dataP);
          fftwTimeShiftMatrix->createC2RFftPlan1DY(dataP);
        }
        size_t zOffset = i / sliceSize;

        mParameters.getOutputFile().readHyperSlab(datasetP, DimensionSizes(0, 0, zOffset, 0), datasetBlockSizes, dataP.getData());
        mParameters.getOutputFile().readHyperSlab(datasetUx, DimensionSizes(0, 0, zOffset, 0), datasetBlockSizes, dataUx.getData());
        mParameters.getOutputFile().readHyperSlab(datasetUy, DimensionSizes(0, 0, zOffset, 0), datasetBlockSizes, dataUy.getData());
        if (simulationDimension == SD::k3D)
        {
          mParameters.getOutputFile().readHyperSlab(datasetUz, DimensionSizes(0, 0, zOffset, 0), datasetBlockSizes, dataUz.getData());
        }
      }

      // Phase shifts using FFT
      const float divider = 1.0f / static_cast<float>(steps);
      FloatComplex* tempFftTimeShift = fftwTimeShiftMatrix->getComplexData();
      fftwTimeShiftMatrix->computeR2CFft1DY(dataUx);
      #pragma omp simd
      for (size_t step = 0; step < stepsComplex; step++)
      {
        for (size_t x = 0; x < blockSize; x++)
        {
          tempFftTimeShift[step * blockSize + x] *= divider * kx[step];
        }
      }
      fftwTimeShiftMatrix->computeC2RFft1DY(dataUx);
      fftwTimeShiftMatrix->computeR2CFft1DY(dataUy);
      #pragma omp simd
      for (size_t step = 0; step < stepsComplex; step++)
      {
        for (size_t x = 0; x < blockSize; x++)
        {
          tempFftTimeShift[step * blockSize + x] *= divider * kx[step];
        }
      }
      fftwTimeShiftMatrix->computeC2RFft1DY(dataUy);
      if (simulationDimension == SD::k3D)
      {
        fftwTimeShiftMatrix->computeR2CFft1DY(dataUz);
        #pragma omp simd
        for (size_t step = 0; step < stepsComplex; step++)
        {
          for (size_t x = 0; x < blockSize; x++)
          {
            tempFftTimeShift[step * blockSize + x] *= divider * kx[step];
          }
        }
        fftwTimeShiftMatrix->computeC2RFft1DY(dataUz);
      }

      // Compute average of intensity
      for (size_t step = 0; step < steps; step++)
      {
        for (size_t x = 0; x < blockSize; x++)
        {
          intensityXAvgData[i + x] += dataUx[step * blockSize + x] * dataP[step * blockSize + x];
          intensityYAvgData[i + x] += dataUy[step * blockSize + x] * dataP[step * blockSize + x];
          if (simulationDimension == SD::k3D)
          {
            intensityZAvgData[i + x] += dataUz[step * blockSize + x] * dataP[step * blockSize + x];
          }
          if (step == steps - 1)
          {
            intensityXAvgData[i + x] /= steps;
            intensityYAvgData[i + x] /= steps;
            if (simulationDimension == SD::k3D)
            {
              intensityZAvgData[i + x] /= steps;
            }
          }
        }
      }
    }

    mParameters.getOutputFile().closeDataset(datasetP);
    mParameters.getOutputFile().closeDataset(datasetUx);
    mParameters.getOutputFile().closeDataset(datasetUy);
    if (simulationDimension == SD::k3D)
    {
      mParameters.getOutputFile().closeDataset(datasetUz);
    }
  }
  _mm_free(kx);
  delete fftwTimeShiftMatrix;

}// end of computeAverageIntensities
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute average intensities from stored p and u compression coefficients.
 *
 * NOTE does not work with 40-bit complex floats
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::computeAverageIntensitiesC()
{
  float* intensityXAvgData = mOutputStreamContainer[OutputStreamContainer::OutputStreamIdx::kIntensityXAvgC].getCurrentStoreBuffer();
  float* intensityYAvgData = mOutputStreamContainer[OutputStreamContainer::OutputStreamIdx::kIntensityYAvgC].getCurrentStoreBuffer();
  float* intensityZAvgData = (simulationDimension == SD::k3D) ? mOutputStreamContainer[OutputStreamContainer::OutputStreamIdx::kIntensityZAvgC].getCurrentStoreBuffer() : nullptr;
  size_t steps = 1;
  if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kIndex)
  {
    steps = mParameters.getOutputFile().getDatasetDimensionSizes(mParameters.getOutputFile().getRootGroup(), kPName + kCompressSuffix).ny;
  }
  else if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kCorners)
  {
    steps = mParameters.getOutputFile().getDatasetDimensionSizes(mParameters.getOutputFile().getRootGroup(), kPName + kCompressSuffix + "/1").nt;
  }

  size_t sensorMaskSize = 0;
  size_t numberOfDatasets = 1;

  if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kIndex)
  {
    sensorMaskSize = getSensorMaskIndex().capacity();
  }
  else if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kCorners)
  {
    sensorMaskSize = getSensorMaskCorners().getSizeOfAllCuboids();
    numberOfDatasets = getSensorMaskCorners().getDimensionSizes().ny;
  }

  hid_t datasetP = 0;
  hid_t datasetUx = 0;
  hid_t datasetUy = 0;
  hid_t datasetUz = 0;

  size_t maxBlockSize = mParameters.getBlockSize();
  DimensionSizes fullDims = mParameters.getFullDimensionSizes();
  fullDims.nx *= CompressHelper::getInstance().getHarmonics() * 2; // TODO size of complex (40-bit complex floats)

  // Compute max block size for dataset reading
  if (maxBlockSize == 0)
  {
    size_t memory = getAvailableMemory() * 1000000;
    // dataP, dataUx, dataUy, dataUz, fftwTimeShiftMatrix -> 5 matrices x 4 (size of float)
    maxBlockSize = size_t(float(memory) / 20 * 0.98f);
  }
  size_t sliceSize = fullDims.nx * fullDims.ny;
  size_t fullSize = fullDims.nx * fullDims.ny * fullDims.nz;
  mBlockSizeDefaultC = maxBlockSize / steps;
  // Minimal size is sliceSize
  mBlockSizeDefaultC = mBlockSizeDefaultC < sliceSize ? sliceSize : mBlockSizeDefaultC;
  // Maximal size is fullSize
  mBlockSizeDefaultC = mBlockSizeDefaultC > fullSize ? fullSize : mBlockSizeDefaultC;
  mBlockSizeDefault *= steps;

  // For every dataset
  // TODO test with more cuboids
  for (size_t indexOfDataset = 1; indexOfDataset <= numberOfDatasets; indexOfDataset++)
  {
    size_t datasetSize = 0;
    DimensionSizes datasetDimensionSizes;
    size_t sliceSize = 0;
    // Block size for given dataset
    size_t blockSize = mBlockSizeDefaultC / steps;
    DimensionSizes datasetBlockSizes;

    // Open datasets
    if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kIndex)
    {
      datasetP = mParameters.getOutputFile().openDataset(mParameters.getOutputFile().getRootGroup(), kPName + kCompressSuffix);
      datasetUx = mParameters.getOutputFile().openDataset(mParameters.getOutputFile().getRootGroup(), kUxNonStaggeredName + kCompressSuffix);
      datasetUy = mParameters.getOutputFile().openDataset(mParameters.getOutputFile().getRootGroup(), kUyNonStaggeredName + kCompressSuffix);
      if (simulationDimension == SD::k3D)
      {
        datasetUz = mParameters.getOutputFile().openDataset(mParameters.getOutputFile().getRootGroup(), kUzNonStaggeredName + kCompressSuffix);
      }
      datasetSize = mParameters.getOutputFile().getDatasetSize(mParameters.getOutputFile().getRootGroup(), kPName + kCompressSuffix) / steps;
      // Maximum size is datasetSize
      if (blockSize > datasetSize)
      {
        blockSize = datasetSize;
      }
      datasetBlockSizes = DimensionSizes(blockSize, steps, 1);
    }
    else if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kCorners)
    {
      datasetP = mParameters.getOutputFile().openDataset(mParameters.getOutputFile().getRootGroup(), kPName + kCompressSuffix + "/" + std::to_string(indexOfDataset));
      datasetUx = mParameters.getOutputFile().openDataset(mParameters.getOutputFile().getRootGroup(), kUxNonStaggeredName + kCompressSuffix + "/" + std::to_string(indexOfDataset));
      datasetUy = mParameters.getOutputFile().openDataset(mParameters.getOutputFile().getRootGroup(), kUyNonStaggeredName + kCompressSuffix + "/" + std::to_string(indexOfDataset));
      if (simulationDimension == SD::k3D)
      {
        datasetUz = mParameters.getOutputFile().openDataset(mParameters.getOutputFile().getRootGroup(), kUzNonStaggeredName + kCompressSuffix + "/" + std::to_string(indexOfDataset));
      }
      datasetSize = mParameters.getOutputFile().getDatasetSize(mParameters.getOutputFile().getRootGroup(), kPName + kCompressSuffix + "/" + std::to_string(indexOfDataset)) / steps;
      datasetDimensionSizes = mParameters.getOutputFile().getDatasetDimensionSizes(mParameters.getOutputFile().getRootGroup(), kPName + kCompressSuffix + "/" + std::to_string(indexOfDataset));
      sliceSize = datasetDimensionSizes.nx * datasetDimensionSizes.ny;
      // Maximum size is datasetSize
      if (blockSize > datasetSize)
      {
        blockSize = datasetSize;
      }
      size_t zCount = blockSize / sliceSize;
      blockSize = sliceSize * zCount;
      datasetBlockSizes = DimensionSizes(datasetDimensionSizes.nx, datasetDimensionSizes.ny, zCount, steps);
    }

    Logger::log(Logger::LogLevel::kBasic, kOutFmtNoDone);
    Logger::log(Logger::LogLevel::kBasic, kOutFmtBlockSizePostProcessing, (Logger::formatMessage(kOutFmt2DDomainSizeFormat, blockSize, steps)).c_str(), (blockSize * steps * 4) / 1000000);
    Logger::log(Logger::LogLevel::kBasic, kOutFmtEmpty);

    RealMatrix dataP(DimensionSizes(blockSize, steps, 1));
    RealMatrix dataUx(DimensionSizes(blockSize, steps, 1));
    RealMatrix dataUy(DimensionSizes(blockSize, steps, 1));
    RealMatrix dataUz(DimensionSizes(blockSize, steps, 1));

    for (size_t i = 0; i < datasetSize; i += blockSize)
    {
      const size_t outputIndex = i / (CompressHelper::getInstance().getHarmonics() * 2); // TODO size of complex (40-bit complex floats)

      // Read block by sensor mask type
      if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kIndex)
      {
        if (i + blockSize > datasetSize)
        {
          blockSize = datasetSize - i;
          datasetBlockSizes = DimensionSizes(blockSize, steps, 1);
          dataP.resize(DimensionSizes(blockSize, steps, 1));
          dataUx.resize(DimensionSizes(blockSize, steps, 1));
          dataUy.resize(DimensionSizes(blockSize, steps, 1));
          dataUz.resize(DimensionSizes(blockSize, steps, 1));
        }

        mParameters.getOutputFile().readHyperSlab(datasetP, DimensionSizes(i, 0, 0), datasetBlockSizes, dataP.getData());
        mParameters.getOutputFile().readHyperSlab(datasetUx, DimensionSizes(i, 0, 0), datasetBlockSizes, dataUx.getData());
        mParameters.getOutputFile().readHyperSlab(datasetUy, DimensionSizes(i, 0, 0), datasetBlockSizes, dataUy.getData());
        if (simulationDimension == SD::k3D)
        {
          mParameters.getOutputFile().readHyperSlab(datasetUz, DimensionSizes(i, 0, 0), datasetBlockSizes, dataUz.getData());
        }
      }
      else if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kCorners)
      {
        if (i + blockSize > datasetSize)
        {
          blockSize = datasetSize - i;
          datasetBlockSizes = DimensionSizes(datasetDimensionSizes.nx, datasetDimensionSizes.ny, blockSize / sliceSize, steps);
          dataP.resize(DimensionSizes(blockSize, steps, 1));
          dataUx.resize(DimensionSizes(blockSize, steps, 1));
          dataUy.resize(DimensionSizes(blockSize, steps, 1));
          dataUz.resize(DimensionSizes(blockSize, steps, 1));
        }
        size_t zOffset = i / sliceSize;

        mParameters.getOutputFile().readHyperSlab(datasetP, DimensionSizes(0, 0, zOffset, 0), datasetBlockSizes, dataP.getData());
        mParameters.getOutputFile().readHyperSlab(datasetUx, DimensionSizes(0, 0, zOffset, 0), datasetBlockSizes, dataUx.getData());
        mParameters.getOutputFile().readHyperSlab(datasetUy, DimensionSizes(0, 0, zOffset, 0), datasetBlockSizes, dataUy.getData());
        if (simulationDimension == SD::k3D)
        {
          mParameters.getOutputFile().readHyperSlab(datasetUz, DimensionSizes(0, 0, zOffset, 0), datasetBlockSizes, dataUz.getData());
        }
      }

      FloatComplex* bufferP = reinterpret_cast<FloatComplex*>(dataP.getData());
      FloatComplex* bufferUx = reinterpret_cast<FloatComplex*>(dataUx.getData());
      FloatComplex* bufferUy = reinterpret_cast<FloatComplex*>(dataUy.getData());
      FloatComplex* bufferUz = reinterpret_cast<FloatComplex*>(dataUz.getData());

      size_t outputBlockSize = blockSize / (CompressHelper::getInstance().getHarmonics() * 2); // TODO size of complex (40-bit complex floats)
      // Compute average of intensity
      for (size_t step = 0; step < steps; step++)
      {
        for (size_t x = 0; x < outputBlockSize; x++)
        {
          size_t offset = step * (blockSize / 2) + CompressHelper::getInstance().getHarmonics() * x;
          //For every harmonics
          for (size_t ih = 0; ih < CompressHelper::getInstance().getHarmonics(); ih++)
          {
            size_t pH = offset + ih;
            intensityXAvgData[outputIndex + x] += real(bufferP[pH] * conj(bufferUx[pH])) / 2.0f;
            intensityYAvgData[outputIndex + x] += real(bufferP[pH] * conj(bufferUy[pH])) / 2.0f;
            if (simulationDimension == SD::k3D)
            {
              intensityZAvgData[outputIndex + x] += real(bufferP[pH] * conj(bufferUz[pH])) / 2.0f;
            }
          }
          if (step == steps - 1)
          {
            intensityXAvgData[outputIndex + x] /= steps;
            intensityYAvgData[outputIndex + x] /= steps;
            if (simulationDimension == SD::k3D)
            {
              intensityZAvgData[outputIndex + x] /= steps;
            }
          }
        }
      }
    }

    mParameters.getOutputFile().closeDataset(datasetP);
    mParameters.getOutputFile().closeDataset(datasetUx);
    mParameters.getOutputFile().closeDataset(datasetUy);
    if (simulationDimension == SD::k3D)
    {
      mParameters.getOutputFile().closeDataset(datasetUz);
    }
  }
}// end of computeAverageIntensitiesC
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute Q term (volume rate of heat deposition) from average intensities.
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::computeQTerm(OutputStreamContainer::OutputStreamIdx intensityXAvgStreamIndex,
                                          OutputStreamContainer::OutputStreamIdx intensityYAvgStreamIndex,
                                          OutputStreamContainer::OutputStreamIdx intensityZAvgStreamIndex,
                                          OutputStreamContainer::OutputStreamIdx qTermStreamIdx)
{
  float* intensityXAvgData = mOutputStreamContainer[intensityXAvgStreamIndex].getCurrentStoreBuffer();
  float* intensityYAvgData = mOutputStreamContainer[intensityYAvgStreamIndex].getCurrentStoreBuffer();
  float* intensityZAvgData = (simulationDimension == SD::k3D) ? mOutputStreamContainer[intensityZAvgStreamIndex].getCurrentStoreBuffer() : nullptr;

  // Full sized matrices for intensities
  RealMatrix intensityXAvg(mParameters.getFullDimensionSizes());
  RealMatrix intensityYAvg(mParameters.getFullDimensionSizes());
  RealMatrix intensityZAvg(mParameters.getFullDimensionSizes());

  const DimensionSizes& fullDimensionSizes = mParameters.getFullDimensionSizes();

  // Copy values from store buffer to positions defined by sensor mask indices
  if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kIndex)
  {
    const size_t sensorMaskSize = getSensorMaskIndex().capacity();
    #pragma omp parallel for simd schedule(static)
    for (size_t i = 0; i < sensorMaskSize; i++)
    {
      intensityXAvg[getSensorMaskIndex()[i]] = intensityXAvgData[i];
      intensityYAvg[getSensorMaskIndex()[i]] = intensityYAvgData[i];
      if (simulationDimension == SD::k3D)
      {
        intensityZAvg[getSensorMaskIndex()[i]] = intensityZAvgData[i];
      }
    }
  }

  // Copy values from store buffer to positions defined by sensor mask corners
  if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kCorners)
  {
    const size_t slabSize = fullDimensionSizes.ny * fullDimensionSizes.nx;
    const size_t rowSize  = fullDimensionSizes.nx;
    const size_t nCuboids = getSensorMaskCorners().getDimensionSizes().ny;
    size_t cuboidInBufferStart = 0;

    #pragma omp parallel
    for (size_t cuboidIdx = 0; cuboidIdx < nCuboids; cuboidIdx++)
    {
      const DimensionSizes topLeftCorner     = getSensorMaskCorners().getTopLeftCorner(cuboidIdx);
      const DimensionSizes bottomRightCorner = getSensorMaskCorners().getBottomRightCorner(cuboidIdx);

      size_t cuboidSlabSize = (bottomRightCorner.ny - topLeftCorner.ny + 1) *
                              (bottomRightCorner.nx - topLeftCorner.nx + 1);
      size_t cuboidRowSize  = (bottomRightCorner.nx - topLeftCorner.nx + 1);

      DimensionSizes cuboidSize(0, 0, 0, 0);      // Size of the cuboid
      cuboidSize    = bottomRightCorner - topLeftCorner;
      cuboidSize.nt = 1;

      #pragma omp for collapse(3)
      for (size_t z = topLeftCorner.nz; z <= bottomRightCorner.nz; z++)
      {
        for (size_t y = topLeftCorner.ny; y <= bottomRightCorner.ny; y++)
        {
          for (size_t x = topLeftCorner.nx; x <= bottomRightCorner.nx; x++)
          {
            const size_t storeBufferIndex = cuboidInBufferStart +
                                            (z - topLeftCorner.nz) * cuboidSlabSize +
                                            (y - topLeftCorner.ny) * cuboidRowSize  +
                                            (x - topLeftCorner.nx);
            const size_t sourceIndex = z * slabSize + y * rowSize + x;
            intensityXAvg[sourceIndex] = intensityXAvgData[storeBufferIndex];
            intensityYAvg[sourceIndex] = intensityYAvgData[storeBufferIndex];
            if (simulationDimension == SD::k3D)
            {
              intensityZAvg[sourceIndex] = intensityZAvgData[storeBufferIndex];
            }
          }
        }
      }
      // must be done only once
      #pragma omp single
      {
        cuboidInBufferStart += cuboidSize.nElements();
      }
    }
  }

  // Helper complex memory
  FloatComplex* tempFftShift = getTempFftwShift().getComplexData();

  const float pi2 = static_cast<float>(M_PI) * 2.0f;

  // Normalization constants for FFT
  const float dividerX = 1.0f / static_cast<float>(fullDimensionSizes.nx);
  const float dividerY = 1.0f / static_cast<float>(fullDimensionSizes.ny);
  const float dividerZ = 1.0f / static_cast<float>(fullDimensionSizes.nz);
  // Helper values
  const size_t nx = fullDimensionSizes.nx;
  const size_t ny = fullDimensionSizes.ny;
  const size_t nz = fullDimensionSizes.nz;
  const float dx  = mParameters.getDx();
  const float dy  = mParameters.getDy();
  const float dz  = mParameters.getDz();

  DimensionSizes xShiftDims = mParameters.getFullDimensionSizes();
  xShiftDims.nx = xShiftDims.nx / 2 + 1;
  DimensionSizes yShiftDims = mParameters.getFullDimensionSizes();
  yShiftDims.ny = yShiftDims.ny / 2 + 1;
  DimensionSizes zShiftDims = mParameters.getFullDimensionSizes();
  zShiftDims.nz = zShiftDims.nz / 2 + 1;

  // Helper memory for shifts
  FloatComplex* kx = reinterpret_cast<FloatComplex*>(_mm_malloc(xShiftDims.nx * sizeof(FloatComplex), kDataAlignment));
  FloatComplex* ky = reinterpret_cast<FloatComplex*>(_mm_malloc(yShiftDims.ny * sizeof(FloatComplex), kDataAlignment));
  FloatComplex* kz = (simulationDimension == SD::k3D) ? reinterpret_cast<FloatComplex*>(_mm_malloc(zShiftDims.nz * sizeof(FloatComplex), kDataAlignment)) : nullptr;

  // Compute shifts for x gradient
  #pragma omp simd
  for (size_t i = 0; i < xShiftDims.nx; i++)
  {
    const ssize_t shift = ssize_t((i + (nx / 2)) % nx - (nx / 2));
    kx[i] = FloatComplex(0.0f, 1.0f) * (pi2 / dx) * (float(shift) / float(nx));
  }
  // Compute shifts for y gradient
  #pragma omp simd
  for (size_t i = 0; i < yShiftDims.ny; i++)
  {
    const ssize_t shift = ssize_t((i + (ny / 2)) % ny - (ny / 2));
    ky[i] = FloatComplex(0.0f, 1.0f) * (pi2 / dy) * (float(shift) / float(ny));
  }
  if (simulationDimension == SD::k3D)
  {
    // Compute shifts for z gradient
    #pragma omp simd
    for (size_t i = 0; i < zShiftDims.nz; i++)
    {
      const ssize_t shift = ssize_t((i + (nz / 2)) % nz - (nz / 2));
      kz[i] = FloatComplex(0.0f, 1.0f) * (pi2 / dz) * (float(shift) / float(nz));
    }
  }

  getTempFftwShift().computeR2CFft1DX(intensityXAvg);
  #pragma omp parallel for schedule(static) if (simulationDimension == SD::k3D)
  for (size_t z = 0; z < xShiftDims.nz; z++)
  {
    #pragma omp parallel for schedule(static) if (simulationDimension == SD::k2D)
    for (size_t y = 0; y < xShiftDims.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < xShiftDims.nx;  x++)
      {
        const size_t i = get1DIndex(z, y, x, xShiftDims);
        tempFftShift[i] *= dividerX * kx[x];
      } // x
    } // y
  } // z
  getTempFftwShift().computeC2RFft1DX(intensityXAvg);

  getTempFftwShift().computeR2CFft1DY(intensityYAvg);
  #pragma omp parallel for schedule(static) if (simulationDimension == SD::k3D)
  for (size_t z = 0; z < yShiftDims.nz; z++)
  {
    #pragma omp parallel for schedule(static) if (simulationDimension == SD::k2D)
    for (size_t y = 0; y < yShiftDims.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < yShiftDims.nx;  x++)
      {
        const size_t i = get1DIndex(z, y, x, yShiftDims);
        tempFftShift[i] *= dividerY * ky[y];
      } // x
    } // y
  } // z
  getTempFftwShift().computeC2RFft1DY(intensityYAvg);

  if (simulationDimension == SD::k3D)
  {
    getTempFftwShift().computeR2CFft1DZ(intensityZAvg);
    #pragma omp parallel for schedule(static)
    for (size_t z = 0; z < zShiftDims.nz; z++)
    {
      #pragma omp parallel for schedule(static)
      for (size_t y = 0; y < zShiftDims.ny; y++)
      {
        #pragma omp simd
        for (size_t x = 0; x < zShiftDims.nx;  x++)
        {
          const size_t i = get1DIndex(z, y, x, zShiftDims);
          tempFftShift[i] *= dividerZ * kz[z];
        } // x
      } // y
    } // z
    getTempFftwShift().computeC2RFft1DZ(intensityZAvg);
  }

  _mm_free(kx);
  _mm_free(ky);
  _mm_free(kz);

  #pragma omp simd
  for (size_t i = 0; i < fullDimensionSizes.nElements(); i++)
  {
    if (simulationDimension == SD::k3D)
    {
      intensityXAvg[i] = -(intensityXAvg[i] + intensityYAvg[i] + intensityZAvg[i]);
    }
    else
    {
      intensityXAvg[i] = -(intensityXAvg[i] + intensityYAvg[i]);
    }
  }

  float* qTermData = mOutputStreamContainer[qTermStreamIdx].getCurrentStoreBuffer();

  // Copy values from positions defined by sensor mask indices to store buffer
  if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kIndex)
  {
    const size_t sensorMaskSize = getSensorMaskIndex().capacity();
    #pragma omp parallel for simd schedule(static)
    for (size_t i = 0; i < sensorMaskSize; i++)
    {
      qTermData[i] = intensityXAvg[getSensorMaskIndex()[i]];
    }
  }

  // Copy values from positions defined by sensor mask corners to store buffer
  if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kCorners)
  {
    const size_t slabSize = fullDimensionSizes.ny * fullDimensionSizes.nx;
    const size_t rowSize  = fullDimensionSizes.nx;
    const size_t nCuboids = getSensorMaskCorners().getDimensionSizes().ny;
    size_t cuboidInBufferStart = 0;

    #pragma omp parallel
    for (size_t cuboidIdx = 0; cuboidIdx < nCuboids; cuboidIdx++)
    {
      const DimensionSizes topLeftCorner     = getSensorMaskCorners().getTopLeftCorner(cuboidIdx);
      const DimensionSizes bottomRightCorner = getSensorMaskCorners().getBottomRightCorner(cuboidIdx);

      size_t cuboidSlabSize = (bottomRightCorner.ny - topLeftCorner.ny + 1) *
                              (bottomRightCorner.nx - topLeftCorner.nx + 1);
      size_t cuboidRowSize  = (bottomRightCorner.nx - topLeftCorner.nx + 1);

      DimensionSizes cuboidSize(0, 0, 0, 0);      // Size of the cuboid
      cuboidSize    = bottomRightCorner - topLeftCorner;
      cuboidSize.nt = 1;

      #pragma omp for collapse(3)
      for (size_t z = topLeftCorner.nz; z <= bottomRightCorner.nz; z++)
      {
        for (size_t y = topLeftCorner.ny; y <= bottomRightCorner.ny; y++)
        {
          for (size_t x = topLeftCorner.nx; x <= bottomRightCorner.nx; x++)
          {
            const size_t storeBufferIndex = cuboidInBufferStart +
                                            (z - topLeftCorner.nz) * cuboidSlabSize +
                                            (y - topLeftCorner.ny) * cuboidRowSize  +
                                            (x - topLeftCorner.nx);
            const size_t sourceIndex = z * slabSize + y * rowSize + x;
            qTermData[storeBufferIndex] = intensityXAvg[sourceIndex];
          }
        }
      }
      // must be done only once
      #pragma omp single
      {
        cuboidInBufferStart += cuboidSize.nElements();
      }
    }
  }
}// end of computeQTerm
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute new values of acoustic velocity in all used dimensions (UxSgx, UySgy, UzSgz).
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::computeVelocity()
{
  // bsxfun(@times, ddx_k_shift_pos, kappa .* fftn(p)), for all 3 dims
  computePressureGradient<simulationDimension>();

  getTempFftwX().computeC2RFftND(getTemp1RealND());
  getTempFftwY().computeC2RFftND(getTemp2RealND());
  if (simulationDimension == SD::k3D)
  {
    getTempFftwZ().computeC2RFftND(getTemp3RealND());
  }

  if (mParameters.getRho0ScalarFlag())
  { // scalars
    if (mParameters.getNonUniformGridFlag())
    {
      computeVelocityHomogeneousNonuniform<simulationDimension>();
     }
    else
    {
      computeVelocityHomogeneousUniform<simulationDimension>();
    }
  }
  else
  {// matrices
    computeVelocityHeterogeneous<simulationDimension>();
  }
}// end of computeVelocity
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute new values for duxdx, duydy, duzdz.
 */
template<Parameters::SimulationDimension simulationDimension>
void  KSpaceFirstOrderSolver::computeVelocityGradient()
{
  getTempFftwX().computeR2CFftND(getUxSgx());
  getTempFftwY().computeR2CFftND(getUySgy());
  if (simulationDimension == SD::k3D)
  {
    getTempFftwZ().computeR2CFftND(getUzSgz());
  }

  const DimensionSizes& reducedDimensionSizes = mParameters.getReducedDimensionSizes();
  const float divider = 1.0f / static_cast<float>(mParameters.getFullDimensionSizes().nElements());

  const float* kappa = getKappa().getData();

  FloatComplex* tempFftX = getTempFftwX().getComplexData();
  FloatComplex* tempFftY = getTempFftwY().getComplexData();
  FloatComplex* tempFftZ = (simulationDimension == SD::k3D) ? getTempFftwZ().getComplexData() : nullptr;

  FloatComplex* ddxKShiftNeg = getDdxKShiftNeg().getComplexData();
  FloatComplex* ddyKShiftNeg = getDdyKShiftNeg().getComplexData();
  FloatComplex* ddzKShiftNeg = (simulationDimension == SD::k3D) ? getDdzKShiftNeg().getComplexData() : nullptr;

  #pragma omp parallel for schedule(static) if (simulationDimension == SD::k3D)
  for (size_t z = 0; z < reducedDimensionSizes.nz; z++)
  {
    #pragma omp parallel for schedule(static) if (simulationDimension == SD::k2D)
    for (size_t y = 0; y < reducedDimensionSizes.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < reducedDimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, reducedDimensionSizes);
        const float eKappa = divider * kappa[i];

        tempFftX[i] *=  ddxKShiftNeg[x] * eKappa;
        tempFftY[i] *=  ddyKShiftNeg[y] * eKappa;
        if (simulationDimension == SD::k3D)
        {
          tempFftZ[i] *=  ddzKShiftNeg[z] * eKappa;
        }
      } // x
    } // y
  } // z


  getTempFftwX().computeC2RFftND(getDuxdx());
  getTempFftwY().computeC2RFftND(getDuydy());
  if (simulationDimension == SD::k3D)
  {
    getTempFftwZ().computeC2RFftND(getDuzdz());
  }

  //------------------------------------------------ Nonuniform grid -------------------------------------------------//
  if (mParameters.getNonUniformGridFlag() != 0)
  {
    float* duxdx = getDuxdx().getData();
    float* duydy = getDuydy().getData();
    float* duzdz = (simulationDimension == SD::k3D) ? getDuzdz().getData() : nullptr;

    const float* duxdxn = getDxudxn().getData();
    const float* duydyn = getDyudyn().getData();
    const float* duzdzn = (simulationDimension == SD::k3D) ? getDzudzn().getData() : nullptr;

    const DimensionSizes& dimensionSizes = mParameters.getFullDimensionSizes();

    #pragma omp parallel for schedule(static) if (simulationDimension == SD::k3D)
    for (size_t z = 0; z < dimensionSizes.nz; z++)
    {
      #pragma omp parallel for schedule(static) if (simulationDimension == SD::k2D)
      for (size_t y = 0; y < dimensionSizes.ny; y++)
      {
        #pragma omp simd
        for (size_t x = 0; x < dimensionSizes.nx; x++)
        {
          const size_t i = get1DIndex(z, y, x, dimensionSizes);
          duxdx[i] *= duxdxn[x];
          duydy[i] *= duydyn[y];
          if (simulationDimension == SD::k3D)
          {
            duzdz[i] *= duzdzn[z];
          }
        } // x
      } // y
    } // z
 }// nonlinear
}// end of computeVelocityGradient
//----------------------------------------------------------------------------------------------------------------------

/**
 * Calculate new values of acoustic density for nonlinear case (rhoX, rhoy and rhoZ).
 *
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::computeDensityNonliner()
{
  const DimensionSizes& dimensionSizes = mParameters.getFullDimensionSizes();

  const float dt  = mParameters.getDt();

  float* rhoX  = getRhoX().getData();
  float* rhoY  = getRhoY().getData();
  float* rhoZ  = (simulationDimension == SD::k3D) ? getRhoZ().getData() : nullptr;

  const float* pmlX  = getPmlX().getData();
  const float* pmlY  = getPmlY().getData();
  const float* pmlZ  = (simulationDimension == SD::k3D) ? getPmlZ().getData() : nullptr;

  const float* duxdx = getDuxdx().getData();
  const float* duydy = getDuydy().getData();
  const float* duzdz = (simulationDimension == SD::k3D) ? getDuzdz().getData() : nullptr;

  //----------------------------------------------- rho0 is scalar -------------------------------------------------//
  if (mParameters.getRho0ScalarFlag())
  {
    const float rho0 = mParameters.getRho0Scalar();

    #pragma omp parallel for schedule(static) if (simulationDimension == SD::k3D)
    for (size_t z = 0; z < dimensionSizes.nz; z++)
    {
      #pragma omp parallel for schedule(static) if (simulationDimension == SD::k2D)
      for (size_t y = 0; y < dimensionSizes.ny; y++)
      {
        #pragma omp simd
        for (size_t x = 0; x < dimensionSizes.nx; x++)
        {
          const size_t i = get1DIndex(z, y, x, dimensionSizes);
          // 3D and 2D summation
          const float sumRhos   = (simulationDimension == SD::k3D) ? (rhoX[i] + rhoY[i] + rhoZ[i])
                                                                   : (rhoX[i] + rhoY[i]);
          const float sumRhosDt = (2.0f * sumRhos + rho0) * dt;

          rhoX[i] = pmlX[x] * ((pmlX[x] * rhoX[i]) - sumRhosDt * duxdx[i]);
          rhoY[i] = pmlY[y] * ((pmlY[y] * rhoY[i]) - sumRhosDt * duydy[i]);
          if (simulationDimension == SD::k3D)
          {
            rhoZ[i] = pmlZ[z] * ((pmlZ[z] * rhoZ[i]) - sumRhosDt * duzdz[i]);
          }
        }// x
      }// y
    }// z
  }
  else
  { //---------------------------------------------- rho0 is matrix ------------------------------------------------//
    // rho0 is a matrix
    const float* rho0  = getRho0().getData();

    #pragma omp parallel for schedule(static) if (simulationDimension == SD::k3D)
    for (size_t z = 0; z < dimensionSizes.nz; z++)
    {
      #pragma omp parallel for schedule(static) if (simulationDimension == SD::k2D)
      for (size_t y = 0; y < dimensionSizes.ny; y++)
      {
        #pragma omp simd
        for (size_t x = 0; x < dimensionSizes.nx; x++)
        {
          const size_t i = get1DIndex(z, y, x, dimensionSizes);
          // 3D and 2D summation
          const float sumRhos   = (simulationDimension == SD::k3D) ? (rhoX[i] + rhoY[i] + rhoZ[i])
                                                                   : (rhoX[i] + rhoY[i]);

          const float sumRhosDt = (2.0f * sumRhos + rho0[i]) * dt;

          rhoX[i] = pmlX[x] * ((pmlX[x] * rhoX[i]) - sumRhosDt * duxdx[i]);
          rhoY[i] = pmlY[y] * ((pmlY[y] * rhoY[i]) - sumRhosDt * duydy[i]);
          if (simulationDimension == SD::k3D)
          {
            rhoZ[i] = pmlZ[z] * ((pmlZ[z] * rhoZ[i]) - sumRhosDt * duzdz[i]);
          }
        } // x
      }// y
    }// z
  } // end rho is matrix
}// end of computeDensityNonliner
//----------------------------------------------------------------------------------------------------------------------

/**
 * Calculate new values of acoustic density for linear case (rhoX, rhoy and rhoZ).
 *
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::computeDensityLinear()
{
  const DimensionSizes& dimensionSizes = mParameters.getFullDimensionSizes();
  const float dt = mParameters.getDt();

  float* rhoX  = getRhoX().getData();
  float* rhoY  = getRhoY().getData();
  float* rhoZ  = (simulationDimension == SD::k3D) ? getRhoZ().getData() : nullptr;

  const float* pmlX  = getPmlX().getData();
  const float* pmlY  = getPmlY().getData();
  const float* pmlZ  = (simulationDimension == SD::k3D) ? getPmlZ().getData() : nullptr;

  const float* duxdx = getDuxdx().getData();
  const float* duydy = getDuydy().getData();
  const float* duzdz = (simulationDimension == SD::k3D) ? getDuzdz().getData() : nullptr;

  //----------------------------------------------- rho0 is scalar -------------------------------------------------//
  if (mParameters.getRho0ScalarFlag())
  { // rho0 is a scalar
    const float dtRho0 = mParameters.getRho0Scalar() * dt;

    #pragma omp parallel for schedule(static) if (simulationDimension == SD::k3D)
    for (size_t z = 0; z < dimensionSizes.nz; z++)
    {
      #pragma omp parallel for schedule(static) if (simulationDimension == SD::k2D)
      for (size_t y = 0; y < dimensionSizes.ny; y++)
      {
        #pragma omp simd
        for (size_t x = 0; x < dimensionSizes.nx; x++)
        {
          const size_t i = get1DIndex(z, y, x, dimensionSizes);

          rhoX[i] = pmlX[x] * (((pmlX[x] * rhoX[i]) - (dtRho0 * duxdx[i])));
          rhoY[i] = pmlY[y] * (((pmlY[y] * rhoY[i]) - (dtRho0 * duydy[i])));
          if (simulationDimension == SD::k3D)
          {
            rhoZ[i] = pmlZ[z] * (((pmlZ[z] * rhoZ[i]) - (dtRho0 * duzdz[i])));
          }
        } // x
      }// y
    }// z
  }
  else
  { //---------------------------------------------- rho0 is matrix ------------------------------------------------//
    // rho0 is a matrix
    const float* rho0  = getRho0().getData();

    #pragma omp parallel for schedule(static) if (simulationDimension == SD::k3D)
    for (size_t z = 0; z < dimensionSizes.nz; z++)
    {
      #pragma omp parallel for schedule(static) if (simulationDimension == SD::k2D)
      for (size_t y = 0; y < dimensionSizes.ny; y++)
      {
        #pragma omp simd
        for (size_t x = 0; x < dimensionSizes.nx; x++)
        {
          const size_t i = get1DIndex(z, y, x, dimensionSizes);
          const float dtRho0 = dt * rho0[i];

          rhoX[i] = pmlX[x] * (((pmlX[x] * rhoX[i]) - (dtRho0 * duxdx[i])));
          rhoY[i] = pmlY[y] * (((pmlY[y] * rhoY[i]) - (dtRho0 * duydy[i])));
          if (simulationDimension == SD::k3D)
          {
            rhoZ[i] = pmlZ[z] * (((pmlZ[z] * rhoZ[i]) - (dtRho0 * duzdz[i])));
          }
        } // x
      }// y
    }// z
  } // end rho is a matrix
}// end of computeDensityLinear
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute acoustic pressure for non-linear case.
 *
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::computePressureNonlinear()
{
  if (mParameters.getAbsorbingFlag())
  { //----------------------------------------------- absorbing case--------------------------------------------------//
    RealMatrix& densitySum         = getTemp1RealND();
    RealMatrix& nonlinearTerm      = getTemp2RealND();
    RealMatrix& velocitGradientSum = getTemp3RealND();

    // reusing of the temp variables
    RealMatrix& absorbTauTerm = velocitGradientSum;
    RealMatrix& absorbEtaTerm = densitySum;

    // different templated variants of computePressureTermsNonlinear
    if ( mParameters.getBOnAScalarFlag())
    {
      if (mParameters.getRho0ScalarFlag())
      {
        computePressureTermsNonlinear<simulationDimension, true, true>(densitySum, nonlinearTerm, velocitGradientSum);
      }
      else
      {
        computePressureTermsNonlinear<simulationDimension, true, false>(densitySum, nonlinearTerm, velocitGradientSum);
      }
    }
    else
    {
      if (mParameters.getRho0ScalarFlag())
      {
        computePressureTermsNonlinear<simulationDimension, false, true>(densitySum, nonlinearTerm, velocitGradientSum);
      }
      else
      {
        computePressureTermsNonlinear<simulationDimension, false, false>(densitySum, nonlinearTerm, velocitGradientSum);
      }
    }

    // ifftn( absorb_nabla1 * fftn (rho0 * (duxdx+duydy+duzdz))
    getTempFftwX().computeR2CFftND(velocitGradientSum);
    getTempFftwY().computeR2CFftND(densitySum);

    computeAbsorbtionTerm(getTempFftwX(), getTempFftwY());

    getTempFftwX().computeC2RFftND(absorbTauTerm);
    getTempFftwY().computeC2RFftND(absorbEtaTerm);

    // different templated variants of sumPressureTermsNonlinear
    if (mParameters.getC0ScalarFlag())
    {
      if (mParameters.getAlphaCoeffScalarFlag())
      {
        sumPressureTermsNonlinear<true, true>(absorbTauTerm, absorbEtaTerm, nonlinearTerm);
      }
      else
      {
        sumPressureTermsNonlinear<true, false>(absorbTauTerm, absorbEtaTerm, nonlinearTerm);
      }
    }
    else
    {
      sumPressureTermsNonlinear<false, false>(absorbTauTerm, absorbEtaTerm, nonlinearTerm);
    }
  }
  else
  { //------------------------------------------------ lossless case--------------------------------------------------//
    if (mParameters.getC0ScalarFlag())
    {
      if (mParameters.getBOnAScalarFlag())
      {
        if (mParameters.getRho0ScalarFlag())
        {
          sumPressureTermsNonlinearLossless<simulationDimension, true, true, true>();
        }
        else
        {
          sumPressureTermsNonlinearLossless<simulationDimension, true, true, false>();
        }
      }
      else
      {
        if (mParameters.getRho0ScalarFlag())
        {
          sumPressureTermsNonlinearLossless<simulationDimension, true, false, true>();
        }
        else
        {
          sumPressureTermsNonlinearLossless<simulationDimension, true, false, false>();
        }
      }
    }
    else
    {
      if (mParameters.getBOnAScalarFlag())
      {
        if (mParameters.getRho0ScalarFlag())
        {
          sumPressureTermsNonlinearLossless<simulationDimension, false, true, true>();
        }
        else
        {
          sumPressureTermsNonlinearLossless<simulationDimension, false, true, false>();
        }
      }
      else
      {
        if (mParameters.getRho0ScalarFlag())
        {
          sumPressureTermsNonlinearLossless<simulationDimension, false, false, true>();
        }
        else
        {
          sumPressureTermsNonlinearLossless<simulationDimension, false, false, false>();
        }
      }
    }
  }
}// end of computePressureNonlinear
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute new p for linear case.
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::computePressureLinear()
{
  // rhox + rhoy + rhoz
  if (mParameters.getAbsorbingFlag())
  { // absorbing case

    RealMatrix& densitySum           = getTemp1RealND();
    RealMatrix& velocityGradientTerm = getTemp2RealND();

    RealMatrix& absorbTauTerm        = getTemp2RealND();
    RealMatrix& absorbEtaTerm        = getTemp3RealND();

    computePressureTermsLinear<simulationDimension>(densitySum, velocityGradientTerm);

    // ifftn ( absorb_nabla1 * fftn (rho0 * (duxdx+duydy+duzdz))

    getTempFftwX().computeR2CFftND(velocityGradientTerm);
    getTempFftwY().computeR2CFftND(densitySum);

    computeAbsorbtionTerm(getTempFftwX(), getTempFftwY());

    getTempFftwX().computeC2RFftND(absorbTauTerm);
    getTempFftwY().computeC2RFftND(absorbEtaTerm);

    if (mParameters.getC0ScalarFlag())
    {
      if (mParameters.getAlphaCoeffScalarFlag())
      {
        sumPressureTermsLinear<true, true>(absorbTauTerm, absorbEtaTerm, densitySum);
      }
      else
      {
        sumPressureTermsLinear<true, false>(absorbTauTerm, absorbEtaTerm, densitySum);
      }
    }
    else
    {
      sumPressureTermsLinear<false, false>(absorbTauTerm, absorbEtaTerm, densitySum);
    }
  }
  else
  {
    // lossless case
    sumPressureTermsLinearLossless<simulationDimension>();
  }
}// end of computePressureLinear
//----------------------------------------------------------------------------------------------------------------------

/**
 * Add u source to the particle velocity.
 */
void KSpaceFirstOrderSolver::addVelocitySource()
{
  const size_t timeIndex = mParameters.getTimeIndex();

  if (mParameters.getVelocityXSourceFlag() > timeIndex)
  {
    computeVelocitySourceTerm(getUxSgx(), GetVelocityXSourceInput(), getVelocitySourceIndex());
  }

  if (mParameters.getVelocityYSourceFlag() > timeIndex)
  {
    computeVelocitySourceTerm(getUySgy(), GetVelocityYSourceInput(), getVelocitySourceIndex());
  }

  if (mParameters.isSimulation3D())
  {
    if (mParameters.getVelocityZSourceFlag() > timeIndex)
    {
      computeVelocitySourceTerm(getUzSgz(), getVelocityZSourceInput(), getVelocitySourceIndex());
    }
  }
}// end of addVelocitySource
//----------------------------------------------------------------------------------------------------------------------

/**
 * Add in velocity source terms.
 */
void KSpaceFirstOrderSolver::computeVelocitySourceTerm(RealMatrix&        velocityMatrix,
                                                       const RealMatrix&  velocitySourceInput,
                                                       const IndexMatrix& velocitySourceIndex)
{
  // time is solved one level up.
  float*        pVelocityMatrix = velocityMatrix.getData();

  const float*  sourceInput = velocitySourceInput.getData();
  const size_t* sourceIndex = velocitySourceIndex.getData();

  const size_t timeIndex  = mParameters.getTimeIndex();

  const bool   isManyFlag = (mParameters.getVelocitySourceMany() != 0);
  const size_t sourceSize = velocitySourceIndex.size();
  const size_t index2D    = (isManyFlag) ? timeIndex * sourceSize : timeIndex;

  switch (mParameters.getVelocitySourceMode())
  {
    case Parameters::SourceMode::kDirichlet:
    {
      #pragma omp parallel for if (sourceSize > 16384)
      for (size_t i = 0; i < sourceSize; i++)
      {
        const size_t signalIndex = (isManyFlag) ? index2D + i : index2D;
        pVelocityMatrix[sourceIndex[i]] = sourceInput[signalIndex];
      }
      break;
    }

    case Parameters::SourceMode::kAdditiveNoCorrection:
    {
      #pragma omp parallel for if (sourceSize > 16384)
      for (size_t i = 0; i < sourceSize; i++)
      {
        const size_t signalIndex = (isManyFlag) ? index2D + i : index2D;
        pVelocityMatrix[sourceIndex[i]] += sourceInput[signalIndex];
      }
      break;
    }

    case Parameters::SourceMode::kAdditive:
    {
      // temp matrix for additive source
      RealMatrix&        scaledSource = getTemp1RealND();
      FftwComplexMatrix& fftMatrix    = getTempFftwX();

      float*        pScaledSource = getTemp1RealND().getData();
      float*        pSourceKappa  = getSourceKappa().getData();
      FloatComplex* pFftMatrix    = getTempFftwX().getComplexData();

      const size_t  nElementsFull    = mParameters.getFullDimensionSizes().nElements();
      const size_t  nElementsReduced = mParameters.getReducedDimensionSizes().nElements();
      const float   divider          = 1.0f / static_cast<float>(nElementsFull);

      // clear scaledSource the matrix
      scaledSource.zeroMatrix();

      // source_mat(u_source_pos_index) = source.u(u_source_sig_index, t_index);
      #pragma omp parallel for simd
      for (size_t i = 0; i < sourceSize; i++)
      {
        const size_t signalIndex = (isManyFlag) ? index2D + i : index2D;
        pScaledSource[sourceIndex[i]] = sourceInput[signalIndex];
      }

      // source_mat = real(ifftn(source_kappa .* fftn(source_mat)));
      fftMatrix.computeR2CFftND(scaledSource);

      #pragma omp parallel for simd
      for (size_t i = 0; i < nElementsReduced; i++)
      {
        pFftMatrix[i] *= divider * pSourceKappa[i];
      }

      // source_mat = real(ifftn(source_kappa .* fftn(source_mat)));
      fftMatrix.computeC2RFftND(scaledSource);

      // add the source values to the existing field values
      #pragma omp parallel for simd
      for (size_t i = 0; i < nElementsFull; i++)
      {
        pVelocityMatrix[i] += pScaledSource[i];
      }

      break;
    }
    default:
    {
      break;
    }
  } // end of switch
}// end of computeVelocitySourceTerm
//----------------------------------------------------------------------------------------------------------------------

 /**
  * Add in pressure source.
  */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::addPressureSource()
{
  const size_t timeIndex = mParameters.getTimeIndex();

  if (mParameters.getPressureSourceFlag() > timeIndex)
  {
    float* rhox = getRhoX().getData();
    float* rhoy = getRhoY().getData();
    float* rhoz = (simulationDimension == SD::k3D) ? getRhoZ().getData() : nullptr;

    const float*  sourceInput = getPressureSourceInput().getData();
    const size_t* sourceIndex = getPressureSourceIndex().getData();

    const bool   isManyFlag  = (mParameters.getPressureSourceMany() != 0);
    const size_t sourceSize  = getPressureSourceIndex().size();
    const size_t index2D     = (isManyFlag) ? timeIndex * sourceSize : timeIndex;

    // different pressure sources
    switch (mParameters.getPressureSourceMode())
    {
      case Parameters::SourceMode::kDirichlet:
      {
        #pragma omp parallel for if (sourceSize > 16384)
        for (size_t i = 0; i < sourceSize; i++)
        {
          const size_t signalIndex = (isManyFlag) ? index2D + i : index2D;

          rhox[sourceIndex[i]] = sourceInput[signalIndex];
          rhoy[sourceIndex[i]] = sourceInput[signalIndex];
          if (simulationDimension == SD::k3D)
          {
            rhoz[sourceIndex[i]] = sourceInput[signalIndex];
          }
        }
        break;
      }

      case Parameters::SourceMode::kAdditiveNoCorrection:
      {
        #pragma omp parallel for if (sourceSize > 16384)
        for (size_t i = 0; i < sourceSize; i++)
        {
          const size_t signalIndex = (isManyFlag) ? index2D + i : index2D;

          rhox[sourceIndex[i]] += sourceInput[signalIndex];
          rhoy[sourceIndex[i]] += sourceInput[signalIndex];
          if (simulationDimension == SD::k3D)
          {
            rhoz[sourceIndex[i]] += sourceInput[signalIndex];
          }
        }
        break;
      }

      case Parameters::SourceMode::kAdditive:
      { // temp matrix for additive source
        RealMatrix&        scaledSource = getTemp1RealND();
        FftwComplexMatrix& fftMatrix    = getTempFftwX();

        float*        pScaledSource = getTemp1RealND().getData();
        float*        pSourceKappa  = getSourceKappa().getData();
        FloatComplex* pFftMatrix    = getTempFftwX().getComplexData();

        const size_t  nElementsFull    = mParameters.getFullDimensionSizes().nElements();
        const size_t  nElementsReduced = mParameters.getReducedDimensionSizes().nElements();
        const float   divider          = 1.0f / static_cast<float>(nElementsFull);

        // clear scaledSource the matrix
        scaledSource.zeroMatrix();

        // source_mat(p_source_pos_index) = source.p(p_source_sig_index, t_index);
        #pragma omp parallel for simd
        for (size_t i = 0; i < sourceSize; i++)
        {
          const size_t signalIndex = (isManyFlag) ? index2D + i : index2D;
          pScaledSource[sourceIndex[i]] = sourceInput[signalIndex];
        }

        // source_mat = real(ifftn(source_kappa .* fftn(source_mat)));
        fftMatrix.computeR2CFftND(scaledSource);

        #pragma omp parallel for simd
        for (size_t i = 0; i < nElementsReduced ; i++)
        {
          pFftMatrix[i] *= divider * pSourceKappa[i];
        }

        // source_mat = real(ifftn(source_kappa .* fftn(source_mat)));
        fftMatrix.computeC2RFftND(scaledSource);

        // add the source values to the existing field values
        #pragma omp parallel for simd
        for (size_t i = 0; i < nElementsFull; i++)
        {
          rhox[i] += pScaledSource[i];
          rhoy[i] += pScaledSource[i];
          if (simulationDimension == SD::k3D)
          {
            rhoz[i] += pScaledSource[i];
          }
        }
        break;
      }
      default:
      {
        break;
      }
    } // switch
  }// if do at all
}// end of addPressureSource
//----------------------------------------------------------------------------------------------------------------------

/**
 * Calculate p0 source when necessary.
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::addInitialPressureSource()
{
  getP().copyData(getInitialPressureSourceInput());

  const float* sourceInput = getInitialPressureSourceInput().getData();

  const bool   c0ScalarFlag = mParameters.getC0ScalarFlag();
  const float  c2Scalar     = (c0ScalarFlag) ? mParameters.getC2Scalar() : 0;
  const float* c2Matrix     = (c0ScalarFlag) ? nullptr : getC2().getData();

  float* rhoX = getRhoX().getData();
  float* rhoY = getRhoY().getData();
  float* rhoZ = (simulationDimension == SD::k3D) ? getRhoZ().getData() : nullptr;

  const size_t nElements = mParameters.getFullDimensionSizes().nElements();
  const float  dimScalingFactor = (simulationDimension == SD::k3D) ? 3.0f : 2.0f;

  #pragma omp parallel for simd schedule(static)
  for (size_t i = 0; i < nElements; i++)
  {
    const float tmp = sourceInput[i] / (dimScalingFactor * ((c0ScalarFlag) ? c2Scalar : c2Matrix[i]));

    rhoX[i] = tmp;
    rhoY[i] = tmp;
    if (simulationDimension == SD::k3D)
    {
      rhoZ[i] = tmp;
    }
  }

  //------------------------------------------------------------------------//
  //--  compute u(t = t1 + dt/2) based on the assumption u(dt/2) = -u(-dt/2) --//
  //--    which forces u(t = t1) = 0 --//
  //------------------------------------------------------------------------//
  computePressureGradient<simulationDimension>();

  if (mParameters.getRho0ScalarFlag())
  {
    if (mParameters.getNonUniformGridFlag())
    { // non uniform grid, homogeneous case
      computeInitialVelocityHomogeneousNonuniform<simulationDimension>();
    }
    else
    { //uniform grid, homogeneous
      computeInitialVelocityHomogeneousUniform<simulationDimension>();
    }
  }
  else
  { // heterogeneous, unifrom grid
    // divide the matrix by 2 and multiply with st./rho0_sg
    computeInitialVelocityHeterogeneous<simulationDimension>();
  }
}// end of addInitialPressureSource
//----------------------------------------------------------------------------------------------------------------------

/**
 * Add transducer data source to velocity x component.
 */
void KSpaceFirstOrderSolver::addTransducerSource()
{
  float* uxSgx = getUxSgx().getData();

  const size_t* velocitySourceIndex   = getVelocitySourceIndex().getData();
  const float*  transducerSourceInput = getTransducerSourceInput().getData();
  const size_t* delayMask             = getDelayMask().getData();

  const size_t timeIndex  = mParameters.getTimeIndex();
  const size_t sourceSize = getVelocitySourceIndex().size();

  #pragma omp parallel for schedule(static) if (sourceSize > 16384)
  for (size_t i = 0; i < sourceSize; i++)
  {
    uxSgx[velocitySourceIndex[i]] += transducerSourceInput[delayMask[i] + timeIndex];
  }
}// end of addTransducerSource
//----------------------------------------------------------------------------------------------------------------------

/**
 * Generate kappa matrix for lossless medium.
 * For 2D simulation, the zPart == 0.
 */
void KSpaceFirstOrderSolver::generateKappa()
{
  const float dx2Rec = 1.0f / (mParameters.getDx() * mParameters.getDx());
  const float dy2Rec = 1.0f / (mParameters.getDy() * mParameters.getDy());
  // For 2D simulation set dz to 0
  const float dz2Rec = (mParameters.isSimulation3D()) ? 1.0f / (mParameters.getDz() * mParameters.getDz()) : 0.0f;

  const float cRefDtPi = mParameters.getCRef() * mParameters.getDt() * static_cast<float>(M_PI);

  const float nxRec = 1.0f / static_cast<float>(mParameters.getFullDimensionSizes().nx);
  const float nyRec = 1.0f / static_cast<float>(mParameters.getFullDimensionSizes().ny);
  // For 2D simulation, nzRec remains 1
  const float nzRec = 1.0f / static_cast<float>(mParameters.getFullDimensionSizes().nz);

  const DimensionSizes& reducedDimensionSizes = mParameters.getReducedDimensionSizes();

  float* kappa = getKappa().getData();

  #pragma omp parallel for schedule(static) if (mParameters.isSimulation3D())
  for (size_t z = 0; z < reducedDimensionSizes.nz; z++)
  {
    const float zf    = static_cast<float>(z);
          float zPart = 0.5f - fabs(0.5f - zf * nzRec);
                zPart = (zPart * zPart) * dz2Rec;

    #pragma omp parallel for schedule(static) if (mParameters.isSimulation2D())
    for (size_t y = 0; y < reducedDimensionSizes.ny; y++)
    {
      const float yf    = static_cast<float>(y);
            float yPart = 0.5f - fabs(0.5f - yf * nyRec);
                  yPart = (yPart * yPart) * dy2Rec;

      const float yzPart = zPart + yPart;
      for (size_t x = 0; x < reducedDimensionSizes.nx; x++)
      {
        const float xf = static_cast<float>(x);
              float xPart = 0.5f - fabs(0.5f - xf * nxRec);
                    xPart = (xPart * xPart) * dx2Rec;

              float k = cRefDtPi * sqrt(xPart + yzPart);

        const size_t i = get1DIndex(z, y, x, reducedDimensionSizes);

        // kappa element
        kappa[i] = (k == 0.0f) ? 1.0f : sin(k) / k;
      }//x
    }//y
  }// z
}// end of generateKappa
//----------------------------------------------------------------------------------------------------------------------

/**
 * Generate sourceKappa matrix for additive sources.
 * For 2D simulation, the zPart == 0.
 */
void KSpaceFirstOrderSolver::generateSourceKappa()
{
  const float dx2Rec = 1.0f / (mParameters.getDx() * mParameters.getDx());
  const float dy2Rec = 1.0f / (mParameters.getDy() * mParameters.getDy());
  const float dz2Rec = (mParameters.isSimulation3D()) ? 1.0f / (mParameters.getDz() * mParameters.getDz()) : 0.0f;

  const float cRefDtPi = mParameters.getCRef() * mParameters.getDt() * static_cast<float>(M_PI);

  const float nxRec = 1.0f / static_cast<float>(mParameters.getFullDimensionSizes().nx);
  const float nyRec = 1.0f / static_cast<float>(mParameters.getFullDimensionSizes().ny);
  const float nzRec = 1.0f / static_cast<float>(mParameters.getFullDimensionSizes().nz);

  const DimensionSizes& reducedDimensionSizes = mParameters.getReducedDimensionSizes();

  float* sourceKappa = getSourceKappa().getData();

  #pragma omp parallel for schedule(static) if (mParameters.isSimulation3D())
  for (size_t z = 0; z < reducedDimensionSizes.nz; z++)
  {
    const float zf    = static_cast<float>(z);
          float zPart = 0.5f - fabs(0.5f - zf * nzRec);
                zPart = (zPart * zPart) * dz2Rec;

    #pragma omp parallel for schedule(static) if (mParameters.isSimulation2D())
    for (size_t y = 0; y < reducedDimensionSizes.ny; y++)
    {
      const float yf    = static_cast<float>(y);
            float yPart = 0.5f - fabs(0.5f - yf * nyRec);
                  yPart = (yPart * yPart) * dy2Rec;

      const float yzPart = zPart + yPart;
      for (size_t x = 0; x < reducedDimensionSizes.nx; x++)
      {
        const float xf = static_cast<float>(x);
              float xPart = 0.5f - fabs(0.5f - xf * nxRec);
                    xPart = (xPart * xPart) * dx2Rec;

              float k = cRefDtPi * sqrt(xPart + yzPart);

        const size_t i = get1DIndex(z, y, x, reducedDimensionSizes);

        // sourceKappa element
        sourceKappa[i] = cos(k);
      }//x
    }//y
  }// z
}// end of generateSourceKappa
//----------------------------------------------------------------------------------------------------------------------

/**
 * Generate kappa matrix, absorbNabla1, absorbNabla2 for absorbing medium.
 * For the 2D simulation the zPart == 0
 */
void KSpaceFirstOrderSolver::generateKappaAndNablas()
{
  const float dxSqRec    = 1.0f / (mParameters.getDx() * mParameters.getDx());
  const float dySqRec    = 1.0f / (mParameters.getDy() * mParameters.getDy());
  const float dzSqRec    = (mParameters.isSimulation3D()) ? 1.0f / (mParameters.getDz() * mParameters.getDz()) : 0.0f;

  const float cRefDt2    = mParameters.getCRef() * mParameters.getDt() * 0.5f;
  const float pi2        = static_cast<float>(M_PI) * 2.0f;

  const size_t nx        = mParameters.getFullDimensionSizes().nx;
  const size_t ny        = mParameters.getFullDimensionSizes().ny;
  const size_t nz        = mParameters.getFullDimensionSizes().nz;

  const float nxRec      = 1.0f / static_cast<float>(nx);
  const float nyRec      = 1.0f / static_cast<float>(ny);
  const float nzRec      = 1.0f / static_cast<float>(nz);

  const DimensionSizes& reducedDimensionSizes = mParameters.getReducedDimensionSizes();

  float* kappa           = getKappa().getData();
  float* absorbNabla1    = getAbsorbNabla1().getData();
  float* absorbNabla2    = getAbsorbNabla2().getData();
  const float alphaPower = mParameters.getAlphaPower();

  #pragma omp parallel for schedule(static) if (mParameters.isSimulation3D())
  for (size_t z = 0; z < reducedDimensionSizes.nz; z++)
  {
    const float zf    = static_cast<float>(z);
          float zPart = 0.5f - fabs(0.5f - zf * nzRec);
                zPart = (zPart * zPart) * dzSqRec;

    #pragma omp parallel for schedule(static) if (mParameters.isSimulation2D())
    for (size_t y = 0; y < reducedDimensionSizes.ny; y++)
    {
      const float yf    = static_cast<float>(y);
            float yPart = 0.5f - fabs(0.5f - yf * nyRec);
                  yPart = (yPart * yPart) * dySqRec;

      const float yzPart = zPart + yPart;

      for (size_t x = 0; x < reducedDimensionSizes.nx; x++)
      {
        const float xf    = static_cast<float>(x);
              float xPart = 0.5f - fabs(0.5f - xf * nxRec);
                    xPart = (xPart * xPart) * dxSqRec;

              float k     = pi2 * sqrt(xPart + yzPart);
              float cRefK = cRefDt2 * k;

        const size_t i = get1DIndex(z, y, x, reducedDimensionSizes);

        kappa[i]          = (cRefK == 0.0f) ? 1.0f : sin(cRefK) / cRefK;

        absorbNabla1[i] = pow(k, alphaPower - 2);
        absorbNabla2[i] = pow(k, alphaPower - 1);

        if (absorbNabla1[i] ==  std::numeric_limits<float>::infinity()) absorbNabla1[i] = 0.0f;
        if (absorbNabla2[i] ==  std::numeric_limits<float>::infinity()) absorbNabla2[i] = 0.0f;
      }//x
    }//y
  }// z

}// end of generateKappaAndNablas
//----------------------------------------------------------------------------------------------------------------------

/**
 * Generate absorbTau and absorbEta in for heterogenous medium.
 */
void KSpaceFirstOrderSolver::generateTauAndEta()
{
  if ((mParameters.getAlphaCoeffScalarFlag()) && (mParameters.getC0ScalarFlag()))
  { // scalar values
    const float alphaPower       = mParameters.getAlphaPower();
    const float tanPi2AlphaPower = tan(static_cast<float> (M_PI_2) * alphaPower);
    const float alphaNeperCoeff  = (100.0f * pow(1.0e-6f / (2.0f * static_cast<float>(M_PI)), alphaPower)) /
                                   (20.0f * static_cast<float>(M_LOG10E));

    const float alphaCoeff2      = 2.0f * mParameters.getAlphaCoeffScalar() * alphaNeperCoeff;

    mParameters.setAbsorbTauScalar((-alphaCoeff2) * pow(mParameters.getC0Scalar(), alphaPower - 1));
    mParameters.setAbsorbEtaScalar(  alphaCoeff2  * pow(mParameters.getC0Scalar(), alphaPower) * tanPi2AlphaPower);
  }
  else
  { // matrix

    const DimensionSizes& dimensionSizes = mParameters.getFullDimensionSizes();

    float* absorbTau = getAbsorbTau().getData();
    float* absorbEta = getAbsorbEta().getData();

    const bool   alphaCoeffScalarFlag = mParameters.getAlphaCoeffScalarFlag();
    const float  alphaCoeffScalar     = (alphaCoeffScalarFlag) ? mParameters.getAlphaCoeffScalar() : 0;
    const float* alphaCoeffMatrix     = (alphaCoeffScalarFlag) ? nullptr : getTemp1RealND().getData();


    const bool   c0ScalarFlag = mParameters.getC0ScalarFlag();
    const float  c0Scalar     = (c0ScalarFlag) ? mParameters.getC0Scalar() : 0;
    // here c2 still holds just c0!
    const float* cOMatrix     = (c0ScalarFlag) ? nullptr : getC2().getData();


    const float alphaPower       = mParameters.getAlphaPower();
    const float tanPi2AlphaPower = tan(static_cast<float>(M_PI_2) * alphaPower);

    //alpha = 100*alpha.*(1e-6/(2*pi)).^y./
    //                  (20*log10(exp(1)));
    const float alphaNeperCoeff = (100.0f * pow(1.0e-6f / (2.0f * static_cast<float>(M_PI)), alphaPower)) /
                                  (20.0f * static_cast<float>(M_LOG10E));


    #pragma omp parallel for schedule(static) if (mParameters.isSimulation3D())
    for (size_t z = 0; z < dimensionSizes.nz; z++)
    {
      #pragma omp parallel for schedule(static) if (mParameters.isSimulation2D())
      for (size_t y = 0; y < dimensionSizes.ny; y++)
      {
        for (size_t x = 0; x < dimensionSizes.nx; x++)
        {
          const size_t i = get1DIndex(z, y, x, dimensionSizes);

          const float alphaCoeff2 = 2.0f * alphaNeperCoeff *
                                    ((alphaCoeffScalarFlag) ? alphaCoeffScalar : alphaCoeffMatrix[i]);

          absorbTau[i] = (-alphaCoeff2) * pow(((c0ScalarFlag) ? c0Scalar : cOMatrix[i]), alphaPower - 1);
          absorbEta[i] =   alphaCoeff2  * pow(((c0ScalarFlag) ? c0Scalar : cOMatrix[i]),
                                                alphaPower) * tanPi2AlphaPower;

        }//x
      }//y
    }// z
  } // matrix
}// end of generateTauAndEta
//----------------------------------------------------------------------------------------------------------------------

/**
 * Prepare dt./ rho0  for non-uniform grid.
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::generateInitialDenisty()
{
  float* dtRho0Sgx   = getDtRho0Sgx().getData();
  float* dtRho0Sgy   = getDtRho0Sgy().getData();
  float* dtRho0Sgz   = (simulationDimension == SD::k3D) ? getDtRho0Sgz().getData() : nullptr;

  const float dt = mParameters.getDt();

  const float* duxdxnSgx = getDxudxnSgx().getData();
  const float* duydynSgy = getDyudynSgy().getData();
  const float* duzdznSgz = (simulationDimension == SD::k3D) ? getDzudznSgz().getData() : nullptr;

  const DimensionSizes& dimensionSizes = mParameters.getFullDimensionSizes();

  #pragma omp parallel for schedule(static) if (simulationDimension == SD::k3D)
  for (size_t z = 0; z < dimensionSizes.nz; z++)
  {
    #pragma omp parallel for schedule(static) if (simulationDimension == SD::k2D)
    for (size_t y = 0; y < dimensionSizes.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < dimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, dimensionSizes);

        dtRho0Sgx[i] = (dt * duxdxnSgx[x]) / dtRho0Sgx[i];
        dtRho0Sgy[i] = (dt * duydynSgy[y]) / dtRho0Sgy[i];
        if (simulationDimension == SD::k3D)
        {
          dtRho0Sgz[i] = (dt * duzdznSgz[z]) / dtRho0Sgz[i];
        }
      } // x
    } // y
  } // z

}// end of generateInitialDenisty
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute c^2.
 */
void KSpaceFirstOrderSolver::computeC2()
{
  if (!mParameters.getC0ScalarFlag())
  {
    float* c2 =  getC2().getData();
    const size_t size = getC2().size();

    #pragma omp parallel for simd schedule(static) aligned(c2)
    for (size_t i=0; i < size; i++)
    {
      c2[i] = c2[i] * c2[i];
    }
  }// matrix
}// computeC2
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute acoustic velocity for initial pressure problem.
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::computeInitialVelocityHeterogeneous()
{
  getTempFftwX().computeC2RFftND(getUxSgx());
  getTempFftwY().computeC2RFftND(getUySgy());
  if (simulationDimension == SD::k3D)
  {
    getTempFftwZ().computeC2RFftND(getUzSgz());
  }

  const size_t nElements = mParameters.getFullDimensionSizes().nElements();
  const float  divider   = 1.0f / (2.0f * static_cast<float>(nElements));

  float* uxSgx = getUxSgx().getData();
  float* uySgy = getUySgy().getData();
  float* uzSgz = (simulationDimension == SD::k3D) ? getUzSgz().getData() : nullptr;

  const float* dtRho0Sgx = getDtRho0Sgx().getData();
  const float* dtRho0Sgy = getDtRho0Sgy().getData();
  const float* dtRho0Sgz = (simulationDimension == SD::k3D) ? getDtRho0Sgz().getData() : nullptr;


  #pragma omp parallel for simd schedule(static) \
          aligned(uxSgx, uySgy, uzSgz, dtRho0Sgx, dtRho0Sgy, dtRho0Sgz)
  for (size_t i = 0; i < nElements; i++)
  {
    uxSgx[i] *= dtRho0Sgx[i] * divider;
    uySgy[i] *= dtRho0Sgy[i] * divider;
    if (simulationDimension == SD::k3D)
    {
      uzSgz[i] *= dtRho0Sgz[i] * divider;
    }
  }
}// end of computeInitialVelocityHeterogeneous
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute velocity for the initial pressure problem, homogeneous medium, uniform grid.
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::computeInitialVelocityHomogeneousUniform()
{
  getTempFftwX().computeC2RFftND(getUxSgx());
  getTempFftwY().computeC2RFftND(getUySgy());
  if (simulationDimension == SD::k3D)
  {
    getTempFftwZ().computeC2RFftND(getUzSgz());
  }

  const size_t nElements = mParameters.getFullDimensionSizes().nElements();
  const float dividerX = 1.0f / (2.0f * static_cast<float>(nElements)) * mParameters.getDtRho0SgxScalar();
  const float dividerY = 1.0f / (2.0f * static_cast<float>(nElements)) * mParameters.getDtRho0SgyScalar();
  const float dividerZ = ((simulationDimension == SD::k3D))
                              ? 1.0f / (2.0f * static_cast<float>(nElements)) * mParameters.getDtRho0SgzScalar()
                              : 1.0f;

  float* uxSgx = getUxSgx().getData();
  float* uySgy = getUySgy().getData();
  float* uzSgz = (simulationDimension == SD::k3D) ? getUzSgz().getData() : nullptr;

  #pragma omp parallel for simd schedule(static) aligned(uxSgx, uySgy, uzSgz)
  for (size_t i = 0; i < nElements; i++)
  {
    uxSgx[i] *= dividerX;
    uySgy[i] *= dividerY;
    if (simulationDimension == SD::k3D)
    {
      uzSgz[i] *= dividerZ;
    }
  }
}// end of computeInitialVelocityHomogeneousUniform
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute acoustic velocity for initial pressure problem, homogenous medium, nonuniform grid.
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::computeInitialVelocityHomogeneousNonuniform()
{
  getTempFftwX().computeC2RFftND(getUxSgx());
  getTempFftwY().computeC2RFftND(getUySgy());
  if (simulationDimension == SD::k3D)
  {
    getTempFftwZ().computeC2RFftND(getUzSgz());
  }

  const DimensionSizes& dimensionSizes = mParameters.getFullDimensionSizes();
  const size_t nElements               = dimensionSizes.nElements();

  const float dividerX = 1.0f / (2.0f * static_cast<float>(nElements)) * mParameters.getDtRho0SgxScalar();
  const float dividerY = 1.0f / (2.0f * static_cast<float>(nElements)) * mParameters.getDtRho0SgyScalar();
  const float dividerZ = (simulationDimension == SD::k3D)
                             ? 1.0f / (2.0f * static_cast<float>(nElements)) * mParameters.getDtRho0SgzScalar()
                             : 1.0f;

  const float* dxudxnSgx = getDxudxnSgx().getData();
  const float* dyudynSgy = getDyudynSgy().getData();
  const float* dzudznSgz = (simulationDimension == SD::k3D) ? getDzudznSgz().getData() : nullptr;

  float* uxSgx = getUxSgx().getData();
  float* uySgy = getUySgy().getData();
  float* uzSgz = (simulationDimension == SD::k3D) ? getUzSgz().getData() : nullptr;


  #pragma omp parallel for schedule(static) if (simulationDimension == SD::k3D)
  for (size_t z = 0; z < dimensionSizes.nz; z++)
  {
    #pragma omp parallel for schedule(static) if (simulationDimension == SD::k2D)
    for (size_t y = 0; y < dimensionSizes.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < dimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, dimensionSizes);
        uxSgx[i] *= dividerX * dxudxnSgx[x];
        uySgy[i] *= dividerY * dyudynSgy[y];
        if ((simulationDimension == SD::k3D))
        {
          uzSgz[i] *= dividerZ * dzudznSgz[z];
        }
      } // x
    } // y
  } // z
}// end of computeInitialVelocityHomogeneousNonuniform
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute acoustic velocity for heterogeneous medium and a uniform grid, x direction.
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::computeVelocityHeterogeneous()
{
  const DimensionSizes& dimensionSizes = mParameters.getFullDimensionSizes();
  const size_t nElements   = dimensionSizes.nElements();
  const float  divider     = 1.0f / static_cast<float>(nElements);

  const float* ifftX = getTemp1RealND().getData();
  const float* ifftY = getTemp2RealND().getData();
  const float* ifftZ = (simulationDimension == SD::k3D) ? getTemp3RealND().getData() : nullptr;

  const float* dtRho0Sgx = getDtRho0Sgx().getData();
  const float* dtRho0Sgy = getDtRho0Sgy().getData();
  const float* dtRho0Sgz = (simulationDimension == SD::k3D) ? getDtRho0Sgz().getData() : nullptr;

  const float* pmlX = getPmlXSgx().getData();
  const float* pmlY = getPmlYSgy().getData();
  const float* pmlZ = (simulationDimension == SD::k3D) ? getPmlZSgz().getData() : nullptr;

  float* uxSgx = getUxSgx().getData();
  float* uySgy = getUySgy().getData();
  float* uzSgz = (simulationDimension == SD::k3D) ? getUzSgz().getData() : nullptr;

  // long loops are replicated for every dimension to save SIMD registers
  #pragma omp parallel for schedule(static) if (simulationDimension == SD::k3D)
  for (size_t z = 0; z < dimensionSizes.nz; z++)
  {
    #pragma omp parallel for schedule(static) if (simulationDimension == SD::k2D)
    for (size_t y = 0; y < dimensionSizes.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < dimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, dimensionSizes);

        uxSgx[i] = (uxSgx[i] * pmlX[x] - divider * ifftX[i] * dtRho0Sgx[i]) * pmlX[x];
      } // x
    } // y
  } // z

  #pragma omp parallel for schedule(static) if (simulationDimension == SD::k3D)
  for (size_t z = 0; z < dimensionSizes.nz; z++)
  {
    #pragma omp parallel for schedule(static) if (simulationDimension == SD::k2D)
    for (size_t y = 0; y < dimensionSizes.ny; y++)
    {
      const float ePmlY = pmlY[y];
      #pragma omp simd
      for (size_t x = 0; x < dimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, dimensionSizes);

        uySgy[i] = (uySgy[i] * ePmlY - divider * ifftY[i] * dtRho0Sgy[i]) * ePmlY;
      } // x
    } // y
  } // z

  if (simulationDimension == SD::k3D)
  {
    #pragma omp parallel for schedule(static)
    for (size_t z = 0; z < dimensionSizes.nz; z++)
    {
      const float ePmlZ = pmlZ[z];
      for (size_t y = 0; y < dimensionSizes.ny; y++)
      {
        #pragma omp simd
        for (size_t x = 0; x < dimensionSizes.nx; x++)
        {
          const size_t i = get1DIndex(z, y, x, dimensionSizes);

          uzSgz[i] = (uzSgz[i] * ePmlZ - divider * ifftZ[i] * dtRho0Sgz[i]) * ePmlZ;
        } // x
      } // y
    } // z
  } // k3D
}// end of computeVelocityHeterogeneous
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute acoustic velocity for homogeneous medium and a uniform grid.
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::computeVelocityHomogeneousUniform()
{
  const DimensionSizes& dimensionSizes = mParameters.getFullDimensionSizes();
  const size_t nElements = dimensionSizes.nElements();

  const float dividerX = mParameters.getDtRho0SgxScalar() / static_cast<float>(nElements);
  const float dividerY = mParameters.getDtRho0SgyScalar() / static_cast<float>(nElements);
  const float dividerZ = (simulationDimension == SD::k3D)
                             ? mParameters.getDtRho0SgzScalar() / static_cast<float>(nElements) : 1.0f;

  const float* ifftX = getTemp1RealND().getData();
  const float* ifftY = getTemp2RealND().getData();
  const float* ifftZ = (simulationDimension == SD::k3D) ? getTemp3RealND().getData() : nullptr;

  const float* pmlX = getPmlXSgx().getData();
  const float* pmlY = getPmlYSgy().getData();
  const float* pmlZ = (simulationDimension == SD::k3D) ? getPmlZSgz().getData() : nullptr;

  float* uxSgx = getUxSgx().getData();
  float* uySgy = getUySgy().getData();
  float* uzSgz = (simulationDimension == SD::k3D) ? getUzSgz().getData() : nullptr;

  // long loops are replicated for every dimension to save SIMD registers
  #pragma omp parallel for schedule(static) if (simulationDimension == SD::k3D)
  for (size_t z = 0; z < dimensionSizes.nz; z++)
  {
    #pragma omp parallel for schedule(static) if (simulationDimension == SD::k2D)
    for (size_t y = 0; y < dimensionSizes.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < dimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, dimensionSizes);

        uxSgx[i] = (uxSgx[i] * pmlX[x] - dividerX * ifftX[i]) * pmlX[x];
      } // x
    } // y
  } // z

  #pragma omp parallel for schedule(static) if (simulationDimension == SD::k3D)
  for (size_t z = 0; z < dimensionSizes.nz; z++)
  {
    #pragma omp parallel for schedule(static) if (simulationDimension == SD::k2D)
    for (size_t y = 0; y < dimensionSizes.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < dimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, dimensionSizes);

        uySgy[i] = (uySgy[i] * pmlY[y] - dividerY * ifftY[i]) * pmlY[y];
      } // x
    } // y
  } // z

  if (simulationDimension == SD::k3D)
  {
  #pragma omp parallel for schedule(static)
    for (size_t z = 0; z < dimensionSizes.nz; z++)
    {
      for (size_t y = 0; y < dimensionSizes.ny; y++)
      {
        #pragma omp simd
        for (size_t x = 0; x < dimensionSizes.nx; x++)
        {
          const size_t i = get1DIndex(z, y, x, dimensionSizes);

          uzSgz[i] = (uzSgz[i] * pmlZ[z] - dividerZ * ifftZ[i]) * pmlZ[z];
        } // x
      } // y
    } // z
  } // k3D
}// end of computeVelocityXHomogeneousUniform
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute acoustic velocity for homogenous medium and nonuniform grid, x direction.
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::computeVelocityHomogeneousNonuniform()
{
  const DimensionSizes& dimensionSizes = mParameters.getFullDimensionSizes();
  const size_t nElements = dimensionSizes.nElements();

  const float dividerX = mParameters.getDtRho0SgxScalar() / static_cast<float>(nElements);
  const float dividerY = mParameters.getDtRho0SgyScalar() / static_cast<float>(nElements);
  const float dividerZ = (simulationDimension == SD::k3D)
                             ? mParameters.getDtRho0SgzScalar() / static_cast<float>(nElements) : 1.0f;


  const float* ifftX = getTemp1RealND().getData();
  const float* ifftY = getTemp2RealND().getData();
  const float* ifftZ = (simulationDimension == SD::k3D) ? getTemp3RealND().getData() : nullptr;

  const float* dxudxnSgx = getDxudxnSgx().getData();
  const float* dyudynSgy = getDyudynSgy().getData();
  const float* dzudznSgz = (simulationDimension == SD::k3D) ? getDzudznSgz().getData() : nullptr;

  const float* pmlX = getPmlXSgx().getData();
  const float* pmlY = getPmlYSgy().getData();
  const float* pmlZ = (simulationDimension == SD::k3D) ? getPmlZSgz().getData() : nullptr;

  float* uxSgx = getUxSgx().getData();
  float* uySgy = getUySgy().getData();
  float* uzSgz = (simulationDimension == SD::k3D) ? getUzSgz().getData() : nullptr;

  // long loops are replicated for every dimension to save SIMD registers
  #pragma omp parallel for schedule(static) if (simulationDimension == SD::k3D)
  for (size_t z = 0; z < dimensionSizes.nz; z++)
  {
    #pragma omp parallel for schedule(static) if (simulationDimension == SD::k2D)
    for (size_t y = 0; y < dimensionSizes.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < dimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, dimensionSizes);

        uxSgx[i] = (uxSgx[i] * pmlX[x] - (dividerX * dxudxnSgx[x] * ifftX[i])) * pmlX[x];
      } // x
    } // y
  } // z

  #pragma omp parallel for schedule(static) if (simulationDimension == SD::k3D)
  for (size_t z = 0; z < dimensionSizes.nz; z++)
  {
    #pragma omp parallel for schedule(static) if (simulationDimension == SD::k2D)
    for (size_t y = 0; y < dimensionSizes.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < dimensionSizes.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, dimensionSizes);

        uySgy[i] = (uySgy[i] * pmlY[y] - (dividerY * dyudynSgy[y] * ifftY[i])) * pmlY[y];
      } // x
    } // y
  } // z

  if (simulationDimension == SD::k3D)
  {
    #pragma omp parallel for schedule(static)
    for (size_t z = 0; z < dimensionSizes.nz; z++)
    {
      for (size_t y = 0; y < dimensionSizes.ny; y++)
      {
        #pragma omp simd
        for (size_t x = 0; x < dimensionSizes.nx; x++)
        {
          const size_t i = get1DIndex(z, y, x, dimensionSizes);

          uzSgz[i] = (uzSgz[i] * pmlZ[z] - (dividerZ * dzudznSgz[z] * ifftZ[i])) * pmlZ[z];
        } // x
      } // y
    } // z
  } // k3D
}// end of computeVelocityHomogeneousNonuniform
//----------------------------------------------------------------------------------------------------------------------

/**
 *  Compute part of the new velocity term - gradient of pressure.
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::computePressureGradient()
{
  // Compute FFT of pressure
  getTempFftwX().computeR2CFftND(getP());

  FloatComplex* ifftX = getTempFftwX().getComplexData();
  FloatComplex* ifftY = getTempFftwY().getComplexData();
  FloatComplex* ifftZ = (simulationDimension == SD::k3D) ? getTempFftwZ().getComplexData() : nullptr;

  const FloatComplex* ddxKShiftPos = getDdxKShiftPos().getComplexData();
  const FloatComplex* ddyKShiftPos = getDdyKShiftPos().getComplexData();
  const FloatComplex* ddzKShiftPos = (simulationDimension == SD::k3D) ? getDdzKShiftPos().getComplexData() : nullptr;

  const float* kappa  = getKappa().getData();

  const DimensionSizes& reducedDimensionSizes= mParameters.getReducedDimensionSizes();

  #pragma omp parallel for schedule(static) if (simulationDimension == SD::k3D)
  for (size_t z = 0; z < reducedDimensionSizes.nz; z++)
  {
    #pragma omp parallel for schedule(static) if (simulationDimension == SD::k2D)
    for (size_t y = 0; y < reducedDimensionSizes.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < reducedDimensionSizes.nx;  x++)
      {
        const size_t i = get1DIndex(z, y, x, reducedDimensionSizes);

        const FloatComplex eKappa = ifftX[i] * kappa[i];

        ifftX[i] = eKappa * ddxKShiftPos[x];
        ifftY[i] = eKappa * ddyKShiftPos[y];
        if (simulationDimension == SD::k3D)
        {
          ifftZ[i] = eKappa * ddzKShiftPos[z];
        }
      } // x
    } // y
  } // z
}// end of computePressureGradient
//----------------------------------------------------------------------------------------------------------------------

/**
 * Calculate three temporary sums in the new pressure formula non-linear absorbing case.
 */
template<Parameters::SimulationDimension simulationDimension,
         bool bOnAScalarFlag,
         bool rho0ScalarFlag>
void KSpaceFirstOrderSolver::computePressureTermsNonlinear(RealMatrix& densitySum,
                                                           RealMatrix& nonlinearTerm,
                                                           RealMatrix& velocityGradientSum)
{
  const float* rhoX = getRhoX().getData();
  const float* rhoY = getRhoY().getData();
  const float* rhoZ = (simulationDimension == SD::k3D) ? getRhoZ().getData() : nullptr;

  const float* duxdx = getDuxdx().getData();
  const float* duydy = getDuydy().getData();
  const float* duzdz = (simulationDimension == SD::k3D) ? getDuzdz().getData() : nullptr;

  const float  bOnAScalar     = (bOnAScalarFlag) ? mParameters.getBOnAScalar() : 0;
  const float* bOnAMatrix     = (bOnAScalarFlag) ? nullptr : getBOnA().getData();

  const float  rho0Scalar     = (rho0ScalarFlag) ? mParameters.getRho0Scalar() : 0;
  const float* rho0Matrix     = (rho0ScalarFlag) ? nullptr : getRho0().getData();


  float* eDensitySum          = densitySum.getData();
  float* eNonlinearTerm       = nonlinearTerm.getData();
  float* eVelocityGradientSum = velocityGradientSum.getData();

  const size_t nElements = mParameters.getFullDimensionSizes().nElements();

  #pragma omp parallel for simd schedule(static) \
          aligned(eDensitySum, eNonlinearTerm, eVelocityGradientSum, \
                  rhoX, rhoY, rhoZ, bOnAMatrix, rho0Matrix, duxdx, duydy, duzdz)
  for (size_t i = 0; i < nElements ; i++)
  {
    const float rhoSum = (simulationDimension == SD::k3D) ? (rhoX[i]  + rhoY[i]  + rhoZ[i])  : (rhoX[i]  + rhoY[i]);
    const float duSum  = (simulationDimension == SD::k3D) ? (duxdx[i] + duydy[i] + duzdz[i]) : (duxdx[i] + duydy[i]);

    const float bOnA   = (bOnAScalarFlag) ? bOnAScalar : bOnAMatrix[i];
    const float rho0   = (rho0ScalarFlag) ? rho0Scalar : rho0Matrix[i];

    eDensitySum[i]          = rhoSum;
    eNonlinearTerm[i]       = (bOnA * rhoSum * rhoSum) / (2.0f * rho0) + rhoSum;
    eVelocityGradientSum[i] = rho0 * duSum;
  }
} // end of computePressureTermsNonlinear
//----------------------------------------------------------------------------------------------------------------------

/**
 * Calculate two temporary sums in the new pressure formula, linear absorbing case.
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::computePressureTermsLinear(RealMatrix& densitySum,
                                                        RealMatrix& velocityGradientSum)
{
  const size_t size = mParameters.getFullDimensionSizes().nElements();

  const float* rhoX = getRhoX().getData();
  const float* rhoY = getRhoY().getData();
  const float* rhoZ = (simulationDimension == SD::k3D) ? getRhoZ().getData() : nullptr;

  const float* duxdx = getDuxdx().getData();
  const float* duydy = getDuydy().getData();
  const float* duzdz = (simulationDimension == SD::k3D) ? getDuzdz().getData() : nullptr;

  float* pDensitySum          = densitySum.getData();
  float* pVelocityGradientSum = velocityGradientSum.getData();

  #pragma omp parallel for simd schedule(static) aligned (pDensitySum, rhoX, rhoY, rhoZ)
  for (size_t i = 0; i < size; i++)
  {
    pDensitySum[i] = (simulationDimension == SD::k3D) ? (rhoX[i] + rhoY[i] + rhoZ[i]) : (rhoX[i] + rhoY[i]);
  }

  if (mParameters.getRho0ScalarFlag())
  { // scalar
    const float eRho0 = mParameters.getRho0Scalar();
    #pragma omp parallel for simd schedule(static) aligned (pDensitySum, duxdx, duydy, duzdz)
    for (size_t i = 0; i < size; i++)
    {
      const float duSum = (simulationDimension == SD::k3D) ? (duxdx[i] + duydy[i] + duzdz[i]) : (duxdx[i] + duydy[i]);
      pVelocityGradientSum[i] = eRho0 * duSum;
    }
  }
  else
  { // matrix
    const float* rho0 = getRho0().getData();
    #pragma omp parallel for simd schedule(static) aligned (pDensitySum, rho0, duxdx, duydy, duzdz)
    for (size_t i = 0; i < size; i++)
    {
      const float duSum = (simulationDimension == SD::k3D) ? (duxdx[i] + duydy[i] + duzdz[i]) : (duxdx[i] + duydy[i]);
      pVelocityGradientSum[i] = rho0[i] * duSum;
    }
  }
}// end of computePressureTermsLinear
//----------------------------------------------------------------------------------------------------------------------


/**
 * Compute absorbing term with abosrbNabla1 and absorbNabla2.
 */
void KSpaceFirstOrderSolver::computeAbsorbtionTerm(FftwComplexMatrix& fftPart1,
                                                   FftwComplexMatrix& fftPart2)
{
  const size_t nElements    = mParameters.getReducedDimensionSizes().nElements();

  FloatComplex* pFftPart1 = fftPart1.getComplexData();
  FloatComplex* pFftPart2 = fftPart2.getComplexData();

  const float* absorbNabla1 = getAbsorbNabla1().getData();
  const float* absorbNabla2 = getAbsorbNabla2().getData();

  #pragma omp parallel for simd schedule(static) aligned(pFftPart1, pFftPart2, absorbNabla1, absorbNabla2)
  for (size_t i = 0; i < nElements; i++)
  {
    pFftPart1[i] *= absorbNabla1[i];
    pFftPart2[i] *= absorbNabla2[i];
  }
} // end of computeAbsorbtionTerm
//----------------------------------------------------------------------------------------------------------------------

/**
 * Sum sub-terms to calculate new pressure, after FFTs, non-linear case.
 */
template<bool c0ScalarFlag, bool areTauAndEtaScalars>
void KSpaceFirstOrderSolver::sumPressureTermsNonlinear(const RealMatrix& absorbTauTerm,
                                                       const RealMatrix& absorbEtaTerm,
                                                       const RealMatrix& nonlinearTerm)
{
  const float* pAbsorbTauTerm = absorbTauTerm.getData();
  const float* pAbsorbEtaTerm = absorbEtaTerm.getData();

  const size_t nElements = mParameters.getFullDimensionSizes().nElements();
  const float  divider = 1.0f / static_cast<float>(nElements);

  const float  c2Scalar     = (c0ScalarFlag) ? mParameters.getC2Scalar() : 0;
  const float* c2Matrix     = (c0ScalarFlag) ? nullptr : getC2().getData();

  const float  absorbTauScalar = (areTauAndEtaScalars) ? mParameters.getAbsorbTauScalar() : 0;
  const float* absorbTauMatrix = (areTauAndEtaScalars) ? nullptr : getAbsorbTau().getData();

  const float  absorbEtaScalar = (areTauAndEtaScalars) ? mParameters.getAbsorbEtaScalar() : 0;
  const float* absorbEtaMatrix = (areTauAndEtaScalars) ? nullptr : getAbsorbEta().getData();;

  const float* bOnA = nonlinearTerm.getData();
  float*       p    = getP().getData();

  #pragma omp parallel for simd schedule(static) \
          aligned(p, c2Matrix, pAbsorbTauTerm, absorbTauMatrix, pAbsorbEtaTerm, absorbEtaMatrix)
  for (size_t i = 0; i < nElements; i++)
  {
    const float c2        = (c0ScalarFlag) ?        c2Scalar        : c2Matrix[i];
    const float absorbTau = (areTauAndEtaScalars) ? absorbTauScalar : absorbTauMatrix[i];
    const float absorbEta = (areTauAndEtaScalars) ? absorbEtaScalar : absorbEtaMatrix[i];

    p[i] = c2 *(bOnA[i] + (divider * ((pAbsorbTauTerm[i] * absorbTau) - (pAbsorbEtaTerm[i] * absorbEta))));
  }
}// end of sumPressureTermsNonlinear
//----------------------------------------------------------------------------------------------------------------------

/**
 * Sum sub-terms to calculate new pressure, after FFTs, linear case.
 */
template<bool c0ScalarFlag, bool areTauAndEtaScalars>
void KSpaceFirstOrderSolver::sumPressureTermsLinear(const RealMatrix& absorbTauTerm,
                                                    const RealMatrix& absorbEtaTerm,
                                                    const RealMatrix& densitySum)
{
  const float* pAbsorbTauTerm = absorbTauTerm.getData();
  const float* pAbsorbEtaTerm = absorbEtaTerm.getData();

  const size_t nElements = mParameters.getFullDimensionSizes().nElements();
  const float  divider = 1.0f / static_cast<float>(nElements);

  const float  c2Scalar     = (c0ScalarFlag) ? mParameters.getC2Scalar() : 0;
  const float* c2Matrix     = (c0ScalarFlag) ? nullptr : getC2().getData();

  const float  absorbTauScalar     = (areTauAndEtaScalars) ? mParameters.getAbsorbTauScalar() : 0;
  const float* absorbTauMatrix     = (areTauAndEtaScalars) ? nullptr : getAbsorbTau().getData();

  const float  absorbEtaScalar     = (areTauAndEtaScalars) ? mParameters.getAbsorbEtaScalar() : 0;
  const float* absorbEtaMatrix     = (areTauAndEtaScalars) ? nullptr : getAbsorbEta().getData();;

  const float* pDenistySum = densitySum.getData();
        float* p           = getP().getData();

  #pragma omp parallel for simd schedule(static) \
          aligned (p, pDenistySum, c2Matrix, absorbTauMatrix, absorbEtaMatrix, pAbsorbTauTerm, pAbsorbEtaTerm)
  for (size_t i = 0; i < nElements; i++)
  {
    const float c2        = (c0ScalarFlag) ?        c2Scalar        : c2Matrix[i];
    const float absorbTau = (areTauAndEtaScalars) ? absorbTauScalar : absorbTauMatrix[i];
    const float absorbEta = (areTauAndEtaScalars) ? absorbEtaScalar : absorbEtaMatrix[i];

    p[i] = c2 * (pDenistySum[i] + (divider * ((pAbsorbTauTerm[i] * absorbTau) - (pAbsorbEtaTerm[i] * absorbEta))));
  }
}// end of sumPressureTermsLinear
//----------------------------------------------------------------------------------------------------------------------

/**
 * Sum sub-terms for new p, nonlinear lossless case.
 */
template<Parameters::SimulationDimension simulationDimension,
         bool c0ScalarFlag,
         bool nonlinearFlag,
         bool rho0ScalarFlag>
void KSpaceFirstOrderSolver::sumPressureTermsNonlinearLossless()
{
  const size_t nElements = mParameters.getFullDimensionSizes().nElements();

  float* p = getP().getData();

  const float* rhoX = getRhoX().getData();
  const float* rhoY = getRhoY().getData();
  const float* rhoZ = (simulationDimension == SD::k3D) ? getRhoZ().getData() : nullptr;

  const float  c2Scalar     = (c0ScalarFlag) ? mParameters.getC2Scalar() : 0;
  const float* c2Matrix     = (c0ScalarFlag) ? nullptr : getC2().getData();

  const float  bOnAScalar   = (nonlinearFlag) ? mParameters.getBOnAScalar(): 0;
  const float* bOnAMatrix   = (nonlinearFlag) ? nullptr : getBOnA().getData();

  const float  rho0Scalar   = (rho0ScalarFlag) ? mParameters.getRho0Scalar() : 0;
  const float* rho0Matrix   = (rho0ScalarFlag) ? nullptr : getRho0().getData();

  #pragma omp parallel for simd schedule (static)
  for (size_t i = 0; i < nElements; i++)
  {
    const float c2   = (c0ScalarFlag)   ? c2Scalar   : c2Matrix[i];
    const float bOnA = (nonlinearFlag)  ? bOnAScalar : bOnAMatrix[i];
    const float rho0 = (rho0ScalarFlag) ? rho0Scalar : rho0Matrix[i];

    const float sumDensity = (simulationDimension == SD::k3D) ? (rhoX[i] + rhoY[i] + rhoZ[i]) : (rhoX[i] + rhoY[i]);

    p[i] = c2 * (sumDensity + (bOnA * (sumDensity * sumDensity) / (2.0f * rho0)));
  }
}// end of sumPressureTermsNonlinearLossless
//----------------------------------------------------------------------------------------------------------------------

/**
 * Sum sub-terms for new pressure, linear lossless case.
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::sumPressureTermsLinearLossless()
{
  const float* rhoX = getRhoX().getData();
  const float* rhoY = getRhoY().getData();
  const float* rhoZ = (simulationDimension == SD::k3D) ? getRhoZ().getData() : nullptr;
        float* p    = getP().getData();

  const size_t nElements = mParameters.getFullDimensionSizes().nElements();

  if (mParameters.getC0ScalarFlag())
  {
    const float c2 = mParameters.getC2Scalar();

    #pragma omp parallel for simd schedule(static) aligned(p, rhoX, rhoY, rhoZ)
    for (size_t i = 0; i < nElements; i++)
    {
      const float sumRhos = (simulationDimension == SD::k3D) ? (rhoX[i] + rhoY[i] + rhoZ[i]) : (rhoX[i] + rhoY[i]);
      p[i] = c2 * sumRhos;
    }
  }
  else
  {
    const float* c2 = getC2().getData();

    #pragma omp parallel for simd schedule(static) aligned(p, c2, rhoX, rhoY, rhoZ)
    for (size_t i = 0; i < nElements; i++)
    {
      const float sumRhos = (simulationDimension == SD::k3D) ? (rhoX[i] + rhoY[i] + rhoZ[i]) : (rhoX[i] + rhoY[i]);
      p[i] = c2[i] * sumRhos;
    }
  }
}// end of sumPressureTermsLinearLossless()
//----------------------------------------------------------------------------------------------------------------------

/**
 * Calculated shifted velocities.
 */
template<Parameters::SimulationDimension simulationDimension>
void KSpaceFirstOrderSolver::computeShiftedVelocity()
{
  const FloatComplex* xShiftNegR  = getXShiftNegR().getComplexData();
  const FloatComplex* yShiftNegR  = getYShiftNegR().getComplexData();
  const FloatComplex* zShiftNegR  = (simulationDimension == SD::k3D) ? getZShiftNegR().getComplexData(): nullptr;

        FloatComplex* tempFftShift = getTempFftwShift().getComplexData();

  // sizes of frequency spaces
  DimensionSizes xShiftDims    = mParameters.getFullDimensionSizes();
                 xShiftDims.nx = xShiftDims.nx / 2 + 1;

  DimensionSizes yShiftDims    = mParameters.getFullDimensionSizes();
                 yShiftDims.ny = yShiftDims.ny / 2 + 1;

  // This remains 1 for 2D simulation
  DimensionSizes zShiftDims    = mParameters.getFullDimensionSizes();
                 zShiftDims.nz = zShiftDims.nz / 2 + 1;

  // normalization constants for FFTs
  const float dividerX = 1.0f / static_cast<float>(mParameters.getFullDimensionSizes().nx);
  const float dividerY = 1.0f / static_cast<float>(mParameters.getFullDimensionSizes().ny);
  // This remains 1 for 2D simulation
  const float dividerZ = 1.0f / static_cast<float>(mParameters.getFullDimensionSizes().nz);

  //-------------------------------------------------- ux_shifted ----------------------------------------------------//
  getTempFftwShift().computeR2CFft1DX(getUxSgx());

  #pragma omp parallel for schedule(static) if (simulationDimension == SD::k3D)
  for (size_t z = 0; z < xShiftDims.nz; z++)
  {
    #pragma omp parallel for schedule(static) if (simulationDimension == SD::k2D)
    for (size_t y = 0; y < xShiftDims.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < xShiftDims.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, xShiftDims);

        tempFftShift[i] = tempFftShift[i] * xShiftNegR[x] * dividerX;
      } // x
    } // y
  }//z*/
  getTempFftwShift().computeC2RFft1DX(getUxShifted());


  //-------------------------------------------------- uy shifted ----------------------------------------------------//
  getTempFftwShift().computeR2CFft1DY(getUySgy());

  #pragma omp parallel for schedule(static) if (simulationDimension == SD::k3D)
  for (size_t z = 0; z < yShiftDims.nz; z++)
  {
    #pragma omp parallel for schedule(static) if (simulationDimension == SD::k2D)
    for (size_t y = 0; y < yShiftDims.ny; y++)
    {
      #pragma omp simd
      for (size_t x = 0; x < yShiftDims.nx; x++)
      {
        const size_t i = get1DIndex(z, y, x, yShiftDims);

        tempFftShift[i] = (tempFftShift[i] * yShiftNegR[y]) * dividerY;
      } // x
    } // y
  }//z
  getTempFftwShift().computeC2RFft1DY(getUyShifted());

  //-------------------------------------------------- uz_shifted ----------------------------------------------------//
  if (simulationDimension == SD::k3D)
  {
    getTempFftwShift().computeR2CFft1DZ(getUzSgz());
    #pragma omp parallel for schedule(static)
    for (size_t z = 0; z < zShiftDims.nz; z++)
    {
      for (size_t y = 0; y < zShiftDims.ny; y++)
      {
        #pragma omp simd
        for (size_t x = 0; x < zShiftDims.nx; x++)
        {
          const size_t i = get1DIndex(z, y, x, zShiftDims);

          tempFftShift[i] = (tempFftShift[i] * zShiftNegR[z]) * dividerZ;
        } // x
      } // y
    }//z
    getTempFftwShift().computeC2RFft1DZ(getUzShifted());
  }
}// end of computeShiftedVelocity
//----------------------------------------------------------------------------------------------------------------------

/**
 * Print progress statistics.
 */
void KSpaceFirstOrderSolver::printStatistics()
{
  const size_t nt =  mParameters.getNt();
  const size_t timeIndex = mParameters.getTimeIndex();


  if (timeIndex > (mActPercent * nt * 0.01f))
  {
    mActPercent += mParameters.getProgressPrintInterval();

    mIterationTime.stop();

    const double elTime = mIterationTime.getElapsedTime();
    const double elTimeWithLegs = mIterationTime.getElapsedTime() + mSimulationTime.getElapsedTimeOverPreviousLegs();
    const double toGo   = ((elTimeWithLegs / static_cast<double>((timeIndex + 1)) *  nt)) - elTimeWithLegs;

    struct tm *current;
    time_t now;
    time(&now);
    now += toGo;
    current = localtime(&now);

    Logger::log(Logger::LogLevel::kBasic,
                kOutFmtSimulationProgress,
                static_cast<size_t>(((timeIndex) / (nt * 0.01f))),'%',
                elTime, toGo,
                current->tm_mday, current->tm_mon+1, current->tm_year-100,
                current->tm_hour, current->tm_min, current->tm_sec);
    Logger::flush(Logger::LogLevel::kBasic);
  }
}// end of printStatistics
//----------------------------------------------------------------------------------------------------------------------

/**
 * Was the loop interrupted to checkpoint?
 */
bool KSpaceFirstOrderSolver::isCheckpointInterruption() const
{
  return (mParameters.getTimeIndex() != mParameters.getNt());
}// end of isCheckpointInterruption
//----------------------------------------------------------------------------------------------------------------------

/**
 * Check the output file has the correct format and version.
 */
void KSpaceFirstOrderSolver::checkOutputFile()
{
  // The header has already been read
  Hdf5FileHeader& fileHeader = mParameters.getFileHeader();
  Hdf5File&       outputFile = mParameters.getOutputFile();

  // test file type
  if (fileHeader.getFileType() != Hdf5FileHeader::FileType::kOutput)
  {
    throw ios::failure(Logger::formatMessage(kErrFmtBadOutputFileFormat, mParameters.getOutputFileName().c_str()));
  }

  // test file major version
  if (!fileHeader.checkMajorFileVersion())
  {
    throw ios::failure(Logger::formatMessage(kErrFmtBadMajorFileVersion,
                                             mParameters.getOutputFileName().c_str(),
                                             fileHeader.getFileMajorVersion().c_str()));
  }

  // test file minor version
  if (!fileHeader.checkMinorFileVersion())
  {
    throw ios::failure(Logger::formatMessage(kErrFmtBadMinorFileVersion,
                                             mParameters.getOutputFileName().c_str(),
                                             fileHeader.getFileMinorVersion().c_str()));
  }


  // Check dimension sizes
  DimensionSizes outputDimSizes;
  outputFile.readScalarValue(outputFile.getRootGroup(), kNxName, outputDimSizes.nx);
  outputFile.readScalarValue(outputFile.getRootGroup(), kNyName, outputDimSizes.ny);
  outputFile.readScalarValue(outputFile.getRootGroup(), kNzName, outputDimSizes.nz);

 if (mParameters.getFullDimensionSizes() != outputDimSizes)
 {
    throw ios::failure(Logger::formatMessage(kErrFmtOutputDimensionsMismatch,
                                             outputDimSizes.nx,
                                             outputDimSizes.ny,
                                             outputDimSizes.nz,
                                             mParameters.getFullDimensionSizes().nx,
                                             mParameters.getFullDimensionSizes().ny,
                                             mParameters.getFullDimensionSizes().nz));
 }
}// end of checkOutputFile
//----------------------------------------------------------------------------------------------------------------------


/**
 * Check the file type and the version of the checkpoint file.
 */
void KSpaceFirstOrderSolver::checkCheckpointFile()
{
  // read the header and check the file version
  Hdf5FileHeader fileHeader;
  Hdf5File&      checkpointFile = mParameters.getCheckpointFile();

  fileHeader.readHeaderFromCheckpointFile(checkpointFile);

  // test file type
  if (fileHeader.getFileType() != Hdf5FileHeader::FileType::kCheckpoint)
  {
    throw ios::failure(Logger::formatMessage(kErrFmtBadCheckpointFileFormat,
                                             mParameters.getCheckpointFileName().c_str()));
  }

  // test file major version
  if (!fileHeader.checkMajorFileVersion())
  {
    throw ios::failure(Logger::formatMessage(kErrFmtBadMajorFileVersion,
                                             mParameters.getCheckpointFileName().c_str(),
                                             fileHeader.getFileMajorVersion().c_str()));
  }

  // test file minor version
  if (!fileHeader.checkMinorFileVersion())
  {
    throw ios::failure(Logger::formatMessage(kErrFmtBadMinorFileVersion,
                                             mParameters.getCheckpointFileName().c_str(),
                                             fileHeader.getFileMinorVersion().c_str()));
  }


  // Check dimension sizes
  DimensionSizes checkpointDimSizes;
  checkpointFile.readScalarValue(checkpointFile.getRootGroup(), kNxName, checkpointDimSizes.nx);
  checkpointFile.readScalarValue(checkpointFile.getRootGroup(), kNyName, checkpointDimSizes.ny);
  checkpointFile.readScalarValue(checkpointFile.getRootGroup(), kNzName, checkpointDimSizes.nz);

 if (mParameters.getFullDimensionSizes() != checkpointDimSizes)
 {
    throw ios::failure(Logger::formatMessage(kErrFmtCheckpointDimensionsMismatch,
                                             checkpointDimSizes.nx,
                                             checkpointDimSizes.ny,
                                             checkpointDimSizes.nz,
                                             mParameters.getFullDimensionSizes().nx,
                                             mParameters.getFullDimensionSizes().ny,
                                             mParameters.getFullDimensionSizes().nz));
 }
}// end of checkCheckpointFile
//----------------------------------------------------------------------------------------------------------------------

/**
 * Restore cumulated elapsed time from the output file.
 */
void KSpaceFirstOrderSolver::loadElapsedTimeFromOutputFile()
{
  double totalTime, dataLoadTime, preProcessingTime, simulationTime, postProcessingTime;

  // Get execution times stored in the output file header
  mParameters.getFileHeader().getExecutionTimes(totalTime,
                                                dataLoadTime,
                                                preProcessingTime,
                                                simulationTime,
                                                postProcessingTime);

  mTotalTime.SetElapsedTimeOverPreviousLegs(totalTime);
  mDataLoadTime.SetElapsedTimeOverPreviousLegs(dataLoadTime);
  mPreProcessingTime.SetElapsedTimeOverPreviousLegs(preProcessingTime);
  mSimulationTime.SetElapsedTimeOverPreviousLegs(simulationTime);
  mPostProcessingTime.SetElapsedTimeOverPreviousLegs(postProcessingTime);

}// end of loadElapsedTimeFromOutputFile
//----------------------------------------------------------------------------------------------------------------------

inline size_t KSpaceFirstOrderSolver::get1DIndex(const size_t          z,
                                                 const size_t          y,
                                                 const size_t          x,
                                                 const DimensionSizes& dimensionSizes)
{
  return (z * dimensionSizes.ny + y) * dimensionSizes.nx + x;
}// end of get1DIndex
//----------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Private methods --------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//
