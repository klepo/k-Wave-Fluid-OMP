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
 * @date        12 July      2012, 10:27 (created)\n
 *              30 August    2017, 18:11 (revised)
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
#endif

// Windows build
#ifdef _WIN64
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

#include <KSpaceSolver/KSpaceFirstOrder3DSolver.h>
#include <Containers/MatrixContainer.h>
#include <Containers/OutputStreamContainer.h>

#include <MatrixClasses/FftwComplexMatrix.h>
#include <Logger/Logger.h>

using std::ios;

//--------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------- Constants -----------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Public methods ---------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Constructor of the class.
 */
KSpaceFirstOrder3DSolver::KSpaceFirstOrder3DSolver():
        mMatrixContainer(), mOutputStreamContainer(),
        mParameters(Parameters::getInstance()),
        mActPercent(0l),
        mTotalTime(), mPreProcessingTime(), mDataLoadTime (), mSimulationTime(),
        mPostProcessingTime(), mIterationTime()
{
  mTotalTime.start();

  //Switch off HDF5 error messages
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);
}// end of KSpaceFirstOrder3DSolver
//----------------------------------------------------------------------------------------------------------------------


/**
 * Destructor of the class.
 */
KSpaceFirstOrder3DSolver::~KSpaceFirstOrder3DSolver()
{
  freeMemory();
}// end of KSpaceFirstOrder3DSolver
//----------------------------------------------------------------------------------------------------------------------

/**
 * The method allocates the matrix container, creates all matrices and creates all output streams
 * (however not allocating memory).
 */
void KSpaceFirstOrder3DSolver::allocateMemory()
{
  Logger::log(Logger::LogLevel::kBasic, kOutFmtMemoryAllocation);
  Logger::flush(Logger::LogLevel::kBasic);

  // create container, then all matrices
  mMatrixContainer.addMatrices();
  mMatrixContainer.createMatrices();

  // add output streams into container
  //@todo Think about moving under LoadInputData routine...
  mOutputStreamContainer.addStreams(mMatrixContainer);

  Logger::log(Logger::LogLevel::kBasic, kOutFmtDone);
}// end of allocateMemory
//----------------------------------------------------------------------------------------------------------------------

/**
 * The method frees all memory allocated by the class.
 */
void KSpaceFirstOrder3DSolver::freeMemory()
{
  mMatrixContainer.freeMatrices();
  mOutputStreamContainer.freeStreams();
}// end of freeMemory
//----------------------------------------------------------------------------------------------------------------------

/**
 * Load data from the input file provided by the Parameter class and creates the output time series streams.
 */
void KSpaceFirstOrder3DSolver::loadInputData()
{
  Logger::log(Logger::LogLevel::kBasic, kOutFmtDataLoading);
  Logger::flush(Logger::LogLevel::kBasic);
  // Start timer
  mDataLoadTime.start();

  // get handles
  Hdf5File& inputFile      = mParameters.getInputFile(); // file is opened (in Parameters)
  Hdf5File& outputFile     = mParameters.getOutputFile();
  Hdf5File& checkpointFile = mParameters.getCheckpointFile();

  // Load data from disk
  Logger::log(Logger::LogLevel::kFull, kOutFmtNoDone);
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

    checkpointFile.close();
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
  }
  else
  { //------------------------------------ First round of multi-leg simulation ---------------------------------------//
    // Create the output file
    Logger::log(Logger::LogLevel::kFull, kOutFmtCreatingOutputFile);
    Logger::flush(Logger::LogLevel::kFull);

    outputFile.create(mParameters.getOutputFileName());
    Logger::log(Logger::LogLevel::kFull, kOutFmtDone);

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
* This method computes k-space First Order 3D simulation.
 */
void KSpaceFirstOrder3DSolver::compute()
{
  // fft initialisation and preprocessing
  try
  {
    mPreProcessingTime.start();

    // initilaise all FFTW plans
    InitializeFftwPlans();

    // preprocessing phase generating necessary variables
    preProcessing();

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
  try
  {
    mSimulationTime.start();

    computeMainLoop();

    mSimulationTime.stop();

    Logger::log(Logger::LogLevel::kBasic,kOutFmtSimulationEndSeparator);
  }
  catch (const std::exception& e)
  {
    Logger::log(Logger::LogLevel::kBasic, kOutFmtSimulatoinFinalSeparator);
    Logger::errorAndTerminate(Logger::wordWrapString(e.what(),kErrFmtPathDelimiters, 9));
  }

  // Post processing region
  mPostProcessingTime.start();
  try
  {
    if (isCheckpointInterruption())
    { // Checkpoint
      Logger::log(Logger::LogLevel::kBasic, kOutFmtElapsedTime, mSimulationTime.getElapsedTime());
      Logger::log(Logger::LogLevel::kBasic, kOutFmtCheckpointTimeSteps, mParameters.getTimeIndex());
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
      Logger::log(Logger::LogLevel::kBasic, kOutFmtElapsedTime, mSimulationTime.getElapsedTime());
      Logger::log(Logger::LogLevel::kBasic, kOutFmtSeparator);
      Logger::log(Logger::LogLevel::kBasic, kOutFmtPostProcessing);
      Logger::flush(Logger::LogLevel::kBasic);

      postProcessing();

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
size_t KSpaceFirstOrder3DSolver::getMemoryUsage() const
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
 * Get release code version.
 */
std::string KSpaceFirstOrder3DSolver::getCodeName() const
{
  return std::string(kOutFmtKWaveVersion);
}// end of getCodeName
//----------------------------------------------------------------------------------------------------------------------


/**
 * Print full code name and the license.
 * @todo - Add __AVX512__ intrinsics
 */
void KSpaceFirstOrder3DSolver::printFullCodeNameAndLicense() const
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
void KSpaceFirstOrder3DSolver::setProcessorAffinity()
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
void KSpaceFirstOrder3DSolver::InitializeFftwPlans()
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

  // create real to complex plans
  getTempFftwX().createR2CFftPlan3D(getP());
  getTempFftwY().createR2CFftPlan3D(getP());
  getTempFftwZ().createR2CFftPlan3D(getP());

  // create real to complex plans
  getTempFftwX().createC2RFftPlan3D(getP());
  getTempFftwY().createC2RFftPlan3D(getP());
  getTempFftwZ().createC2RFftPlan3D(getP());

  // if necessary, create 1D shift plans.
  // in this case, the matrix has a bit bigger dimensions to be able to store
  // shifted matrices.
  if (Parameters::getInstance().getStoreVelocityNonStaggeredRawFlag())
  {
    // X shifts
    getTempFftwShift().createR2CFftPlan1DX(getP());
    getTempFftwShift().createC2RFftPlan1DX(getP());

    // Y shifts
    getTempFftwShift().createR2CFftPlan1DY(getP());
    getTempFftwShift().createC2RFftPlan1DY(getP());

    // Z shifts
    getTempFftwShift().createR2CFftPlan1DZ(getP());
    getTempFftwShift().createC2RFftPlan1DZ(getP());
  }// end u_non_staggered

  Logger::log(Logger::LogLevel::kBasic, kOutFmtDone);
}// end of InitializeFftwPlans
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute pre-processing phase.
 */
void KSpaceFirstOrder3DSolver::preProcessing()
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
      generateInitialDenisty();
    }
    else
    {
      getDtRho0Sgx().scalarDividedBy(mParameters.getDt());
      getDtRho0Sgy().scalarDividedBy(mParameters.getDt());
      getDtRho0Sgz().scalarDividedBy(mParameters.getDt());
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

  // calculate c^2. It has to be after kappa gen... because of c modification
  computeC2();

  Logger::log(Logger::LogLevel::kBasic, kOutFmtDone);
}// end of preProcessing
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute the main time loop of KSpaceFirstOrder3D.
 */
void KSpaceFirstOrder3DSolver::computeMainLoop()
{
  mActPercent = 0;
  // set ActPercent to correspond the t_index after recovery
  if (mParameters.getTimeIndex() > 0)
  {
    mActPercent = (100 * mParameters.getTimeIndex()) / mParameters.getNt();
  }

  // Progress header
  Logger::log(Logger::LogLevel::kBasic,kOutFmtSimulationHeader);

  mIterationTime.start();

  while ((mParameters.getTimeIndex() < mParameters.getNt()) && (!isTimeToCheckpoint()))
  {
    const size_t timeIndex = mParameters.getTimeIndex();

    // compute velocity
    computeVelocity();
    // add in the velocity source term
    addVelocitySource();

    // add in the transducer source term (t = t1) to ux
    if (mParameters.getTransducerSourceFlag() > timeIndex)
    {
     getUxSgx().addTransducerSource(getVelocitySourceIndex(),
                                    getTransducerSourceInput(),
                                    getDelayMask(),
                                    mParameters.getTimeIndex());
    }

    // compute gradient of velocity
    computeVelocityGradient();

    if (mParameters.getNonLinearFlag())
    {
      computeDensityNonliner();
    }
    else
    {
      computeDensityLinear();
    }


     // add in the source pressure term
     addPressureSource();

    if (mParameters.getNonLinearFlag())
    {
      computePressureNonlinear();
    }
    else
    {
      computePressureLinear();
    }

    // calculate initial pressure
    if ((timeIndex == 0) && (mParameters.getInitialPressureSourceFlag() == 1)) addInitialPressureSource();

    storeSensorData();
    printStatistics();
    mParameters.incrementTimeIndex();
  }// time loop
}// end of computeMainLoop
//----------------------------------------------------------------------------------------------------------------------

/**
 * Post processing the quantities, closing the output streams and storing the sensor mask.
 */
void KSpaceFirstOrder3DSolver::postProcessing()
{
  if (mParameters.getStorePressureFinalAllFlag())
  {
    getP().writeData(mParameters.getOutputFile(), kPressureFinalName, mParameters.getCompressionLevel());
  }// p_final

  if (mParameters.getStoreVelocityFinalAllFlag())
  {
    getUxSgx().writeData(mParameters.getOutputFile(), kUxFinalName, mParameters.getCompressionLevel());
    getUySgy().writeData(mParameters.getOutputFile(), kUyFinalName, mParameters.getCompressionLevel());
    getUzSgz().writeData(mParameters.getOutputFile(), kUzFinalName, mParameters.getCompressionLevel());
  }// u_final

  // Apply post-processing and close
  mOutputStreamContainer.postProcessStreams();
  mOutputStreamContainer.closeStreams();


  // store sensor mask if wanted
  if (mParameters.getCopySensorMaskFlag())
  {
    if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kIndex)
    {
      getSensorMaskIndex().recomputeIndicesToMatlab();
      getSensorMaskIndex().writeData(mParameters.getOutputFile(),kSensorMaskIndexName,
                                     mParameters.getCompressionLevel());
    }
    if (mParameters.getSensorMaskType() == Parameters::SensorMaskType::kCorners)
    {
      getSensorMaskCorners().recomputeIndicesToMatlab();
      getSensorMaskCorners().writeData(mParameters.getOutputFile(),kSensorMaskCornersName,
                                       mParameters.getCompressionLevel());
    }
  }
}// end of postProcessing
//----------------------------------------------------------------------------------------------------------------------

/**
 * Store sensor data.
 */
void KSpaceFirstOrder3DSolver::storeSensorData()
{
  // Unless the time for sampling has come, exit
  if (mParameters.getTimeIndex() >= mParameters.getSamplingStartTimeIndex())
  {
    if (mParameters.getStoreVelocityNonStaggeredRawFlag())
    {
      computeShiftedVelocity();
    }
    mOutputStreamContainer.sampleStreams();
  }
}// end of storeSensorData
//----------------------------------------------------------------------------------------------------------------------

/**
 * Write statistics and the header into the output file.
 */
void KSpaceFirstOrder3DSolver::writeOutputDataInfo()
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
}// end of writeOutputDataInfo
//----------------------------------------------------------------------------------------------------------------------

/**
 * Save checkpoint data into the checkpoint file, flush aggregated outputs into the output file.
 */
void KSpaceFirstOrder3DSolver::saveCheckpointData()
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

  // Close the checkpoint file
  checkpointFile.close();
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
}// end of saveCheckpointData
//----------------------------------------------------------------------------------------------------------------------


 /**
 * Compute new values of acoustic velocity in all three dimensions (UxSgx, UySgy, UzSgz).
 *
 * <b>Matlab code:</b> \n
 *
 * \verbatim
   p_k = fftn(p);
   ux_sgx = bsxfun(@times, pml_x_sgx, ...
       bsxfun(@times, pml_x_sgx, ux_sgx) ...
       - dt .* rho0_sgx_inv .* real(ifftn( bsxfun(@times, ddx_k_shift_pos, kappa .* fftn(p)) )) ...
       );
   uy_sgy = bsxfun(@times, pml_y_sgy, ...
       bsxfun(@times, pml_y_sgy, uy_sgy) ...
       - dt .* rho0_sgy_inv .* real(ifftn( bsxfun(@times, ddy_k_shift_pos, kappa .* fftn(p)) )) ...
       );
   uz_sgz = bsxfun(@times, pml_z_sgz, ...
       bsxfun(@times, pml_z_sgz, uz_sgz) ...
       - dt .* rho0_sgz_inv .* real(ifftn( bsxfun(@times, ddz_k_shift_pos, kappa .* fftn(p)) )) ...
       );
 \endverbatim
 */
 void KSpaceFirstOrder3DSolver::computeVelocity()
 {
  // bsxfun(@times, ddx_k_shift_pos, kappa .* fftn(p)), for all 3 dims
  computePressureGradient();

   getTempFftwX().computeC2RFft3D(getTemp1Real3D());
   getTempFftwY().computeC2RFft3D(getTemp2Real3D());
   getTempFftwZ().computeC2RFft3D(getTemp3Real3D());

  #pragma omp parallel
  {
    if (mParameters.getRho0ScalarFlag())
    { // scalars
      if (mParameters.getNonUniformGridFlag())
      {
        getUxSgx().computeVelocityXHomogeneousNonuniform(getTemp1Real3D(),
                                                         mParameters.getDtRho0SgxScalar(),
                                                         getDxudxnSgx(),
                                                         getPmlXSgx());
        getUySgy().computeVelocityYHomogeneousNonuniform(getTemp2Real3D(),
                                                         mParameters.getDtRho0SgyScalar(),
                                                         getDyudynSgy(),
                                                         getPmlYSgy());
        getUzSgz().computeVelocityZHomogeneousNonuniform(getTemp3Real3D(),
                                                         mParameters.getDtRho0SgzScalar(),
                                                         getDzudznSgz(),
                                                         getPmlZSgz());
       }
      else
      {
        getUxSgx().computeVelocityXHomogeneousUniform(getTemp1Real3D(),
                                                      mParameters.getDtRho0SgxScalar(),
                                                      getPmlXSgx());
        getUySgy().computeVelocityYHomogeneousUniform(getTemp2Real3D(),
                                                      mParameters.getDtRho0SgyScalar(),
                                                      getPmlYSgy());
        getUzSgz().computeVelocityZHomogeneousUniform(getTemp3Real3D(),
                                                      mParameters.getDtRho0SgzScalar(),
                                                      getPmlZSgz());
      }
    }
    else
    {// matrices
      getUxSgx().computeVelocityX(getTemp1Real3D(),
                                  getDtRho0Sgx(),
                                  getPmlXSgx());
      getUySgy().computeVelocityY(getTemp2Real3D(),
                                  getDtRho0Sgy(),
                                  getPmlYSgy());
      getUzSgz().computeVelocityZ(getTemp3Real3D(),
                                  getDtRho0Sgz(),
                                  getPmlZSgz());
    }
  } // parallel
}// end of computeVelocity
//----------------------------------------------------------------------------------------------------------------------

 /**
 * Compute new values for duxdx, duydy, duzdz.
 */
void  KSpaceFirstOrder3DSolver::computeVelocityGradient()
{
  getTempFftwX().computeR2CFft3D(getUxSgx());
  getTempFftwY().computeR2CFft3D(getUySgy());
  getTempFftwZ().computeR2CFft3D(getUzSgz());

  #pragma omp parallel
  {
    float* tempFftX = getTempFftwX().getData();
    float* tempFftY = getTempFftwY().getData();
    float* tempFftZ = getTempFftwZ().getData();

    const float* kappa   = getKappa().getData();

    const size_t fftDimZ = getTempFftwX().getDimensionSizes().nz;
    const size_t fftDimY = getTempFftwX().getDimensionSizes().ny;
    const size_t fftDimX = getTempFftwX().getDimensionSizes().nx;

    const size_t slabSize = (fftDimX * fftDimY) << 1;
    const float  divider = 1.0f / static_cast<float>(getUxSgx().size());

    const FloatComplex* ddx = reinterpret_cast<FloatComplex*>(getDdxKShiftNeg().getData());
    const FloatComplex* ddy = reinterpret_cast<FloatComplex*>(getDdyKShiftNeg().getData());
    const FloatComplex* ddz = reinterpret_cast<FloatComplex*>(getDdzKShiftNeg().getData());


    #pragma omp for schedule (static)
    for (size_t z = 0; z < fftDimZ; z++)
    {
      register size_t i = z * slabSize;

      const float ddzNegRe = ddz[z].real;
      const float ddzNegIm = ddz[z].imag;
      for (size_t y = 0; y < fftDimY; y++)
      {
        const float ddyNegRe = ddy[y].real;
        const float ddyNegIm = ddy[y].imag;
        for (size_t x = 0; x < fftDimX; x++)
        {
          const float eKappa = kappa[i >> 1];

          const float fftXRe = tempFftX[i]   *= eKappa;
          const float ffyXIm = tempFftX[i+1] *= eKappa;

          const float fftYre = tempFftY[i]   *= eKappa;
          const float fftYIm = tempFftY[i+1] *= eKappa;

          const float fftZRe = tempFftZ[i]   *= eKappa;
          const float fftZIm = tempFftZ[i+1] *= eKappa;

          tempFftX[i]     = ((fftXRe * ddx[x].real) - (ffyXIm * ddx[x].imag) ) * divider;
          tempFftX[i + 1] = ((ffyXIm * ddx[x].real) + (fftXRe * ddx[x].imag) ) * divider;

          tempFftY[i]     = ((fftYre * ddyNegRe) - (fftYIm * ddyNegIm)) * divider;
          tempFftY[i + 1] = ((fftYIm * ddyNegRe) + (fftYre * ddyNegIm)) * divider;

          tempFftZ[i]     = ((fftZRe * ddzNegRe) - (fftZIm * ddzNegIm)) * divider;
          tempFftZ[i + 1] = ((fftZIm * ddzNegRe) + (fftZRe * ddzNegIm)) * divider;

          i+=2;
        } // x
      } // y
    } // z
  } // parallel;

  getTempFftwX().computeC2RFft3D(getDuxdx());
  getTempFftwY().computeC2RFft3D(getDuydy());
  getTempFftwZ().computeC2RFft3D(getDuzdz());

 //------------------------------------------------- Non linear grid -------------------------------------------------//
  if (mParameters.getNonUniformGridFlag() != 0)
  {
    #pragma omp parallel
    {
      float* duxdx = getDuxdx().getData();
      float* duydy = getDuydy().getData();
      float* duzdz = getDuzdz().getData();

      const float* duxdxn = getDxudxn().getData();
      const float* duydyn = getDyudyn().getData();
      const float* duzdzn = getDzudzn().getData();

      const size_t nz = getDuxdx().getDimensionSizes().nz;
      const size_t ny = getDuxdx().getDimensionSizes().ny;
      const size_t nx = getDuxdx().getDimensionSizes().nx;

      const size_t SliceSize = (nx * ny);

      #pragma omp for schedule (static)
      for (size_t z = 0; z < nz; z++)
      {
        register size_t i = z * SliceSize;
        for (size_t y = 0; y < ny; y++)
        {
          for (size_t x = 0; x < nx; x++)
          {
            duxdx[i] *= duxdxn[x];
            i++;
          } // x
        } // y
      } // z

      #pragma omp for schedule (static)
      for (size_t z = 0; z < nz; z++)
      {
        register size_t i = z * SliceSize;
        for (size_t y = 0; y < ny; y++)
        {
          const float eDyudyn = duydyn[y];
          for (size_t x = 0; x < nx; x++)
          {
            duydy[i] *= eDyudyn;
            i++;
          } // x
        } // y
      } // z

      #pragma omp for schedule (static)
      for (size_t z = 0; z < nz; z++)
      {
        const float eDuzdzn = duzdzn[z];
        register size_t i = z * SliceSize;
        for (size_t y = 0; y < ny; y++)
        {
          for (size_t x = 0; x < nx; x++)
          {
            duzdz[i] *=  eDuzdzn;
            i++;
          } // x
        } // y
      } // z
    } // parallel
 }// nonlinear
}// end of computeVelocityGradient
//----------------------------------------------------------------------------------------------------------------------

/**
 * Calculate new values of acoustic density for nonlinear case (rhoX, rhoy and rhoZ).
 *
 * <b>Matlab code:</b> \n
 *
 *\verbatim
    rho0_plus_rho = 2 .* (rhox + rhoy + rhoz) + rho0;
    rhox = bsxfun(@times, pml_x, bsxfun(@times, pml_x, rhox) - dt .* rho0_plus_rho .* duxdx);
    rhoy = bsxfun(@times, pml_y, bsxfun(@times, pml_y, rhoy) - dt .* rho0_plus_rho .* duydy);
    rhoz = bsxfun(@times, pml_z, bsxfun(@times, pml_z, rhoz) - dt .* rho0_plus_rho .* duzdz);
 \endverbatim
 */
void KSpaceFirstOrder3DSolver::computeDensityNonliner()
{
  const size_t nz = getRhoX().getDimensionSizes().nz;
  const size_t ny = getRhoX().getDimensionSizes().ny;
  const size_t nx = getRhoX().getDimensionSizes().nx;

  const float dt  = mParameters.getDt();
  const size_t slabSize = ny * nx;

  #pragma omp parallel
  {
    float* rhoX  = getRhoX().getData();
    float* rhoY  = getRhoY().getData();
    float* rhoZ  = getRhoZ().getData();

    const float* pmlX  = getPmlX().getData();
    const float* pmlY  = getPmlY().getData();
    const float* pmlZ  = getPmlZ().getData();

    const float* duxdx = getDuxdx().getData();
    const float* duydy = getDuydy().getData();
    const float* duzdz = getDuzdz().getData();

    //----------------------------------------------- rho0 is scalar -------------------------------------------------//
    if (mParameters.getRho0ScalarFlag())
    {
      const float rho0 = mParameters.getRho0Scalar();

      #pragma omp for schedule(static)
      for (size_t z = 0; z < nz; z++)
      {
        size_t i = z * slabSize;
        for (size_t y = 0; y < ny; y++)
        {
          for (size_t x = 0; x < nx; x++)
          {
            const float sumRhosDt = (2.0f * (rhoX[i] + rhoY[i] + rhoZ[i]) + rho0) * dt;

            rhoX[i] = pmlX[x] * ((pmlX[x] * rhoX[i]) - sumRhosDt * duxdx[i]);
            rhoY[i] = pmlY[y] * ((pmlY[y] * rhoY[i]) - sumRhosDt * duydy[i]);
            rhoZ[i] = pmlZ[z] * ((pmlZ[z] * rhoZ[i]) - sumRhosDt * duzdz[i]);

            i++;
          }// x
        }// y
      }// z
    }
    else
    { //---------------------------------------------- rho0 is matrix ------------------------------------------------//
      // rho0 is a matrix
      const float* rho0  = getRho0().getData();

      #pragma omp for schedule(static)
      for (size_t z = 0; z < nz; z++)
      {
        size_t i = z * slabSize;
        for (size_t y = 0; y < ny; y++)
        {
          for (size_t x = 0; x < nx; x++)
          {
            const float sumRhosDt = (2.0f * (rhoX[i] + rhoY[i] + rhoZ[i]) + rho0[i]) * dt;

            rhoX[i] = pmlX[x] * ((pmlX[x] * rhoX[i]) - sumRhosDt * duxdx[i]);
            rhoY[i] = pmlY[y] * ((pmlY[y] * rhoY[i]) - sumRhosDt * duydy[i]);
            rhoZ[i] = pmlZ[z] * ((pmlZ[z] * rhoZ[i]) - sumRhosDt * duzdz[i]);

            i++;
          } // x
        }// y
      }// z
    } // end rho is matrix
  }// parallel
}// end of computeDensityNonliner
//----------------------------------------------------------------------------------------------------------------------




/**
 * Calculate new values of acoustic density for linear case (rhoX, rhoy and rhoZ).
 *
 * <b>Matlab code:</b> \n
 *
 *\verbatim
    rhox = bsxfun(@times, pml_x, bsxfun(@times, pml_x, rhox) - dt .* rho0 .* duxdx);
    rhoy = bsxfun(@times, pml_y, bsxfun(@times, pml_y, rhoy) - dt .* rho0 .* duydy);
    rhoz = bsxfun(@times, pml_z, bsxfun(@times, pml_z, rhoz) - dt .* rho0 .* duzdz);
\endverbatim
 *
 */
void KSpaceFirstOrder3DSolver::computeDensityLinear()
{
  const size_t nz = getRhoX().getDimensionSizes().nz;
  const size_t ny = getRhoX().getDimensionSizes().ny;
  const size_t nx = getRhoX().getDimensionSizes().nx;

  const float dt        = mParameters.getDt();
  const size_t slabSize =  ny * nx;

  #pragma omp parallel
  {
    float* rhox  = getRhoX().getData();
    float* rhoy  = getRhoY().getData();
    float* rhoz  = getRhoZ().getData();

    const float* pmlX  = getPmlX().getData();
    const float* pmlY  = getPmlY().getData();
    const float* pmlZ  = getPmlZ().getData();

    const float* duxdx = getDuxdx().getData();
    const float* duydy = getDuydy().getData();
    const float* duzdz = getDuzdz().getData();

    //----------------------------------------------- rho0 is scalar -------------------------------------------------//
    if (mParameters.getRho0ScalarFlag())
    { // rho0 is a scalar
      const float dtRho0 = mParameters.getRho0Scalar() * dt;

      #pragma omp for schedule (static)
      for (size_t z = 0; z < nz; z++)
      {
        register size_t i = z * slabSize;
        for (size_t y = 0; y < ny; y++)
        {
          for (size_t x = 0; x < nx; x++)
          {
            const float ePmlX   = pmlX[x];

            rhox[i] = ePmlX * (((ePmlX * rhox[i]) - (dtRho0 * duxdx[i])) );
            i++;
          } // x
        }// y
      }// z

      #pragma omp for schedule (static)
      for (size_t z = 0; z < nz; z++)
      {
        register size_t i = z * slabSize;
        for (size_t y = 0; y < ny; y++)
        {
          const float ePmlY = pmlY[y];
          for (size_t x = 0; x < nx; x++)
          {
            rhoy[i] = ePmlY * (((ePmlY * rhoy[i]) - (dtRho0 * duydy[i])));
            i++;
          } // x
        }// y
      }// z

      #pragma omp for schedule (static)
      for (size_t z = 0; z < nz; z++)
      {
        register size_t i = z * slabSize;
        const float ePmlZ = pmlZ[z];

        for (size_t y = 0; y < ny; y++)
        {
          for (size_t x = 0; x < nx; x++)
          {
            rhoz[i] = ePmlZ * (((ePmlZ * rhoz[i]) - (dtRho0 * duzdz[i])));
            i++;
          } // x
        }// y
      }// z

    }
    else
    { //---------------------------------------------- rho0 is matrix ------------------------------------------------//
      // rho0 is a matrix
      const float* rho0  = getRho0().getData();

      #pragma omp for schedule (static)
      for (size_t z = 0; z < nz; z++)
      {
        register size_t i = z * slabSize;
        for (size_t y = 0; y < ny; y++)
        {
          for (size_t x = 0; x < nx; x++)
          {
            const float ePmlX   = pmlX[x];
            const float dtRho0 = dt * rho0[i];

            rhox[i] = ePmlX * (((ePmlX * rhox[i]) - (dtRho0 * duxdx[i])));

            i++;
          } // x
        }// y
      }// z

      #pragma omp for schedule (static)
      for (size_t z = 0; z < nz; z++)
      {
        register size_t i = z * slabSize;
        for (size_t y = 0; y < ny; y++)
        {
          const float ePmlY = pmlY[y];
          for (size_t x = 0; x < nx; x++)
          {
            const float dtRho0 = dt * rho0[i];

            rhoy[i] = ePmlY * (((ePmlY * rhoy[i]) - (dtRho0 * duydy[i])));
            i++;

          } // x
        }// y
      }// z


      #pragma omp for schedule (static)
      for (size_t z = 0; z < nz; z++)
      {
        register size_t i = z * slabSize;
        const float ePmlZ = pmlZ[z];

        for (size_t y = 0; y < ny; y++)
        {
          for (size_t x = 0; x < nx; x++)
          {
            const float dtRho0 = dt * rho0[i];

            rhoz[i] = ePmlZ * (((ePmlZ * rhoz[i]) - (dtRho0 * duzdz[i])));
            i++;
          } // x
        }// y
      }// z

   } // end rho is a matrix
  }// parallel
}// end of computeDensityLinear
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute acoustic pressure for non-linear case.
 *
 * <b>Matlab code:</b> \n
 *
 *\verbatim
    case 'lossless'
        % calculate p using a nonlinear adiabatic equation of state
        p = c.^2 .* (rhox + rhoy + rhoz + medium.BonA .* (rhox + rhoy + rhoz).^2 ./ (2 .* rho0));

    case 'absorbing'
        % calculate p using a nonlinear absorbing equation of state
        p = c.^2 .* (...
            (rhox + rhoy + rhoz) ...
            + absorb_tau .* real(ifftn( absorb_nabla1 .* fftn(rho0 .* (duxdx + duydy + duzdz)) ))...
            - absorb_eta .* real(ifftn( absorb_nabla2 .* fftn(rhox + rhoy + rhoz) ))...
            + medium.BonA .*(rhox + rhoy + rhoz).^2 ./ (2 .* rho0) ...
            );

 \endverbatim
 */
 void KSpaceFirstOrder3DSolver::computePressureNonlinear()
{
  if (mParameters.getAbsorbingFlag())
  { // absorbing case

    RealMatrix& densitySum         = getTemp1Real3D();
    RealMatrix& nonlinearTerm      = getTemp2Real3D();
    RealMatrix& velocitGradientSum = getTemp3Real3D();

    // reusing of the temp variables
    RealMatrix& absorbTauTerm = velocitGradientSum;
    RealMatrix& absorbEtaTerm = densitySum;


    computePressureTermsNonlinearSSE2(densitySum, nonlinearTerm, velocitGradientSum);

    // ifftn( absorb_nabla1 * fftn (rho0 * (duxdx+duydy+duzdz))
    getTempFftwX().computeR2CFft3D(velocitGradientSum);
    getTempFftwY().computeR2CFft3D(densitySum);

    computeAbsorbtionTermSSE2(getTempFftwX(), getTempFftwY());

    getTempFftwX().computeC2RFft3D(absorbTauTerm);
    getTempFftwY().computeC2RFft3D(absorbEtaTerm);

    sumPressureTermsNonlinear(absorbTauTerm, absorbEtaTerm, nonlinearTerm);
  }
  else
  {
    sumPressureTermsNonlinearLossless();
  }
}// end of computePressureNonlinear
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute new p for linear case.
 *
 * <b>Matlab code:</b> \n
 *
 *\verbatim
    case 'lossless'

        % calculate p using a linear adiabatic equation of state
        p = c.^2 .* (rhox + rhoy + rhoz);

    case 'absorbing'

        % calculate p using a linear absorbing equation of state
        p = c.^2 .* ( ...
            (rhox + rhoy + rhoz) ...
            + absorb_tau .* real(ifftn( absorb_nabla1 .* fftn(rho0 .* (duxdx + duydy + duzdz)) )) ...
            - absorb_eta .* real(ifftn( absorb_nabla2 .* fftn(rhox + rhoy + rhoz) )) ...
            );
 \endverbatim
 */
 void KSpaceFirstOrder3DSolver::computePressureLinear()
 {
  // rhox + rhoy + rhoz
  if (mParameters.getAbsorbingFlag())
  { // absorbing case

    RealMatrix& densitySum           = getTemp1Real3D();
    RealMatrix& velocityGradientTerm = getTemp2Real3D();

    RealMatrix& absorbTauTerm        = getTemp2Real3D();
    RealMatrix& absorbEtaTerm        = getTemp3Real3D();

    computePressureTermsLinear(densitySum, velocityGradientTerm);

    // ifftn ( absorb_nabla1 * fftn (rho0 * (duxdx+duydy+duzdz))

    getTempFftwX().computeR2CFft3D(velocityGradientTerm);
    getTempFftwY().computeR2CFft3D(densitySum);

    computeAbsorbtionTermSSE2(getTempFftwX(), getTempFftwY());

    getTempFftwX().computeC2RFft3D(absorbTauTerm);
    getTempFftwY().computeC2RFft3D(absorbEtaTerm);

    sumPressureTermsLinear(absorbTauTerm, absorbEtaTerm, densitySum);
  }
  else
  {
    // lossless case
    sumPressureTermsLinearLossless();
  }
 }// end of computePressureLinear
//----------------------------------------------------------------------------------------------------------------------

/**
 * Add u source to the particle velocity.
 */
void KSpaceFirstOrder3DSolver::addVelocitySource()
{
  const size_t timeIndex = mParameters.getTimeIndex();

  if (mParameters.getVelocityXSourceFlag() > timeIndex)
  {
    getUxSgx().addVelocitySource(GetVelocityXSourceInput(),
                                 getVelocitySourceIndex(),
                                 timeIndex,
                                 mParameters.getVelocitySourceMode(),
                                 mParameters.getVelocitySourceMany());
  }

  if (mParameters.getVelocityYSourceFlag() > timeIndex)
  {
    getUySgy().addVelocitySource(GetVelocityYSourceInput(),
                                 getVelocitySourceIndex(),
                                 timeIndex,
                                 mParameters.getVelocitySourceMode(),
                                 mParameters.getVelocitySourceMany());
  }

  if (mParameters.getVelocityZSourceFlag() > timeIndex)
  {
    getUzSgz().addVelocitySource(getVelocityZSourceInput(),
                                 getVelocitySourceIndex(),
                                 timeIndex,
                                 mParameters.getVelocitySourceMode(),
                                 mParameters.getVelocitySourceMany());
  }
}// end of addVelocitySource
//----------------------------------------------------------------------------------------------------------------------

 /**
  * Add in pressure source.
  */
void KSpaceFirstOrder3DSolver::addPressureSource()
{
  const size_t timeIndex = mParameters.getTimeIndex();

  if (mParameters.getPressureSourceFlag() > timeIndex)
  {
    float* rhox = getRhoX().getData();
    float* rhoy = getRhoY().getData();
    float* rhoz = getRhoZ().getData();

    const float*  sourceInput = getPressureSourceInput().getData();
    const size_t* sourceIndex = getPressureSourceIndex().getData();

    const bool   isManyFlag  = (mParameters.getPressureSourceMany() != 0);
    const size_t sourceSize  = getPressureSourceIndex().size();
    const size_t index2D     = (isManyFlag) ? timeIndex * sourceSize : timeIndex;

    // replacement
    if (mParameters.getPressureSourceMode() == 0)
    {
      #pragma omp parallel for if (sourceSize > 16384)
      for (size_t i = 0; i < sourceSize; i++)
      {
        const size_t signalIndex = (isManyFlag) ? index2D + i : index2D;

        rhox[sourceIndex[i]] = sourceInput[signalIndex];
        rhoy[sourceIndex[i]] = sourceInput[signalIndex];
        rhoz[sourceIndex[i]] = sourceInput[signalIndex];
      }
    }
    // Addition
    else
    {
      #pragma omp parallel for if (sourceSize > 16384)
      for (size_t i = 0; i < sourceSize; i++)
      {
        const size_t signalIndex = (isManyFlag) ? index2D + i : index2D;

        rhox[sourceIndex[i]] += sourceInput[signalIndex];
        rhoy[sourceIndex[i]] += sourceInput[signalIndex];
        rhoz[sourceIndex[i]] += sourceInput[signalIndex];
      }
    }// type of replacement
  }// if do at all
}// end of addPressureSource
//----------------------------------------------------------------------------------------------------------------------

/**
 * Calculate p0 source when necessary.
 *
 * <b>Matlab code:</b> \n
 *
 *\verbatim
    % add the initial pressure to rho as a mass source
    p = source.p0;
    rhox = source.p0 ./ (3 .* c.^2);
    rhoy = source.p0 ./ (3 .* c.^2);
    rhoz = source.p0 ./ (3 .* c.^2);

    % compute u(t = t1 + dt/2) based on the assumption u(dt/2) = -u(-dt/2)
    % which forces u(t = t1) = 0
    ux_sgx = dt .* rho0_sgx_inv .* real(ifftn( bsxfun(@times, ddx_k_shift_pos, kappa .* fftn(p)) )) / 2;
    uy_sgy = dt .* rho0_sgy_inv .* real(ifftn( bsxfun(@times, ddy_k_shift_pos, kappa .* fftn(p)) )) / 2;
    uz_sgz = dt .* rho0_sgz_inv .* real(ifftn( bsxfun(@times, ddz_k_shift_pos, kappa .* fftn(p)) )) / 2;
 \endverbatim
 */
void KSpaceFirstOrder3DSolver::addInitialPressureSource()
{
  getP().copyData(getInitialPressureSourceInput());

  const float* sourceInput = getInitialPressureSourceInput().getData();

  const bool   c0ScalarFlag = mParameters.getC0ScalarFlag();
  const float  c2Scalar     = (c0ScalarFlag) ? mParameters.getC2Scalar() : 0;
  const float* c2Matrix     = (c0ScalarFlag) ? nullptr : getC2().getData();

  float* rhox = getRhoX().getData();
  float* rhoy = getRhoY().getData();
  float* rhoz = getRhoZ().getData();

  const size_t size = getRhoX().size();

  #pragma omp parallel for schedule (static)
  for (size_t i = 0; i < size; i++)
  {
    const float tmp = sourceInput[i] / (3.0f * ((c0ScalarFlag) ? c2Scalar : c2Matrix[i]));
    rhox[i] = tmp;
    rhoy[i] = tmp;
    rhoz[i] = tmp;
  }

  //------------------------------------------------------------------------//
  //--  compute u(t = t1 + dt/2) based on the assumption u(dt/2) = -u(-dt/2) --//
  //--    which forces u(t = t1) = 0 --//
  //------------------------------------------------------------------------//
  computePressureGradient();

  if (mParameters.getRho0ScalarFlag())
  {
    if (mParameters.getNonUniformGridFlag())
    { // non uniform grid
      getUxSgx().computeInitialVelocityXHomogeneousNonuniform(mParameters.getDtRho0SgxScalar(),
                                                              getDxudxnSgx(),
                                                              getTempFftwX());
      getUySgy().computeInitialVelocityYHomogeneousNonuniform(mParameters.getDtRho0SgyScalar(),
                                                              getDyudynSgy(),
                                                              getTempFftwY());
      getUzSgz().computeInitialVelocityZHomogeneousNonuniform(mParameters.getDtRho0SgzScalar(),
                                                              getDzudznSgz(),
                                                              getTempFftwZ());
    }
    else
    { //uniform grid, heterogeneous
      getUxSgx().computeInitialVelocityHomogeneousUniform(mParameters.getDtRho0SgxScalar(), getTempFftwX());
      getUySgy().computeInitialVelocityHomogeneousUniform(mParameters.getDtRho0SgyScalar(), getTempFftwY());
      getUzSgz().computeInitialVelocityHomogeneousUniform(mParameters.getDtRho0SgzScalar(), getTempFftwZ());
    }
  }
  else
  { // homogeneous, unifrom grid
    // divide the matrix by 2 and multiply with st./rho0_sg
    getUxSgx().computeInitialVelocity(getDtRho0Sgx(), getTempFftwX());
    getUySgy().computeInitialVelocity(getDtRho0Sgy(), getTempFftwY());
    getUzSgz().computeInitialVelocity(getDtRho0Sgz(), getTempFftwZ());
  }
}// end of addInitialPressureSource
//----------------------------------------------------------------------------------------------------------------------

/**
 * Generate kappa matrix for lossless medium.
 */
void KSpaceFirstOrder3DSolver::generateKappa()
{
  #pragma omp parallel
  {
    const float dx2Rec = 1.0f / (mParameters.getDx() * mParameters.getDx());
    const float dy2Rec = 1.0f / (mParameters.getDy() * mParameters.getDy());
    const float dz2Rec = 1.0f / (mParameters.getDz() * mParameters.getDz());

    const float cRefDtPi = mParameters.getCRef() * mParameters.getDt() * static_cast<float>(M_PI);

    const float nxRec = 1.0f / static_cast<float>(mParameters.getFullDimensionSizes().nx);
    const float nyRec = 1.0f / static_cast<float>(mParameters.getFullDimensionSizes().ny);
    const float nzRec = 1.0f / static_cast<float>(mParameters.getFullDimensionSizes().nz);

    const size_t nx = mParameters.getReducedDimensionSizes().nx;
    const size_t ny = mParameters.getReducedDimensionSizes().ny;
    const size_t nz = mParameters.getReducedDimensionSizes().nz;

    float* kappa = getKappa().getData();

    #pragma omp for schedule (static)
    for (size_t z = 0; z < nz; z++)
    {
      const float zf    = static_cast<float>(z);
            float zPart = 0.5f - fabs(0.5f - zf * nzRec);
                  zPart = (zPart * zPart) * dz2Rec;

      for (size_t y = 0; y < ny; y++)
      {
        const float yf    = static_cast<float>(y);
              float yPart = 0.5f - fabs(0.5f - yf * nyRec);
                    yPart = (yPart * yPart) * dy2Rec;

        const float yzPart = zPart + yPart;
        for (size_t x = 0; x < nx; x++)
        {
          const float xf = static_cast<float>(x);
                float xPart = 0.5f - fabs(0.5f - xf * nxRec);
                      xPart = (xPart * xPart) * dx2Rec;

                float k = cRefDtPi * sqrt(xPart + yzPart);

          // kappa element
          kappa[(z * ny + y) * nx + x ] = (k == 0.0f) ? 1.0f : sin(k) / k;
        }//x
      }//y
    }// z
  }// parallel
}// end of generateKappa
//----------------------------------------------------------------------------------------------------------------------

/**
 * Generate kappa, absorb_nabla1, absorb_nabla2 for absorbing medium.
 */
void KSpaceFirstOrder3DSolver::generateKappaAndNablas()
{
  #pragma omp parallel
  {
    const float dxSqRec    = 1.0f / (mParameters.getDx() * mParameters.getDx());
    const float dySqRec    = 1.0f / (mParameters.getDy() * mParameters.getDy());
    const float dzSqRec    = 1.0f / (mParameters.getDz() * mParameters.getDz());

    const float cRefDt2    = mParameters.getCRef() * mParameters.getDt() * 0.5f;
    const float pi2        = static_cast<float>(M_PI) * 2.0f;

    const size_t nx        = mParameters.getFullDimensionSizes().nx;
    const size_t ny        = mParameters.getFullDimensionSizes().ny;
    const size_t nz        = mParameters.getFullDimensionSizes().nz;

    const float nxRec      = 1.0f / static_cast<float>(nx);
    const float nyRec      = 1.0f / static_cast<float>(ny);
    const float nzRec      = 1.0f / static_cast<float>(nz);

    const size_t nxComplex = mParameters.getReducedDimensionSizes().nx;
    const size_t nyComplex = mParameters.getReducedDimensionSizes().ny;
    const size_t nzComplex = mParameters.getReducedDimensionSizes().nz;

    float* kappa           = getKappa().getData();
    float* absorbNabla1    = getAbsorbNabla1().getData();
    float* absorbNabla2    = getAbsorbNabla2().getData();
    const float alphaPower = mParameters.getAlphaPower();

    #pragma omp for schedule (static)
    for (size_t z = 0; z < nzComplex; z++)
    {
      const float zf    = static_cast<float>(z);
            float zPart = 0.5f - fabs(0.5f - zf * nzRec);
                  zPart = (zPart * zPart) * dzSqRec;

      for (size_t y = 0; y < nyComplex; y++)
      {
        const float yf    = static_cast<float>(y);
              float yPart = 0.5f - fabs(0.5f - yf * nyRec);
                    yPart = (yPart * yPart) * dySqRec;

        const float yzPart = zPart + yPart;

        size_t i = (z * nyComplex + y) * nxComplex;

        for (size_t x = 0; x < nxComplex; x++)
        {
          const float xf    = static_cast<float>(x);
                float xPart = 0.5f - fabs(0.5f - xf * nxRec);
                      xPart = (xPart * xPart) * dxSqRec;

                float k     = pi2 * sqrt(xPart + yzPart);
                float cRefK = cRefDt2 * k;

          kappa[i]          = (cRefK == 0.0f) ? 1.0f : sin(cRefK) / cRefK;

          absorbNabla1[i] = pow(k, alphaPower - 2);
          absorbNabla2[i] = pow(k, alphaPower - 1);

          if (absorbNabla1[i] ==  std::numeric_limits<float>::infinity()) absorbNabla1[i] = 0.0f;
          if (absorbNabla2[i] ==  std::numeric_limits<float>::infinity()) absorbNabla2[i] = 0.0f;

          i++;
        }//x
      }//y
    }// z
  }// parallel
}// end of generateKappaAndNablas
//----------------------------------------------------------------------------------------------------------------------

/**
 * Generate absorbTau and absorbEta in for heterogenous medium.
 */
void KSpaceFirstOrder3DSolver::generateTauAndEta()
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
    #pragma omp parallel
    {
      const size_t nx  = mParameters.getFullDimensionSizes().nx;
      const size_t ny  = mParameters.getFullDimensionSizes().ny;
      const size_t nz  = mParameters.getFullDimensionSizes().nz;

      float* absorbTau = getAbsorbTau().getData();
      float* absorbEta = getAbsorbEta().getData();

      const bool   alphaCoeffScalarFlag = mParameters.getAlphaCoeffScalarFlag();
      const float  alphaCoeffScalar     = (alphaCoeffScalarFlag) ? mParameters.getAlphaCoeffScalar() : 0;
      const float* alphaCoeffMatrix     = (alphaCoeffScalarFlag) ? nullptr : getTemp1Real3D().getData();


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


      #pragma omp for schedule (static)
      for (size_t z = 0; z < nz; z++)
      {
        for (size_t y = 0; y < ny; y++)
        {
          size_t i = (z * ny + y) * nx;
          for (size_t x = 0; x < nx; x++)
          {
            const float alphaCoeff2 = 2.0f * alphaNeperCoeff *
                                      ((alphaCoeffScalarFlag) ? alphaCoeffScalar : alphaCoeffMatrix[i]);

            absorbTau[i] = (-alphaCoeff2) * pow(((c0ScalarFlag) ? c0Scalar : cOMatrix[i]), alphaPower - 1);
            absorbEta[i] =   alphaCoeff2  * pow(((c0ScalarFlag) ? c0Scalar : cOMatrix[i]),
                                                  alphaPower) * tanPi2AlphaPower;

            i++;
          }//x
        }//y
      }// z
    }// parallel
  } // matrix
}// end of generateTauAndEta
//----------------------------------------------------------------------------------------------------------------------

/**
 * Prepare dt./ rho0  for non-uniform grid.
 */
void KSpaceFirstOrder3DSolver::generateInitialDenisty()
{
  #pragma omp parallel
  {
    float* dtRho0Sgx   = getDtRho0Sgx().getData();
    float* dtRho0Sgy   = getDtRho0Sgy().getData();
    float* dtRho0Sgz   = getDtRho0Sgz().getData();

    const float dt = mParameters.getDt();

    const float* duxdxnSgx = getDxudxnSgx().getData();
    const float* duydynSgy = getDyudynSgy().getData();
    const float* duzdznSgz = getDzudznSgz().getData();

    const size_t nz = getDtRho0Sgx().getDimensionSizes().nz;
    const size_t ny = getDtRho0Sgx().getDimensionSizes().ny;
    const size_t nx = getDtRho0Sgx().getDimensionSizes().nx;

    const size_t sliceSize = (nx * ny);

    #pragma omp for schedule (static)
    for (size_t z = 0; z < nz; z++)
    {
      register size_t i = z * sliceSize;
      for (size_t y = 0; y < ny; y++)
      {
        for (size_t x = 0; x < nx; x++)
        {
          dtRho0Sgx[i] = (dt * duxdxnSgx[x]) / dtRho0Sgx[i];
          i++;
        } // x
      } // y
    } // z

    #pragma omp for schedule (static)
    for (size_t z = 0; z < nz; z++)
    {
      register size_t i = z * sliceSize;
      for (size_t y = 0; y < ny; y++)
      {
        const float duydynEl = duydynSgy[y];
        for (size_t x = 0; x < nx; x++)
        {
          dtRho0Sgy[i] = (dt * duydynEl) / dtRho0Sgy[i];
          i++;
        } // x
      } // y
    } // z

    #pragma omp for schedule (static)
    for (size_t z = 0; z < nz; z++)
    {
      register size_t i = z * sliceSize;
      const float duzdznEl = duzdznSgz[z];
      for (size_t y = 0; y < ny; y++)
      {
        for (size_t x = 0; x < nx; x++)
        {
          dtRho0Sgz[i] = (dt * duzdznEl) / dtRho0Sgz[i];
          i++;
        } // x
      } // y
    } // z
  } // parallel
}// end of generateInitialDenisty
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute c^2.
 */
void KSpaceFirstOrder3DSolver::computeC2()
{
  if (!mParameters.getC0ScalarFlag())
  {
    float* c2 =  getC2().getData();

    #pragma omp parallel for schedule (static)
    for (size_t i=0; i< getC2().size(); i++)
    {
      c2[i] = c2[i] * c2[i];
    }
  }// matrix
}// computeC2
//----------------------------------------------------------------------------------------------------------------------

/**
 *  Compute part of the new velocity term - gradient of pressure.
 * <b>Matlab code:</b> \n
 *
 *\verbatim
    bsxfun(\@times, ddx_k_shift_pos, kappa .* fftn(p))
  \endverbatim
 */
void KSpaceFirstOrder3DSolver::computePressureGradient()
{
  // Compute FFT of pressure
  getTempFftwX().computeR2CFft3D(getP());

  #pragma omp parallel
  {
    float* pkX = getTempFftwX().getData();
    float* pkY = getTempFftwY().getData();
    float* pkZ = getTempFftwZ().getData();

    const float* kappa        = getKappa().getData();
    const float* ddxKShiftPos = getDdxKShiftPos().getData();
    const float* ddyKShiftPos = getDdyKShiftPos().getData();
    const float* ddzKShiftPos = getDdzKShiftPos().getData();

    const size_t nz = getTempFftwX().getDimensionSizes().nz;
    const size_t ny = getTempFftwX().getDimensionSizes().ny;
    const size_t nx = getTempFftwX().getDimensionSizes().nx;

    const size_t slabSize = (ny * nx) << 1;

    #pragma omp for schedule (static)
    for (size_t z = 0; z < nz; z++)
    {
      register size_t i = z * slabSize;
      const size_t z2 = z<<1;
      for (size_t y = 0; y < ny; y++)
      {
        const size_t y2 = y<<1;
        for (size_t x = 0; x < nx;  x++)
        {
          // kappa ./ p_k
          const float eKappa  = kappa[i>>1];
          const float ePkXRe = pkX[i]   * eKappa;
          const float ePkXIm = pkX[i+1] * eKappa;
          const size_t x2 = x<<1;

          //bsxfun(ddx...)
          pkX[i]   = ePkXRe * ddxKShiftPos[x2]   - ePkXIm * ddxKShiftPos[x2+1];
          pkX[i+1] = ePkXRe * ddxKShiftPos[x2+1] + ePkXIm * ddxKShiftPos[x2];

          //bsxfun(ddy...)
          pkY[i]   = ePkXRe * ddyKShiftPos[y2]   - ePkXIm * ddyKShiftPos[y2+1];
          pkY[i+1] = ePkXRe * ddyKShiftPos[y2+1] + ePkXIm * ddyKShiftPos[y2];

          //bsxfun(ddz...)
          pkZ[i]   = ePkXRe * ddzKShiftPos[z2]   - ePkXIm * ddzKShiftPos[z2+1];
          pkZ[i+1] = ePkXRe * ddzKShiftPos[z2+1] + ePkXIm * ddzKShiftPos[z2];

          i +=2;
        } // x
      } // y
    } // z
  }// parallel
}// end of computePressureGradient
//----------------------------------------------------------------------------------------------------------------------

/**
 * Calculate three temporary sums in the new pressure formula non-linear absorbing case, SSE2 version.
 */
void KSpaceFirstOrder3DSolver::computePressureTermsNonlinearSSE2(RealMatrix& densitySum,
                                                                 RealMatrix& nonlinearTerm,
                                                                 RealMatrix& velocityGradientSum)
{
  // step of 4
  const size_t nElementsSSE = (densitySum.size() >> 2) << 2;

  const float* rhoX = getRhoX().getData();
  const float* rhoY = getRhoY().getData();
  const float* rhoZ = getRhoZ().getData();

  const float* duxdx = getDuxdx().getData();
  const float* duydy = getDuydy().getData();
  const float* duzdz = getDuzdz().getData();

  float  bOnAScalaValue  = mParameters.getBOnAScalar();
  float  rho0ScalarValue = mParameters.getRho0Scalar();
  // set BonA to be either scalar or a matrix
        float* bOnA;
        size_t bOnAShift;
  const bool   bOnAFlag = mParameters.getBOnAScalarFlag();

  if (bOnAFlag)
  {
    bOnA = &bOnAScalaValue;
    bOnAShift = 0;
  }
  else
  {
    bOnA = getBOnA().getData();
    bOnAShift = 1;
  }


  // set rho0 to be either scalar or a matrix
        float* rho0Data;
        size_t rho0Shift;
  const bool   rho0Flag = mParameters.getRho0ScalarFlag();


  if (rho0Flag)
  {
    rho0Data = &rho0ScalarValue;
    rho0Shift = 0;
  }
  else
  {
    rho0Data = getRho0().getData();
    rho0Shift = 1;
  }

  // compute loop
  #pragma  omp parallel
  {
    float* eDensitySum          = densitySum.getData();
    float* eNonlinearTerm       = nonlinearTerm.getData();
    float* eVelocityGradientSum = velocityGradientSum.getData();


    const __m128 sseTwo   = _mm_set1_ps(2.0f);
          __m128 sseBOnA  = _mm_set1_ps(mParameters.getBOnAScalar());
          __m128 sseRho0  = _mm_set1_ps(mParameters.getRho0Scalar());


   #pragma omp for schedule (static) nowait
   for (size_t i = 0; i < nElementsSSE; i +=4 )
   {
      if (!bOnAFlag) sseBOnA = _mm_load_ps(&bOnA[i]);

      __m128 xmm1 = _mm_load_ps(&rhoX[i]);
      __m128 xmm2 = _mm_load_ps(&rhoY[i]);
      __m128 xmm3 = _mm_load_ps(&rhoZ[i]);

      if (!rho0Flag)  sseRho0 = _mm_load_ps(&rho0Data[i]);

      __m128 sseRhoSumSq;
      __m128 sseRhoSum;

      //  register const float rho_xyz_el = rhox_data[i] + rhoy_data[i] + rhoz_data[i];
      sseRhoSum = _mm_add_ps(xmm1, xmm2);
      sseRhoSum = _mm_add_ps(xmm3, sseRhoSum);

      // RHO_Temp_Data[i]  = rho_xyz_el;
      _mm_stream_ps(&eDensitySum[i], sseRhoSum);

      //  BonA_Temp_Data[i] =  ((BonA * (rho_xyz_el * rho_xyz_el)) / (2.0f * rho0_data[i])) + rho_xyz_el;
      sseRhoSumSq = _mm_mul_ps(sseRhoSum, sseRhoSum);// (rho_xyz_el * rho_xyz_el)

      xmm1           = _mm_mul_ps(sseRhoSumSq, sseBOnA);      //((BonA * (rho_xyz_el * rho_xyz_el))
      xmm2           = _mm_mul_ps(sseTwo, sseRho0);             // (2.0f * rho0_data[i])
      xmm3           = _mm_div_ps(xmm1, xmm2);                    // (BonA * (rho_xyz_el * rho_xyz_el)) /  (2.0f * rho0_data[i])

      xmm1           = _mm_add_ps(xmm3, sseRhoSum);          // + rho_xyz_el

      _mm_stream_ps(&eNonlinearTerm[i], xmm1);   //bypass cache

      xmm1       = _mm_load_ps(&duxdx[i]); //dudx
      xmm2       = _mm_load_ps(&duydy[i]); //dudu
      xmm3       = _mm_load_ps(&duzdz[i]); //dudz

      __m128 xmmAcc = _mm_add_ps(xmm1, xmm2);
      xmmAcc        = _mm_add_ps(xmmAcc, xmm3);
      xmmAcc        = _mm_mul_ps(xmmAcc, sseRho0);

      _mm_stream_ps(&eVelocityGradientSum[i],xmmAcc);

    // BonA_Temp_Data[i] =  ((BonA * (rho_xyz_el * rho_xyz_el)) / (2.0f * rho0_data[i])) + rho_xyz_el;
    }

    // non SSE code, in OpenMP only the last thread does this
    #ifdef _OPENMP
      if (omp_get_thread_num() == omp_get_num_threads() -1)
    #endif
    {
      for (size_t i = nElementsSSE; i < densitySum.size() ; i++)
      {
        register const float rhoSum = rhoX[i] + rhoY[i] + rhoZ[i];

        eDensitySum[i]          = rhoSum;
        eNonlinearTerm[i]       = ((bOnA[i * bOnAShift] * (rhoSum * rhoSum)) /
                                   (2.0f * rho0Data[i* rho0Shift])) + rhoSum;
        eVelocityGradientSum[i] = rho0Data[i * rho0Shift] * (duxdx[i] + duydy[i] + duzdz[i]);
      }
    }
  }// parallel
 } // end of computePressureTermsNonlinearSSE2
//----------------------------------------------------------------------------------------------------------------------


 /**
  * Calculate two temporary sums in the new pressure formula, linear absorbing case.
  */
void KSpaceFirstOrder3DSolver::computePressureTermsLinear(RealMatrix& densitySum,
                                                          RealMatrix& velocityGradientSum)
 {
  const size_t size = mParameters.getFullDimensionSizes().nElements();

  #pragma  omp parallel
  {
    const float* rhoX = getRhoX().getData();
    const float* rhoY = getRhoY().getData();
    const float* rhoZ = getRhoZ().getData();

    const float* duxdx = getDuxdx().getData();
    const float* duydy = getDuydy().getData();
    const float* duzdz = getDuzdz().getData();

    const float* rho0 = nullptr;

    const float eRho0 = mParameters.getRho0Scalar();
    if (!mParameters.getRho0ScalarFlag())
    {
      rho0 = getRho0().getData();
    }

    float* eDensitySum  = densitySum.getData();
    float* eVelocityGradientSum = velocityGradientSum.getData();

    #pragma omp for schedule (static)
    for (size_t i = 0; i < size; i++)
    {
      eDensitySum[i] = rhoX[i] + rhoY[i] + rhoZ[i];
    }

    if (mParameters.getRho0ScalarFlag())
    { // scalar
      #pragma omp for schedule (static)
      for (size_t i = 0; i < size; i++)
      {
        eVelocityGradientSum[i] = eRho0 * (duxdx[i] + duydy[i] + duzdz[i]);
      }
    }
    else
    { // matrix
      #pragma omp for schedule (static)
      for (size_t i = 0; i < size; i++)
      {
        eVelocityGradientSum[i] = rho0[i] * (duxdx[i] + duydy[i] + duzdz[i]);
      }
    }
  } // parallel
}// end of computePressureTermsLinear
//----------------------------------------------------------------------------------------------------------------------


 /**
  * Compute absorbing term with abosrbNabla1 and absorbNabla2, SSE2 version
  */
void KSpaceFirstOrder3DSolver::computeAbsorbtionTermSSE2(FftwComplexMatrix& fftPart1,
                                                         FftwComplexMatrix& fftPart2)
{
  const float* absorbNabla1 = getAbsorbNabla1().getData();
  const float* absorbNabla2 = getAbsorbNabla2().getData();

  const size_t nElements     = fftPart1.size();
  const size_t nElementsSSE = (fftPart1.size() >> 1) << 1;

  #pragma omp parallel
  {
    float * fft1Data  = fftPart1.getData();
    float * fft2Data  = fftPart2.getData();

    #pragma omp for schedule (static) nowait
    for (size_t i = 0; i < nElementsSSE; i += 2)
    {
       __m128 sseAbsorbNabla1 = _mm_set_ps(absorbNabla1[i + 1], absorbNabla1[i + 1], absorbNabla1[i], absorbNabla1[i]);
       __m128 sseFft1  = _mm_load_ps(&fft1Data[2 * i]);

              sseFft1  = _mm_mul_ps(sseAbsorbNabla1, sseFft1);
                         _mm_store_ps(&fft1Data[2 * i],    sseFft1);
    }

    #pragma omp for schedule (static)
    for (size_t i = 0; i < nElements; i += 2)
    {
      __m128 sseAbsorbNabla2 = _mm_set_ps(absorbNabla2[i + 1 ], absorbNabla2[i + 1], absorbNabla2[i], absorbNabla2[i]);
      __m128 sseFft2  = _mm_load_ps(&fft2Data[2 * i]);

             sseFft2  = _mm_mul_ps(sseAbsorbNabla2, sseFft2);
                        _mm_store_ps(&fft2Data[2 * i], sseFft2);
      }

    //-- non SSE code --//
    #ifdef _OPENMP
      if (omp_get_thread_num() == omp_get_num_threads() -1)
    #endif
    {
      for (size_t i = nElementsSSE; i < nElements ; i++)
      {
        fft1Data[(i << 1) ]    *= absorbNabla1[i];
        fft1Data[(i << 1) + 1] *= absorbNabla1[i];

        fft2Data[(i << 1) ]    *=  absorbNabla2[i];
        fft2Data[(i << 1) + 1] *=  absorbNabla2[i];
      }
    }

  }// parallel
 } // end of computeAbsorbtionTermSSE2
//----------------------------------------------------------------------------------------------------------------------


 /**
  * @brief Sum sub-terms to calculate new pressure, after FFTs, non-linear case.
  */
void KSpaceFirstOrder3DSolver::sumPressureTermsNonlinear(const RealMatrix& absorbTauTerm,
                                                         const RealMatrix& absorbEtaTerm,
                                                         const RealMatrix& nonlinearTerm)
{
  const float* eAbsorbTauTerm = absorbTauTerm.getData();
  const float* eAbsorbEtaTerm = absorbEtaTerm.getData();

  const size_t nElements = mParameters.getFullDimensionSizes().nElements();
  const float  divider = 1.0f / static_cast<float>(nElements);

  const bool   c0ScalarFlag = mParameters.getC0ScalarFlag();
  const float  c2Scalar     = (c0ScalarFlag) ? mParameters.getC2Scalar() : 0;
  const float* c2Matrix     = (c0ScalarFlag) ? nullptr : getC2().getData();

  const bool   areTauAndEtaScalars = mParameters.getC0ScalarFlag() && mParameters.getAlphaCoeffScalarFlag();
  const float  absorbTauScalar     = (areTauAndEtaScalars) ? mParameters.getAbsorbTauScalar() : 0;
  const float* absorbTauMatrix     = (areTauAndEtaScalars) ? nullptr : getAbsorbTau().getData();

  const float  absorbEtaScalar     = (areTauAndEtaScalars) ? mParameters.getAbsorbEtaScalar() : 0;
  const float* absorbEtaMatrix     = (areTauAndEtaScalars) ? nullptr : getAbsorbEta().getData();;

  #pragma omp parallel
  {
    const float* bOnA = nonlinearTerm.getData();
    float*       p    = getP().getData();

    #pragma omp for schedule (static)
    for (size_t i = 0; i < nElements; i++)
    {
      const float c2        = (c0ScalarFlag) ?        c2Scalar        : c2Matrix[i];
      const float absorbTau = (areTauAndEtaScalars) ? absorbTauScalar : absorbTauMatrix[i];
      const float absorbEta = (areTauAndEtaScalars) ? absorbEtaScalar : absorbEtaMatrix[i];

      p[i] = c2 *(bOnA[i] + (divider * ((eAbsorbTauTerm[i] * absorbTau) - (eAbsorbEtaTerm[i] * absorbEta))));
    }
  }// parallel
}// end of sumPressureTermsNonlinear
//----------------------------------------------------------------------------------------------------------------------

 /**
  * Sum sub-terms to calculate new pressure, after FFTs, linear case.
  */
void KSpaceFirstOrder3DSolver::sumPressureTermsLinear(const RealMatrix& absorbTauTerm,
                                                      const RealMatrix& absorbEtaTerm,
                                                      const RealMatrix& densitySum)
{
  const float* eAbsorbTauTerm = absorbTauTerm.getData();
  const float* eAbsorbEtaTerm = absorbEtaTerm.getData();

  const size_t nElements = mParameters.getFullDimensionSizes().nElements();
  const float  divider = 1.0f / static_cast<float>(nElements);

  const bool   c0ScalarFlag = mParameters.getC0ScalarFlag();
  const float  c2Scalar     = (c0ScalarFlag) ? mParameters.getC2Scalar() : 0;
  const float* c2Matrix     = (c0ScalarFlag) ? nullptr : getC2().getData();

  const bool   areTauAndEtaScalars = mParameters.getC0ScalarFlag() && mParameters.getAlphaCoeffScalarFlag();
  const float  absorbTauScalar     = (areTauAndEtaScalars) ? mParameters.getAbsorbTauScalar() : 0;
  const float* absorbTauMatrix     = (areTauAndEtaScalars) ? nullptr : getAbsorbTau().getData();

  const float  absorbEtaScalar     = (areTauAndEtaScalars) ? mParameters.getAbsorbEtaScalar() : 0;
  const float* absorbEtaMatrix     = (areTauAndEtaScalars) ? nullptr : getAbsorbEta().getData();;

  #pragma omp parallel
  {
    const float* eDenistySum = densitySum.getData();
          float* p           = getP().getData();

    #pragma omp for schedule (static)
    for (size_t i = 0; i < nElements; i++)
    {
      const float c2        = (c0ScalarFlag) ?        c2Scalar        : c2Matrix[i];
      const float absorbTau = (areTauAndEtaScalars) ? absorbTauScalar : absorbTauMatrix[i];
      const float absorbEta = (areTauAndEtaScalars) ? absorbEtaScalar : absorbEtaMatrix[i];

      p[i] = c2 * (eDenistySum[i] + (divider * ((eAbsorbTauTerm[i] * absorbTau) - (eAbsorbEtaTerm[i] * absorbEta))));
    }
  }// parallel
}// end of sumPressureTermsLinear
//----------------------------------------------------------------------------------------------------------------------

/**
 * Sum sub-terms for new p, nonlinear lossless case.
 */
 void KSpaceFirstOrder3DSolver::sumPressureTermsNonlinearLossless()
 {
  #pragma omp parallel
  {
    const size_t nElements = mParameters.getFullDimensionSizes().nElements();
    float* p = getP().getData();

    const float* rhoX = getRhoX().getData();
    const float* rhoY = getRhoY().getData();
    const float* rhoZ = getRhoZ().getData();

    const bool   c0ScalarFlag = mParameters.getC0ScalarFlag();
    const float  c2Scalar     = (c0ScalarFlag) ? mParameters.getC2Scalar() : 0;
    const float* c2Matrix     = (c0ScalarFlag) ? nullptr : getC2().getData();

    const bool   nonlinearFlag = mParameters.getBOnAScalarFlag();
    const float  bOnAScalar    = (nonlinearFlag) ? mParameters.getBOnAScalar(): 0;
    const float* bOnAMatrix    = (nonlinearFlag) ? nullptr : getBOnA().getData();

    const bool   rho0ScalarFlag = mParameters.getRho0ScalarFlag();
    const float  rho0Scalar     = (rho0ScalarFlag) ? mParameters.getRho0Scalar() : 0;
    const float* rho0Matrix     = (rho0ScalarFlag) ? nullptr : getRho0().getData();

    #pragma omp for schedule (static)
    for (size_t i = 0; i < nElements; i++)
    {
      const float c2   = (c0ScalarFlag)   ? c2Scalar   : c2Matrix[i];
      const float bOnA = (nonlinearFlag)  ? bOnAScalar : bOnAMatrix[i];
      const float rho0 = (rho0ScalarFlag) ? rho0Scalar : rho0Matrix[i];

      const float sumDensity = rhoX[i] + rhoY[i] + rhoZ[i];

      p[i] = c2 * (sumDensity + (bOnA * (sumDensity * sumDensity) / (2.0f * rho0)));
    }
  }// parallel
 }// end of sumPressureTermsNonlinearLossless
//----------------------------------------------------------------------------------------------------------------------

 /**
  * Sum sub-terms for new pressure, linear lossless case.
  */
 void KSpaceFirstOrder3DSolver::sumPressureTermsLinearLossless()
{
  #pragma omp parallel
  {
    const float* rhoX = getRhoX().getData();
    const float* rhoY = getRhoY().getData();
    const float* rhoZ = getRhoZ().getData();
          float* p    = getP().getData();

    const size_t nElements = mParameters.getFullDimensionSizes().nElements();

    if (mParameters.getC0ScalarFlag())
    {
      const float c2 = mParameters.getC2Scalar();

      #pragma omp for schedule (static)
      for (size_t i = 0; i < nElements; i++)
      {
        p[i] = c2 * (rhoX[i] + rhoY[i] + rhoZ[i]);
      }
    }
    else
    {
      const float* c2 = getC2().getData();

      #pragma omp for schedule (static)
      for (size_t i = 0; i < nElements; i++)
      {
        p[i] = c2[i] * (rhoX[i] + rhoY[i] + rhoZ[i]);
      }
    }
  }// parallel
}// end of sumPressureTermsLinearLossless()
//----------------------------------------------------------------------------------------------------------------------

/**
 * Calculated shifted velocities.
 *
 * <b>Matlab code:</b> \n
 *
 * \verbatim
    ux_shifted = real(ifft(bsxfun(\@times, x_shift_neg, fft(ux_sgx, [], 1)), [], 1));
    uy_shifted = real(ifft(bsxfun(\@times, y_shift_neg, fft(uy_sgy, [], 2)), [], 2));
    uz_shifted = real(ifft(bsxfun(\@times, z_shift_neg, fft(uz_sgz, [], 3)), [], 3));
  \endverbatim
 */

void KSpaceFirstOrder3DSolver::computeShiftedVelocity()
{
  const FloatComplex* xShiftNegR  = reinterpret_cast<FloatComplex*>(getXShiftNegR().getData());
  const FloatComplex* yShiftNegR  = reinterpret_cast<FloatComplex*>(getYShiftNegR().getData());
  const FloatComplex* zShiftNegR  = reinterpret_cast<FloatComplex*>(getZShiftNegR().getData());

        FloatComplex* tempFftShift = reinterpret_cast<FloatComplex*>(getTempFftwShift().getData());


  // sizes of frequency spaces
  DimensionSizes xShiftDims    = mParameters.getFullDimensionSizes();
                 xShiftDims.nx = xShiftDims.nx / 2 + 1;

  DimensionSizes yShiftDims    = mParameters.getFullDimensionSizes();
                 yShiftDims.ny = yShiftDims.ny / 2 + 1;

  DimensionSizes zShiftDims    = mParameters.getFullDimensionSizes();
                 zShiftDims.nz = zShiftDims.nz / 2 + 1;

  // normalization constants for FFTs
  const float dividerX = 1.0f / static_cast<float>(mParameters.getFullDimensionSizes().nx);
  const float dividerY = 1.0f / static_cast<float>(mParameters.getFullDimensionSizes().ny);
  const float dividerZ = 1.0f / static_cast<float>(mParameters.getFullDimensionSizes().nz);

  //-------------------------------------------------- ux_shifted ----------------------------------------------------//
  getTempFftwShift().computeR2CFft1DX(getUxSgx());

  #pragma omp parallel for schedule (static)
  for (size_t z = 0; z < xShiftDims.nz; z++)
  {
    register size_t i = z *  xShiftDims.ny * xShiftDims.nx;
    for (size_t y = 0; y < xShiftDims.ny; y++)
    {
      for (size_t x = 0; x < xShiftDims.nx; x++)
      {
        FloatComplex temp;

        temp.real = ((tempFftShift[i].real * xShiftNegR[x].real) -
                     (tempFftShift[i].imag * xShiftNegR[x].imag)
                    ) * dividerX;


        temp.imag = ((tempFftShift[i].imag * xShiftNegR[x].real) +
                     (tempFftShift[i].real * xShiftNegR[x].imag)
                    ) * dividerX;

        tempFftShift[i] = temp;

        i++;
      } // x
    } // y
  }//z*/
  getTempFftwShift().computeC2RFft1DX(getUxShifted());


  //-------------------------------------------------- uy shifted ----------------------------------------------------//
  getTempFftwShift().computeR2CFft1DY(getUySgy());

  #pragma omp parallel for schedule (static)
  for (size_t z = 0; z < yShiftDims.nz; z++)
  {
    register size_t i = z *  yShiftDims.ny * yShiftDims.nx;
    for (size_t y = 0; y < yShiftDims.ny; y++)
    {
      for (size_t x = 0; x < yShiftDims.nx; x++)
      {
        FloatComplex temp;

        temp.real = ((tempFftShift[i].real * yShiftNegR[y].real) -
                     (tempFftShift[i].imag * yShiftNegR[y].imag)) *
                      dividerY;


        temp.imag = ((tempFftShift[i].imag * yShiftNegR[y].real) +
                     (tempFftShift[i].real * yShiftNegR[y].imag)
                    ) * dividerY;

        tempFftShift[i] = temp;

        i++;
      } // x
    } // y
  }//z
  getTempFftwShift().computeC2RFft1DY(getUyShifted());


  //-------------------------------------------------- uz_shifted ----------------------------------------------------//
  getTempFftwShift().computeR2CFft1DZ(getUzSgz());
  #pragma omp parallel for schedule (static)
  for (size_t z = 0; z < zShiftDims.nz; z++)
  {
    register size_t i = z *  zShiftDims.ny * zShiftDims.nx;
    for (size_t y = 0; y < zShiftDims.ny; y++)
    {
      for (size_t x = 0; x < zShiftDims.nx; x++)
      {
        FloatComplex temp;

        temp.real = ((tempFftShift[i].real * zShiftNegR[z].real) -
                     (tempFftShift[i].imag * zShiftNegR[z].imag)) *
                      dividerZ;


        temp.imag = ((tempFftShift[i].imag * zShiftNegR[z].real) +
                     (tempFftShift[i].real * zShiftNegR[z].imag)
                    ) * dividerZ;

        tempFftShift[i] = temp;

        i++;
      } // x
    } // y
  }//z
  getTempFftwShift().computeC2RFft1DZ(getUzShifted());
}// end of computeShiftedVelocity
//----------------------------------------------------------------------------------------------------------------------

/**
 * Print progress statistics.
 */
void KSpaceFirstOrder3DSolver::printStatistics()
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
 * Is time to checkpoint?
 */
bool KSpaceFirstOrder3DSolver::isTimeToCheckpoint()
{
  if (!mParameters.isCheckpointEnabled()) return false;

  mTotalTime.stop();

  return (mTotalTime.getElapsedTime() > static_cast<float>(mParameters.getCheckpointInterval()));

}// end of isTimeToCheckpoint
//----------------------------------------------------------------------------------------------------------------------

/**
 * Was the loop interrupted to checkpoint?
 */
bool KSpaceFirstOrder3DSolver::isCheckpointInterruption() const
{
  return (mParameters.getTimeIndex() != mParameters.getNt());
}// end of isCheckpointInterruption
//----------------------------------------------------------------------------------------------------------------------

/**
 * Check the output file has the correct format and version.
 */
void KSpaceFirstOrder3DSolver::checkOutputFile()
{
  // The header has already been read
  Hdf5FileHeader& fileHeader = mParameters.getFileHeader();
  Hdf5File&        outputFile = mParameters.getOutputFile();

  // test file type
  if (fileHeader.getFileType() != Hdf5FileHeader::FileType::kOutput)
  {
    throw ios::failure(Logger::formatMessage(kErrFmtBadOutputFileFormat, mParameters.getOutputFileName().c_str()));
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
  DimensionSizes outputDimSizes;
  outputFile.readScalarValue(outputFile.getRootGroup(),
                             kNxName,
                             outputDimSizes.nx);

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
void KSpaceFirstOrder3DSolver::checkCheckpointFile()
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
void KSpaceFirstOrder3DSolver::loadElapsedTimeFromOutputFile()
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

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Private methods --------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

