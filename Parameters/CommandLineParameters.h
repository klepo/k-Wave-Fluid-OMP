/**
 * @file      CommandLineParameters.h
 *
 * @author    Jiri Jaros \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            jarosjir@fit.vutbr.cz
 *
 * @brief     The header file containing the command line parameters.
 *
 * @version   kspaceFirstOrder 2.17
 *
 * @date      29 August    2012, 11:25 (created) \n
 *            14 March     2019, 11:02 (revised)
 *
 * @section   Params Command Line Parameters
 *
 * The C++ code requires two mandatory parameters and accepts a few optional parameters and  flags. Ill parameters,
 * bad simulation files, and runtime errors such as out-of-memory problems, lead to an exception followed by an error
 * message shown and execution termination.
 *
 * The mandatory parameters <tt>-i</tt> and <tt>-o</tt> specify the input and output file. The file names respect the
 * path conventions for particular operating system. If any of the files is not specified, cannot be found or created,
 * an error message is shown and the code terminates.
 *
 * The <tt>-t</tt> parameter sets the number of threads used, which defaults the system maximum. If the system support
 * hyperthreading, it is recommended to use only a half of the threads to prevent cache overloading. If possible enable
 * tread binding and placement using export OMP_PROC_BIND=true.
 *
 *
 * The <tt>-r</tt> parameter specifies how often information about the simulation progress is printed out to the command
 * line. By default, the C++ code prints out the  progress of the simulation, the elapsed time, and the estimated
 * time of completion in intervals corresponding to 5% of the total number of times steps.
 *
 * The <tt>-c</tt> parameter specifies the compression level used by the ZIP library to reduce the size of the output
 * file. The actual compression rate is highly dependent on the shape of the sensor mask and the range of stored
 * quantities and may be computationally expensive. In general, the output data is very hard to compress, and using
 * higher compression levels can greatly increase the time to save data while not having a large impact on the final
 * file size. That's why we decided to disable compression in default settings.
 *
 * The <tt>\--benchmark</tt> parameter enables the total length of simulation (i.e., the number of time steps) to be
 * overridden by setting a new number of time  steps to simulate. This is particularly useful for performance evaluation
 * and benchmarking. As the code performance is relatively stable, 50-100 time steps is usually enough to predict the
 * simulation duration. This parameter can also be used to quickly check the simulation is set up correctly.
 *
 * The <tt>\--verbose</tt> parameter enables to select between three levels of verbosity. For  routine simulations, the
 * verbose level of 0 (the default one) is usually sufficient. For more information about the simulation, checking the
 * parameters of the simulation, code version, GPU used, file paths, and debugging running scripts, verbose levels
 * 1 and 2 may be very useful.
 *
 * The <tt>-h</tt> and <tt>\--help</tt> parameters print all the parameters of the C++ code. The <tt>\--version </tt>
 * parameter reports detail information about the code useful for  debugging and bug reports. It prints out the internal
 * version, the build date and time, the git hash allowing us to track the version of the source code, the operating
 * system, the compiler name and version and the instruction set used.
 *
 * For jobs that are expected to run for a very long time, it may be useful to  checkpoint and restart the execution.
 * One motivation is the wall clock limit  per task on clusters where jobs must fit within a given time span (e.g. 24
 * hours). The second motivation is a level of fault-tolerance, where you can back up the state of the simulation after
 * a predefined period. To enable checkpoint-restart, the user is asked to specify a file to store the actual state of
 * the simulation by  <tt>\--checkpoint_file</tt> and the period in seconds after which the simulation will be
 * interrupted by <tt>\--checkpoint_interval</tt>.  When running on a cluster, please allocate enough time for the
 * checkpoint procedure  that can take a non-negligible amount of time (7 matrices have to be stored in  the
 * checkpoint file and all aggregated quantities are flushed into the output file).
 * Alternatively, the user can specify the number of time steps by <tt>\--checkpoint_timesteps</tt> after which the
 * simulation is interrupted. The time step interval is calculated from the beginning of current leg, not from the
 * beginning of the whole simulation. The user can combine both approaches, seconds and time steps. In this case the
 * first condition met triggers the checkpoint.
 * Please note, that the checkpoint
 * file name and path is not checked at the beginning of the simulation, but at the time the code starts
 * checkpointing. Thus make sure the file path was correctly specified (otherwise you will not find out the simulation
 * crashed until the first leg of the simulation finishes). The rationale behind this is that to keep as high level of
 * fault tolerance as possible, the checkpoint file should be touched even when really necessary.
 * When controlling a multi-leg simulation by a script loop, the parameters of the code remains the same in all legs.
 * The first leg of the simulation creates a checkpoint  file while the last one deletes it. If the checkpoint file is
 * not found the simulation starts from the beginning. In order to find out how many steps have been finished, please
 * open the output file and read the variable <tt>t_index</tt> and compare it with <tt>Nt</tt> (e.g. by the h5dump
 * command).
 *
 *
 * The remaining flags specify the output quantities to be recorded during the  simulation and stored on disk analogous
 * to  the sensor.record input. If the <tt>-p</tt> or <tt>\--p\_raw</tt> flags are set (these are equivalent), a time
 * series of  the acoustic pressure at the grid points specified by  the sensor mask is recorded. If the
 * <tt>\--p_rms</tt>, <tt>\--p_max</tt>, <tt>\--p_min</tt> flags  are set, the root mean square and/or maximum and/or
 * minimum values of the pressure at the grid points specified by  the sensor mask are recorded. If the
 * <tt>\--p_final</tt> flag is set, the values for the entire acoustic pressure field in the final time step of the
 * simulation is stored (this will always include the PML, regardless of  the setting for <tt> 'PMLInside'</tt>).
 * The flags <tt>\--p_max_all</tt> and <tt>\--p_min_all</tt> allow to calculate the maximum and  minimum values over the
 * entire acoustic pressure field, regardless on the shape of the sensor mask. Flags to record the acoustic particle
 * velocity are defined in an analogous fashion. For proper calculation of acoustic intensity, the particle velocity
 * has to be shifted onto the same grid as the acoustic  pressure. This can be done by setting
 * <tt>\--u_non_staggered_raw</tt> flag, that first shifts the  particle velocity and then samples the grid points
 * specified by the sensor mask. Since the  shift operation requires additional FFTs, the impact on the simulation time
 * may be significant.
 *
 * Any combination of <tt>p</tt> and <tt>u</tt> flags is admissible. If no output flag is set, a time-series for the
 * acoustic pressure is recorded. If it is not necessary to collect the output quantities over the entire simulation,
 * the starting time step when the collection begins can be specified using the -s parameter.  Note, the index for the
 * first time step is 1 (this follows the MATLAB indexing convention).
 *
 * The <tt>\--copy_sensor_mask</tt> flag will copy the sensor from the input file to the output  one at the end of the
 * simulation. This helps in post-processing and visualisation of the outputs.
 *
 *
 *
\verbatim
┌───────────────────────────────────────────────────────────────┐
│                  kspaceFirstOrder3D-OMP v1.3                  │
├───────────────────────────────────────────────────────────────┤
│                             Usage                             │
├───────────────────────────────────────────────────────────────┤
│                     Mandatory parameters                      │
├───────────────────────────────────────────────────────────────┤
│ -i <file_name>                │ HDF5 input file               │
│ -o <file_name>                │ HDF5 output file              │
├───────────────────────────────┴───────────────────────────────┤
│                      Optional parameters                      │
├───────────────────────────────┬───────────────────────────────┤
│ -t <num_threads>              │ Number of CPU threads         │
│                               │  (default =  4)               │
│ -g <device_number>            │ GPU device to run on          │
│                               │   (default = the first free)  │
│ -r <interval_in_%>            │ Progress print interval       │
│                               │   (default =  5%)             │
│ -c <compression_level>        │ Compression level <0,9>       │
│                               │   (default = 0)               │
│ --benchmark <time_steps>      │ Run only a specified number   │
│                               │   of time steps               │
│ --verbose <level>             │ Level of verbosity <0,2>      │
│                               │   0 - basic, 1 - advanced,    │
│                               │   2 - full                    │
│                               │   (default = basic)           │
│ -h, --help                    │ Print help                    │
│ --version                     │ Print version and build info  │
├───────────────────────────────┼───────────────────────────────┤
│ --checkpoint_file <file_name> │ HDF5 Checkpoint file          │
│ --checkpoint_interval <sec>   │ Checkpoint after a given      │
│                               │   number of seconds           │
│ --checkpoint_timesteps <step> │ Checkpoint after a given      │
│                               │   number of time steps        │
├───────────────────────────────┴───────────────────────────────┤
│                          Output flags                         │
├───────────────────────────────┬───────────────────────────────┤
│ -p                            │ Store acoustic pressure       │
│                               │   (default output flag)       │
│                               │   (the same as --p_raw)       │
│ --p_raw                       │ Store raw time series of p    │
│ --p_c                         │ Store compressed              │
│                               │    time series of p           │
│ --p_rms                       │ Store rms of p                │
│ --p_max                       │ Store max of p                │
│ --p_min                       │ Store min of p                │
│ --p_max_all                   │ Store max of p (whole domain) │
│ --p_min_all                   │ Store min of p (whole domain) │
│ --p_final                     │ Store final pressure field    │
├───────────────────────────────┼───────────────────────────────┤
│ -u                            │ Store ux, uy, uz              │
│                               │    (the same as --u_raw)      │
│ --u_raw                       │ Store raw time series of      │
│                               │    ux, uy, uz                 │
│ --u_c                         │ Store compressed              │
│                               │    time series of ux, uy, uz  │
│ --u_non_staggered_raw         │ Store non-staggered raw time  │
│                               │   series of ux, uy, uz        │
│ --u_non_staggered_c           │ Store non-staggered compressed│
│                               │   time series of ux, uy, uz   │
│ --u_rms                       │ Store rms of ux, uy, uz       │
│ --u_max                       │ Store max of ux, uy, uz       │
│ --u_min                       │ Store min of ux, uy, uz       │
│ --u_max_all                   │ Store max of ux, uy, uz       │
│                               │   (whole domain)              │
│ --u_min_all                   │ Store min of ux, uy, uz       │
│                               │   (whole domain)              │
│ --u_final                     │ Store final acoustic velocity │
├───────────────────────────────┼───────────────────────────────┤
│ --I_avg                       │ Store average intensity       │
│ --I_avg_c                     │ Store average intensity       │
│                               │   computed using compression  │
│ --Q_term                      │ Store Q term (volume rate of  │
│                               │   heat deposition)            │
│ --Q_term_c                    │ Store Q term (volume rate of  │
│                               │   heat deposition) computed   │
│                               │   using compression           │
│ --block_size                  │ Maximum block size for        │
│                               │   dataset reading (computing  │
│                               │   average intensity without   │
│                               │   compression)                │
├───────────────────────────────┴───────────────────────────────┤
│                Time series compression flags                  │
├───────────────────────────────┬───────────────────────────────┤
│ --frequency <frequency>       │ Frequency for time series     │
│                               │   compression (needs dt for   │
│                               │   computing period)           │
│ --period <period>             │ Period for time series        │
│                               │   compression                 │
│                               │   (default = computed from    │
│                               │   p_source_input dataset)     │
│ --mos <mos>                   │ Multiple of overlap size  for │
│                               │   compression                 │
│                               │   (default = 1)               │
│ --harmonics <harmonics>       │ Number of hamornics for       │
│                               │   compression (default = 1)   │
├───────────────────────────────┼───────────────────────────────┤
│ -s <time_step>                │ When data collection begins   │
│                               │   (default = 1)               │
│ --copy_sensor_mask            │ Copy sensor mask to the       │
│                               │    output file                │
└───────────────────────────────┴───────────────────────────────┘
\endverbatim
 *
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


#ifndef COMMAND_LINE_PARAMETERS_H
#define COMMAND_LINE_PARAMETERS_H

#include <string>


/**
 * @class   CommandLineParameters
 * @brief   The class to parse and store command line parameters.
 * @details The class to parse and store command line parameters.
 */
class CommandLineParameters
{
  public:
    /// Only Parameters can create this class.
    friend class Parameters;

    /// Copy constructor not allowed.
    CommandLineParameters(const CommandLineParameters&) = delete;

    /// Destructor.
    virtual ~CommandLineParameters() {};

    /// operator= not allowed.
    CommandLineParameters& operator=(const CommandLineParameters&) = delete;

    /**
     * @brief Get input file name.
     * @return Input file name.
     */
    const std::string& getInputFileName()      const { return mInputFileName; };
    /**
     * @brief  Get output file name.
     * @return Output file name.
     */
    const std::string& getOutputFileName()     const { return mOutputFileName; };
    /**
     * @brief  Get Checkpoint file name.
     * @return Checkpoint file name.
     */
    const std::string& getCheckpointFileName() const { return mCheckpointFileName; };

    /**
     * @brief  Get number of threads.
     * @return Number of CPU threads value.
     */
    size_t getNumberOfThreads()         const { return mNumberOfThreads; };

    /**
     * @brief Get progress print interval.
     * @return How often to print progress.
     */
    size_t getProgressPrintInterval()   const { return mProgressPrintInterval; };

    /**
     * @brief  Get compression level.
     * @return Compression level value for output and checkpoint files.
     */
    size_t getCompressionLevel()        const { return mCompressionLevel; };

    /**
     * @brief  Is --benchmark set?
     * @return true if the flag is set.
     */
    bool   isBenchmarkEnabled()         const { return mBenchmarkFlag; };

    /**
     * @brief  Get benchmark time step count.
     * @return Number of time steps used to benchmark the code.
     */
    size_t getBenchmarkTimeStepsCount() const { return mBenchmarkTimeStepCount; };

 /**
     * @brief  Is checkpoint enabled?
     * @return true if checkpointing is enabled.
     */
    bool   isCheckpointEnabled()        const { return ((mCheckpointInterval > 0) || (mCheckpointTimeSteps > 0)); };

    /**
     * @brief  Get checkpoint interval.
     * @return Checkpoint interval in seconds.
     */
    size_t getCheckpointInterval()      const { return mCheckpointInterval; };
    /**
     * @brief  Get checkpoint interval in time steps.
     * @return Checkpoint interval in time steps.
     */
    size_t getCheckpointTimeSteps()     const { return mCheckpointTimeSteps; };

    /**
     * @brief  Is --version set?
     * @return true if the flag is set.
     */
    bool   isPrintVersionOnly()         const { return mPrintVersionFlag; };


    //------------------------------------------------ Output flags --------------------------------------------------//
    /**
     * @brief Is --p_raw set?
     * @return true if the flag is set.
     */
    bool   getStorePressureRawFlag()      const { return mStorePressureRawFlag; };
    /**
     * @brief  Is --p_c set?
     * @return true if the flag is set.
     */
    bool   getStorePressureCFlag()        const { return mStorePressureCFlag; };
    /**
     * @brief  Is --p_rms set?
     * @return true if the flag is set.
     */
    bool   getStorePressureRmsFlag()      const { return mStorePressureRmsFlag; };
    /**
     * @brief  Is --p_max set?
     * @return true if the flag is set.
     */
    bool   getStorePressureMaxFlag()      const { return mStorePressureMaxFlag; };
    /**
     * @brief  Is --p_min set?
     * @return true if the flag is set.
     */
    bool   getStorePressureMinFlag()      const { return mStorePressureMinFlag; };
    /**
     * @brief  Is --p_max_all set?
     * @return true if the flag is set.
     */
    bool   getStorePressureMaxAllFlag()   const { return mStorePressureMaxAllFlag; };
    /**
     * @brief  Is --p_min_all set?
     * @return true if the flag is set.
     */
    bool   getStorePressureMinAllFlag()   const { return mStorePressureMinAllFlag; };
    /**
     * @brief  Is --p_final set?
     * @return true if the flag is set.
     */
    bool   getStorePressureFinalAllFlag() const { return mStorePressureFinalAllFlag; };


    /**
     * @brief  Is --u_raw set?
     * @return true if the flag is set.
     */
    bool   getStoreVelocityRawFlag()             const { return mStoreVelocityRawFlag; };
    /**
     * @brief  Is --u_c set?
     * @return true if the flag is set.
     */
    bool   getStoreVelocityCFlag()               const { return mStoreVelocityCFlag; };
    /**
     * @brief  Is --u_non_staggered_raw set?
     * @return true if the flag is set.
     */
    bool   getStoreVelocityNonStaggeredRawFlag() const { return mStoreVelocityNonStaggeredRawFlag; };
    /**
     * @brief  Is --u_non_staggered_c set?
     * @return true if the flag is set.
     */
    bool   getStoreVelocityNonStaggeredCFlag()   const { return mStoreVelocityNonStaggeredCFlag; };
    /**
     * @brief  Is --u_rms set?
     * @return true if the flag is set.
     */
    bool   getStoreVelocityRmsFlag()             const { return mStoreVelocityRmsFlag; };
    /**
     * @brief  Is --u_max set?
     * @return true if the flag is set.
     */
    bool   getStoreVelocityMaxFlag()             const { return mStoreVelocityMaxFlag; };
    /**
     * @brief  Is --u_min set?
     * @return true if the flag is set.
     */
    bool   getStoreVelocityMinFlag()             const { return mStoreVelocityMinFlag; };
    /**
     * @brief  Is --u_max_all set?
     * @return true if the flag is set.
     */
    bool   getStoreVelocityMaxAllFlag()          const { return mStoreVelocityMaxAllFlag; };
    /**
     * @brief  Is --u_min set?
     * @return true if the flag is set.
     */
    bool  getStoreVelocityMinAllFlag()           const { return mStoreVelocityMinAllFlag; };
    /**
     * @brief  Is --u_final set?
     * @return true if the flag is set.
     */
    bool   getStoreVelocityFinalAllFlag()        const { return mStoreVelocityFinalAllFlag; };
    /**
     * @brief  Is --I_avg set?
     * @return true if the flag is set.
     */
    bool   getStoreIntensityAvgFlag()     const { return mStoreIntensityAvgFlag; };
    /**
     * @brief  Is --I_avg_c set?
     * @return true if the flag is set.
     */
    bool   getStoreIntensityAvgCFlag()    const { return mStoreIntensityAvgCFlag; };
    /**
     * @brief  Is --Q_term set?
     * @return true if the flag is set.
     */
    bool   getStoreQTermFlag()            const { return mStoreQTermFlag; };
    /**
     * @brief  Is --Q_term_c set?
     * @return true if the flag is set.
     */
    bool   getStoreQTermCFlag()           const { return mStoreQTermCFlag; };
    /**
     * @brief  Is --copy_mask set set?
     * @return true if the flag is set.
     */
    bool   getCopySensorMaskFlag()        const { return mCopySensorMaskFlag; };

    /**
     * @brief  Get start time index when sensor data collection begins.
     * @return When to start sampling data.
     */
    size_t getSamplingStartTimeIndex()    const { return mSamplingStartTimeStep; };

    /**
     * @brief  Get compression frequency.
     * @return Compression frequency.
     */
    float getFrequency()                  const { return mFrequency; }

    /**
     * @brief  Get compression period.
     * @return Compression period.
     */
    float getPeriod()                     const { return mPeriod; }

    /**
     * @brief  Get Multiple of overlap size for compression.
     * @return Compression multiple of overlap size.
     */
    size_t getMOS()                       const { return mMOS; }

    /**
     * @brief  Get number of harmonics for compression.
     * @return Number of harmonics for compression.
     */
    size_t getHarmonics()                 const { return mHarmonics; }

    /**
     * @brief  Get maximum block size for dataset reading (computing average intensity).
     * @return Block size for dataset reading.
     */
    size_t getBlockSize()                 const { return mBlockSize; }

    /// Print usage of the code
    void printUsage();
    /// Print setup commandline parameters.
    void printComandlineParamers();

    /**
     * @brief Parse commandline parameters.
     * @param [in, out] argc - number of commandline parameters.
     * @param [in, out] argv - commandline parameters.
     *
     * @throw call exit when error in commandline.
     */
    void parseCommandLine(int argc, char** argv);

  protected:
    /// Default constructor - only friend class can create an instance.
    CommandLineParameters();

  private:
    /// Input file name.
    std::string mInputFileName;
    /// Output file name.
    std::string mOutputFileName;
    /// Checkpoint file name.
    std::string mCheckpointFileName;

    /// Number of CPU threads value.
    size_t mNumberOfThreads;
    /// Progress interval value.
    size_t mProgressPrintInterval;
    /// Compression level value for output and checkpoint files.
    size_t mCompressionLevel;

    /// BenchmarkFlag value.
    bool   mBenchmarkFlag;
    /// Number of time steps used to benchmark the code
    size_t mBenchmarkTimeStepCount;
    /// Checkpoint interval in seconds.
    size_t mCheckpointInterval;
    /// Checkpoint interval in time steps.
    size_t mCheckpointTimeSteps;

    /// Print version of the code and exit.
    bool mPrintVersionFlag;


    /// Store raw time-series of pressure over the sensor mask?
    bool mStorePressureRawFlag;
    /// Store compressed time-series of pressure over the sensor mask?
    bool mStorePressureCFlag;
    /// Store RMS of pressure over the the sensor mask?
    bool mStorePressureRmsFlag;
    /// Store maximum of pressure over the sensor mask?
    bool mStorePressureMaxFlag;
    /// Store minimum of pressure over the sensor mask?
    bool mStorePressureMinFlag;
    /// Store maximum of pressure over the whole domain?
    bool mStorePressureMaxAllFlag;
    /// Store minimum of pressure over the whole domain?
    bool mStorePressureMinAllFlag;
    /// Store pressure in the final time step over the whole domain?
    bool mStorePressureFinalAllFlag;

    /// Store raw time-series of velocity over the sensor mask?
    bool mStoreVelocityRawFlag;
    /// Store compressed time-series of velocity over the sensor mask?
    bool mStoreVelocityCFlag;
    /// Store un staggered raw time-series of velocity over the sensor mask?
    bool mStoreVelocityNonStaggeredRawFlag;
    /// Store compressed un staggered raw time-series of velocity over the sensor mask?
    bool mStoreVelocityNonStaggeredCFlag;
    /// Store RMS of velocity over the the sensor mask?
    bool mStoreVelocityRmsFlag;
    /// Store maximum of velocity over the sensor mask?
    bool mStoreVelocityMaxFlag;
    /// Store minimum of velocity over the sensor mask?
    bool mStoreVelocityMinFlag;
    /// Store maximum of velocity over the whole domain?
    bool mStoreVelocityMaxAllFlag;
    /// Store minimum of velocity over the whole domain?
    bool mStoreVelocityMinAllFlag;
    /// Store velocity in the final time step over the whole domain?
    bool mStoreVelocityFinalAllFlag;

    /// Store average intensity?
    bool mStoreIntensityAvgFlag;
    /// Store average intensity using compression?
    bool mStoreIntensityAvgCFlag;
    /// Store Q term (volume rate of heat deposition)?
    bool mStoreQTermFlag;
    /// Store Q term (volume rate of heat deposition) using compression?
    bool mStoreQTermCFlag;

    /// Copy sensor mask to the output file.
    bool   mCopySensorMaskFlag;
    /// StartTimeStep value.
    size_t mSamplingStartTimeStep;

    /// Period for compression.
    float mPeriod = 0.0f;
    /// Frequency for compression.
    float mFrequency = 0.0f;
    /// Multiple of overlap size for compression.
    size_t mMOS = 1;
    /// Number of harmonics for compression.
    size_t mHarmonics = 1;

    /// Maximum block size for dataset reading (computing average intensity).
    /// Default value is computed according to free RAM memory.
    size_t mBlockSize = 0;

    /// Default compression level.
    static constexpr size_t kDefaultCompressionLevel      = 0;
    /// Default progress print interval.
    static constexpr size_t kDefaultProgressPrintInterval = 5;
};// end of class CommandLineParameters
//----------------------------------------------------------------------------------------------------------------------

#endif	/* COMMAND_LINE_PARAMETERS_H */

