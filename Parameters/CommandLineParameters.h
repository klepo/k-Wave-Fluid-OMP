/**
 * @file        CommandLineParameters.h
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au
 * @brief       The header file containing the command line parameters
 * 
 * @version     kspaceFirstOrder3D 2.13
 * @date        29 August 2012, 11:25 (created) \n        
 *              11 October 2012, 17:05 (revised) 
 * 
 * @section Params Command Line Parameters
 * 
 * The  C++ code requires two mandatory parameters and accepts a few optional parameters and flags. 
 * The mandatory parameters \c -i and \c -o specify the input and output file. The file names respect the path conventions for particular operating system.  
 * If any of the files is not specified, cannot be found or created, an error message is shown.
 * 
 * The \c -t parameter sets the number of threads used, which defaults the system maximum. On CPUs with Intel Hyper-Threading (HT), 
 * performance will generally be better if HT is disabled in the BIOS settings. If HT is switched on, the default will be to create 
 * twice as many threads as there are physical processor cores, which might slow down the code execution. Therefore, if the HT is on,
 *  try specifying the number of threads manually for best performance (e.g., 4 for Intel Quad-Core). We recommend experimenting 
 * with this parameter to find the best configuration. Note, if there are other tasks being executed on the system, 
 * it might be useful to further limit the number of threads to prevent system overload.
 * 
 * The \c -r parameter specifies how often information about the simulation progress is printed to the command line. 
 * By default, the C++ code prints out the progress of the simulation, the elapsed time, and the estimated time of 
 * completion in intervals corresponding to 5% of the total number of times steps.
 * 
 * The \c -c parameter specifies the compression level used by the ZIP library to reduce the 
 * size of the output file. The actual compression rate is highly dependent on the shape of 
 * the sensor mask and the range of stored quantities. In general, the output data is very 
 * hard to compress, and using higher compression levels can greatly increase the time to 
 * save data while not having a large impact on the final file size. The default compression 
 * level of 3 represents a balance between compression ratio and performance that is suitable in most cases.
 * 
 * The \c --benchmark parameter enables the total length of simulation (i.e., the number of time steps) 
 * to be overwritten by setting a new number of time steps to simulate. This is particularly useful 
 * for performance evaluation and benchmarking. As the code performance is relatively stable, 50-100 time steps is 
 * usually enough to predict the simulation duration. This parameter can also be used to quickly find the ideal number of CPU threads to use.
 *
 * The \c -h and \c --help parameters print all the parameters of the C++ code, while the \c --version parameter reports the code version and internal build number.
 * 
 * 
 * The remaining flags specify the output quantities to be recorded during the simulation and stored on disk analogous to 
 * the sensor.record input. If the \c -p or \c --p raw flags are set (these are equivalent), a time series of the acoustic pressure at the 
 * grid points specified by the sensor mask is recorded. If the \c --p rms and/or \c --p max flags are set, the root mean square and/or maximum values 
 * of the pressure at the grid points specified by the sensor mask are recorded. Finally, if the \c --p final flag is set, the values for the entire acoustic pressure 
 * field in the final time step of the simulation is stored (this will always include the PML, regardless of the setting for <tt> `PMLInside' </tt>). Flags to record the acoustic 
 * particular velocity are defined in an analogous fashion.
 * 
 * In addition to returning the acoustic pressure and particle velocity, the acoustic intensity at the grid points specified by 
 * the sensor mask can also be calculated and stored. To account for the staggered grid scheme, approximate values for the particle 
 * velocity at the unstaggered grid points are automatically calculated using linear interpolation before multiplication by the acoustic pressure. 
 * Two means of aggregation are possible: \c -I or \c --I_avg calculates and stores the average acoustic intensity, while \c --I max calculates the maximum acoustic intensity.
 * 
 * Any combination of \c p, \c u and \c I fags is admissible. If no output flag is set, a time-series for
 * the acoustic pressure is recorded. If it is not necessary to collect the output quantities over 
 * the entire simulation, the starting time step when the collection begins can be specified 
 * using the -s parameter. Note, the index for the first time step is 1 (this follows the MATLAB indexing convention).
 *
 *
\verbatim
---------------------------------- Usage ---------------------------------  
Mandatory parameters:
  -i <input_file_name>            : HDF5 input file 
  -o <output_file_name>           : HDF5 output file
 
Optional parameters: 
  -t <num_threads>                : Number of CPU threads (default = MAX)
  -r <interval_in_%>              : Progress print interval (default = 5%)
  -c <comp_level>                 : Output file compression level <0,9>
                                      (default = 3)
  --benchmark <steps>             : Run a specified number of time steps
   
  -h                              : Print help
  --help                          : Print help
  --version                       : Print version
   
Output flags:   
  -p                              : Store acoustic pressure 
                                      (default if nothing else is on)
                                      (the same as --p_raw)
  --p_raw                         : Store raw time series of p (default)
  --p_rms                         : Store rms of p
  --p_max                         : Store max of p
  --p_final                       : Store final pressure field 
   
  -u                              : Store ux, uy, uz
                                      (the same as --u_raw)
  --u_raw                         : Store raw time series of ux, uy, uz
  --u_rms                         : Store rms of ux, uy, uz
  --u_max                         : Store max of ux, uy, uz
  --u_final                       : Store final acoustic velocity
   
  -I                              : Store intensity
                                      (the same as --I_avg) 
  --I_avg                         : Store avg of intensity
  --I_max                         : Store max of intensity
   
  -s <timestep>                   : Time step when data collection begins
                                      (default = 1)
--------------------------------------------------------------------------  
\endverbatim
 * 
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
 



#ifndef TCOMMANDLINESPARAMETERS_H
#define	TCOMMANDLINESPARAMETERS_H

#include <cstdlib>
#include <string>


/**
 * @class TCommandLineParameters
 * @brief The class to parse and store command line parameters
 */
class TCommandLineParameters {
public:

    /// Constructor
    TCommandLineParameters();
    /// Destructor
    virtual ~TCommandLineParameters() {};
    
    /// Get input file name
    std::string GetInputFileName()      const {return InputFileName;};    
    /// Get output file name
    std::string GetOutputFileName()     const {return OutputFileName;};
    
    /// Is --benchmark flag set?
    bool IsBenchmarkFlag()              const {return BenchmarkFlag;};
    /// Is --version flag set
    bool IsVersion()                    const {return PrintVersion; };
    /// Get benchmark time step count
    int  GetBenchmarkTimeStepsCount()   const {return BenchmarkTimeStepsCount;};
    
    /// Get compression level
    int  GetCompressionLevel()          const {return CompressionLevel;};
    /// Get number of threads 
    int  GetNumberOfThreads()           const {return NumberOfThreads;};
    /// Get verbose interval
    int  GetVerboseInterval()           const {return VerboseInterval;};
    /// Get start time index when sensor data collection begins
    int GetStartTimeIndex()             const {return StartTimeStep;};
   
    /// Is --p_raw set?
    bool IsStore_p_raw()                const {return Store_p_raw;};
    /// Is --p_rms set?
    bool IsStore_p_rms()                const {return Store_p_rms;};
    /// Is --p_max set?
    bool IsStore_p_max()                const {return Store_p_max;};
    /// Is --p_final set?
    bool IsStore_p_final()              const {return Store_p_final;};
    
    /// Is --u_raw set?
    bool IsStore_u_raw()                const {return Store_u_raw;};
    /// Is --u_rms set?
    bool IsStore_u_rms()                const {return Store_u_rms;};
    /// Is --u_max set?
    bool IsStore_u_max()                const {return Store_u_max;};    
    /// Is --u_final set?
    bool IsStore_u_final()              const {return Store_u_final;};
    
    /// Is --I_avg set
    bool IsStore_I_avg()                const {return Store_I_avg;};
    /// Is --I_max set
    bool IsStore_I_max()                const {return Store_I_max;};

    /// Print usage and exit
    void PrintUsageAndExit();   
    /// Print setup
    void PrintSetup();          
    /// Parse command line
    void ParseCommandLine(int argc, char** argv);    
    
    
protected:
    /// Copy constructor not allowed for public
    TCommandLineParameters(const TCommandLineParameters& src);
    
    /// operator = not allowed for public
    TCommandLineParameters& operator = (const TCommandLineParameters& src);

private:
    /// Input file name
    std::string InputFileName;
    /// Output file name
    std::string OutputFileName;
    
    /// NumberOfThreads value
    int         NumberOfThreads;
    /// VerboseInterval value
    int         VerboseInterval;
    /// CompressionLevel value
    int         CompressionLevel;
    
    /// BenchmarkFlag value
    bool        BenchmarkFlag;    
    /// BenchmarkTimeStepsCount value
    int         BenchmarkTimeStepsCount;
    /// PrintVersion value
    bool        PrintVersion;    
    
    /// Store_p_raw value
    bool        Store_p_raw;
    /// Store_p_rms value
    bool        Store_p_rms;
    /// Store_p_max value
    bool        Store_p_max;
    /// Store_p_final value
    bool        Store_p_final;
    
    /// Store_u_raw value
    bool        Store_u_raw;
    /// Store_u_rms value
    bool        Store_u_rms;
    /// Store_u_max value
    bool        Store_u_max;
    /// Store_u_final value
    bool        Store_u_final;
    
    /// Store_I_avg value
    bool        Store_I_avg;
    /// Store_I_max value
    bool        Store_I_max;
    /// StartTimeStep value
    int         StartTimeStep;
    
    

    /// Default compression level 
    static const int DefaultCompressionLevel = 3;
    /// Default verbose interval
    static const int DefaultVerboseInterval  = 5;
    
    
};// end of class TCommandLineParameters

#endif	/* TCOMMANDLINESPARAMETERS_H */

