/**
 * @file        CommandLineParameters.cpp
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The implementation file containing the command line parameters.
 *
 * @version     kspaceFirstOrder3D 2.16
 *
 * @date        29 August    2012, 11:25 (created) \n
 *              24 August    2017, 15:30 (revised)
 *
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


//Linux build
#ifdef __linux__
  #include <getopt.h>
#endif

//Windows build
#ifdef _WIN64
  #include <GetoptWin64/Getopt.h>
#endif

#include <stdio.h>
#include <string.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

#include <stdexcept>
#include <Parameters/CommandLineParameters.h>

#include <Utils/ErrorMessages.h>
//--------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------- Constants -----------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//



//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Public methods ---------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Print usage and exit.
 */
void CommandLineParameters::printUsage()
{
  printf("---------------------------------- Usage ---------------------------------\n");
  printf("Mandatory parameters:\n");
  printf("  -i <input_file_name>            : HDF5 input file\n");
  printf("  -o <output_file_name>           : HDF5 output file\n");
  printf("\n");
  printf("Optional parameters: \n");

#ifdef _OPENMP
  printf("  -t <num_threads>                : Number of CPU threads\n");
  printf("                                      (default = %d)\n", omp_get_num_procs());
#endif

  printf("  -r <interval_in_%%>              : Progress print interval\n");
  printf("                                      (default = %ld%%)\n", kDefaultProgressPrintInterval);
  printf("  -c <comp_level>                 : Output file compression level <0,9>\n");
  printf("                                      (default = %ld)\n", kDefaultCompressionLevel);
  printf("  --benchmark <steps>             : Run a specified number of time steps\n");
  printf("\n");
  printf("  --checkpoint_file <file_name>   : HDF5 checkpoint file\n");
  printf("  --checkpoint_interval <seconds> : Stop after a given number of seconds and\n");
  printf("                                      store the actual state\n");
  printf("\n");
  printf("  -h                              : Print help\n");
  printf("  --help                          : Print help\n");
  printf("  --version                       : Print version\n");
  printf("\n");
  printf("Output flags:\n");
  printf("  -p                              : Store acoustic pressure \n");
  printf("                                      (default if nothing else is on)\n");
  printf("                                      (the same as --p_raw)\n");
  printf("  --p_raw                         : Store raw time series of p (default)\n");
  printf("  --p_rms                         : Store rms of p\n");
  printf("  --p_max                         : Store max of p\n");
  printf("  --p_min                         : Store min of p\n");
  printf("  --p_max_all                     : Store max of p (whole domain)\n");
  printf("  --p_min_all                     : Store min of p (whole domain)\n");
  printf("  --p_final                       : Store final pressure field \n");
  printf("\n");
  printf("  -u                              : Store ux, uy, uz\n");
  printf("                                      (the same as --u_raw)\n");
  printf("  --u_raw                         : Store raw time series of ux, uy, uz\n");
  printf("  --u_non_staggered_raw           : Store non-staggered raw time series of\n");
  printf("                                      ux, uy, uz \n");
  printf("  --u_rms                         : Store rms of ux, uy, uz\n");
  printf("  --u_max                         : Store max of ux, uy, uz\n");
  printf("  --u_min                         : Store min of ux, uy, uz\n");
  printf("  --u_max_all                     : Store max of ux, uy, uz (whole domain)\n");
  printf("  --u_min_all                     : Store min of ux, uy, uz (whole domain)\n");
  printf("  --u_final                       : Store final acoustic velocity\n");
  printf("\n");
  printf("  --copy_sensor_mask              : Copy sensor mask to the output file\n");
  printf("\n");
  printf("  -s <timestep>                   : Time step when data collection begins\n");
  printf("                                      (default = 1)\n");
  printf("--------------------------------------------------------------------------\n");
  printf("\n");

  exit(EXIT_FAILURE);
}// end of printUsage
//----------------------------------------------------------------------------------------------------------------------

/**
 * Print out commandline parameters.
 */
void CommandLineParameters::printComandlineParamers()
{
  printf("List of enabled parameters:\n");

  printf("  Input  file               %s\n", mInputFileName.c_str());
  printf("  Output file               %s\n", mOutputFileName.c_str());
  printf("\n");
  printf("  Number of threads         %ld\n", mNumberOfThreads);
  printf("  Verbose interval[%%]       %ld\n", mProgressPrintInterval);
  printf("  Compression level         %ld\n", mCompressionLevel);
  printf("\n");
  printf("  Benchmark flag            %d\n", mBenchmarkFlag);
  printf("  Benchmark time steps      %ld\n", mBenchmarkTimeStepCount);
  printf("\n");
  printf("  Checkpoint_file           %s\n", mCheckpointFileName.c_str());
  printf("  Checkpoint_interval       %ld\n", mCheckpointInterval);
  printf("\n");
  printf("  Store p_raw               %d\n", mStorePressureRawFlag);
  printf("  Store p_rms               %d\n", mStorePressureRmsFlag);
  printf("  Store p_max               %d\n", mStorePressureMaxFlag);
  printf("  Store p_min               %d\n", mStorePressureMinFlag);
  printf("  Store p_max_all           %d\n", mStorePressureMaxAllFlag);
  printf("  Store p_min_all           %d\n", mStorePressureMinAllFlag);
  printf("  Store p_final             %d\n", mStorePressureFinalAllFlag);
  printf("\n");
  printf("  Store u_raw               %d\n", mStoreVelocityRawFlag);
  printf("  Store u_non_staggered_raw %d\n", mStoreVelocityNonStaggeredRawFlag);
  printf("  Store u_rms               %d\n", mStoreVelocityRmsFlag);
  printf("  Store u_max               %d\n", mStoreVelocityMaxFlag);
  printf("  Store u_min               %d\n", mStoreVelocityMinFlag);
  printf("  Store u_max_all           %d\n", mStoreVelocityMaxAllFlag);
  printf("  Store u_min_all           %d\n", mStoreVelocityMinAllFlag);
  printf("  Store u_final             %d\n", mStoreVelocityFinalAllFlag);
  printf("\n");
  printf("  Copy sensor mask          %d\n", mCopySensorMaskFlag);
  printf("\n");
  printf("  Collection begins at      %ld\n", mSamplingStartTimeStep + 1);
}// end of  printComandlineParamers
//----------------------------------------------------------------------------------------------------------------------

/**
 * Parse command line.
 */
void CommandLineParameters::parseCommandLine(int argc, char** argv)
{
  char c;
  int longIndex;
  bool CheckpointFlag = false;

#ifdef _OPENMP
  const char * shortOpts = "i:o:r:c:t:puhs:";
#else
  const char * shortOpts = "i:o:r:c:puhs:";
#endif

  const struct option longOpts[] =
  {
    { "benchmark",           required_argument, NULL, 0},
    { "help",                no_argument, NULL, 'h'},
    { "version",             no_argument, NULL, 0},
    { "checkpoint_file"    , required_argument, NULL, 0 },
    { "checkpoint_interval", required_argument, NULL, 0 },

    { "p_raw",               no_argument, NULL, 'p'},
    { "p_rms",               no_argument, NULL, 0},
    { "p_max",               no_argument, NULL, 0},
    { "p_min",               no_argument, NULL, 0},
    { "p_max_all",           no_argument, NULL, 0},
    { "p_min_all",           no_argument, NULL, 0},
    { "p_final",             no_argument, NULL, 0},

    { "u_raw",               no_argument, NULL, 'u'},
    { "u_non_staggered_raw", no_argument, NULL, 0},
    { "u_rms",               no_argument, NULL, 0},
    { "u_max",               no_argument, NULL, 0},
    { "u_min",               no_argument, NULL, 0},
    { "u_max_all",           no_argument, NULL, 0},
    { "u_min_all",           no_argument, NULL, 0},
    { "u_final",             no_argument, NULL, 0},

    { "copy_sensor_mask",    no_argument, NULL, 0},
    { NULL,                  no_argument, NULL, 0}
  };


   // Short parameters //
  while ((c = getopt_long (argc, argv, shortOpts, longOpts, &longIndex )) != -1)
  {
    switch (c)
    {
      case 'i':
      {
         mInputFileName = optarg;
         break;
      }

      case 'o':
      {
         mOutputFileName = optarg;
         break;
      }

      case 'r':
      {
        if ((optarg == NULL) || (atol(optarg) <= 0))
        {
          fprintf(stderr,"%s", kErrFmtNoProgressPrintInterval);
           printUsage();
        }
        else
        {
          mProgressPrintInterval = atol(optarg);
        }
        break;
      }

      case 't':
      {
        if ((optarg == NULL) || (atol(optarg) <= 0))
        {
          fprintf(stderr,"%s", kErrFmtInvalidNumberOfThreads);
          printUsage();
        }
        else
        {
          mNumberOfThreads = atol(optarg);
        }
        break;
      }

      case 'c':
      {
        if ((optarg == NULL) || (atol(optarg) < 0) || atol(optarg) > 9)
        {
          fprintf(stderr,"%s", kErrFmtNoCompressionLevel);
          printUsage();
        }
        else
        {
          mCompressionLevel = atol(optarg);
        }
         break;
      }

      case 'p':
      {
        mStorePressureRawFlag = true;
        break;
      }

      case 'u':
      {
        mStoreVelocityRawFlag = true;
        break;
      }

      case 'h':
      {
        printUsage();
        break;
      }

      case 's':
      {
        if ((optarg == NULL) || (atol(optarg) < 1))
        {
          fprintf(stderr,"%s", kErrFmtNoSamplingStartTimeStep);
          printUsage();
        }
        mSamplingStartTimeStep = (size_t) (atol(optarg) - 1);
        break;
      }

      // long option without a short arg
      case 0:
      {
        if( strcmp( "benchmark", longOpts[longIndex].name ) == 0 )
        {
          mBenchmarkFlag = true;
          if ((optarg == NULL) || (atol(optarg) <= 0))
          {
            fprintf(stderr,"%s", kErrFmtNoBenchmarkTimeStep);
            printUsage();
          }
          else
          {
             mBenchmarkTimeStepCount = atol(optarg);
          }
        }
        else if( strcmp( "checkpoint_file", longOpts[longIndex].name ) == 0 )
        {
          CheckpointFlag = true;
          if ((optarg == NULL))
          {
            fprintf(stderr,"%s", kErrFmtNoCheckpointFile);
            printUsage();
          }
          else
          {
            mCheckpointFileName = optarg;
          }
        }
        else if( strcmp( "checkpoint_interval", longOpts[longIndex].name ) == 0 )
        {
          CheckpointFlag = true ;
          if ((optarg == NULL) || (atol(optarg) <= 0))
          {
            fprintf(stderr,"%s", kErrFmtNoCheckpointInterval);
            printUsage();
          }
          else
          {
            mCheckpointInterval = atol(optarg);
          }
        }
        else if (strcmp("version", longOpts[longIndex].name) == 0)
        {
          mPrintVersionFlag = true;
          return;
        }

        //-- pressure related flags
        else if (strcmp("p_rms", longOpts[longIndex].name) == 0)
        {
          mStorePressureRmsFlag = true;
        }
        else if (strcmp("p_max", longOpts[longIndex].name) == 0)
        {
          mStorePressureMaxFlag = true;
        }
        else if (strcmp("p_min", longOpts[longIndex].name) == 0)
        {
          mStorePressureMinFlag = true;
        }
        else if (strcmp("p_max_all", longOpts[longIndex].name) == 0)
        {
          mStorePressureMaxAllFlag = true;
        }
        else if (strcmp("p_min_all", longOpts[longIndex].name) == 0)
        {
          mStorePressureMinAllFlag = true;
        }
        else if (strcmp("p_final", longOpts[longIndex].name) == 0)
        {
          mStorePressureFinalAllFlag = true;
        }

        //-- velocity related flags
        else if (strcmp("u_non_staggered_raw", longOpts[longIndex].name) == 0)
        {
          mStoreVelocityNonStaggeredRawFlag = true;
        }
        else if (strcmp("u_rms", longOpts[longIndex].name) == 0)
        {
          mStoreVelocityRmsFlag = true;
        }
        else if (strcmp("u_max", longOpts[longIndex].name) == 0)
        {
          mStoreVelocityMaxFlag = true;
        }
        else if (strcmp("u_min", longOpts[longIndex].name) == 0)
        {
          mStoreVelocityMinFlag = true;
        }
        else if (strcmp("u_max_all", longOpts[longIndex].name) == 0)
        {
          mStoreVelocityMaxAllFlag = true;
        }
        else if (strcmp("u_min_all", longOpts[longIndex].name) == 0)
        {
          mStoreVelocityMinAllFlag = true;
        }
        else if (strcmp("u_final", longOpts[longIndex].name) == 0)
        {
          mStoreVelocityFinalAllFlag = true;
        }

        else if (strcmp("copy_sensor_mask", longOpts[longIndex].name) == 0)
        {
          mCopySensorMaskFlag = true;
        }
        else
        {
          printUsage();
        }

        break;
      }
      default:
      {
        printUsage();
      }
    }
  }


  //-- Post checks --//
  if (mInputFileName == "")
  {
    fprintf(stderr, "%s", kErrFmtNoInputFile);
    printUsage();
  }


  if (mOutputFileName == "")
  {
    fprintf(stderr, "%s", kErrFmtNoOutputFile);
    printUsage();
  }

  if (CheckpointFlag)
  {
    if (mCheckpointFileName == "")
    {
      fprintf(stderr, "%s", kErrFmtNoCheckpointFile);
      printUsage();
    }
    if (mCheckpointInterval == 0)
    {
      fprintf(stderr, "%s", kErrFmtNoCheckpointInterval);
      printUsage();
    }
  }

  if (!(mStorePressureRawFlag     || mStorePressureRmsFlag     || mStorePressureMaxFlag   || mStorePressureMinFlag ||
        mStorePressureMaxAllFlag || mStorePressureMinAllFlag || mStorePressureFinalAllFlag ||
        mStoreVelocityRawFlag     || mStoreVelocityNonStaggeredRawFlag        ||
        mStoreVelocityRmsFlag     || mStoreVelocityMaxFlag     || mStoreVelocityMinFlag   ||
        mStoreVelocityMaxAllFlag || mStoreVelocityMinAllFlag || mStoreVelocityFinalAllFlag ))
  {
    mStorePressureRawFlag = true;
  }

}// end of parseCommandLine
//----------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------ Protected methods -------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Constructor.
 */
CommandLineParameters::CommandLineParameters() :
        mInputFileName(""), mOutputFileName (""), mCheckpointFileName(""),
        #ifdef _OPENMP
          mNumberOfThreads(omp_get_num_procs()),
        #else
          mNumberOfThreads(1),
        #endif

        mProgressPrintInterval(kDefaultProgressPrintInterval),
        mCompressionLevel (kDefaultCompressionLevel),
        mBenchmarkFlag (false), mBenchmarkTimeStepCount(0),
        mCheckpointInterval(0),
        mPrintVersionFlag (false),
        // output flags
        mStorePressureRawFlag(false), mStorePressureRmsFlag(false),
        mStorePressureMaxFlag(false), mStorePressureMinFlag(false),
        mStorePressureMaxAllFlag(false), mStorePressureMinAllFlag(false), mStorePressureFinalAllFlag(false),
        mStoreVelocityRawFlag(false), mStoreVelocityNonStaggeredRawFlag(false),
        mStoreVelocityRmsFlag(false), mStoreVelocityMaxFlag(false), mStoreVelocityMinFlag(false),
        mStoreVelocityMaxAllFlag(false), mStoreVelocityMinAllFlag(false), mStoreVelocityFinalAllFlag(false),
        mCopySensorMaskFlag(false),
        mSamplingStartTimeStep(0)
{

}// end of constructor
//----------------------------------------------------------------------------------------------------------------------
