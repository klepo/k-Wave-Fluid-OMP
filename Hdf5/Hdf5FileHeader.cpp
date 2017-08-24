/**
 * @file      Hdf5FileHeader.cpp
 *
 * @author    Jiri Jaros \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            jarosjir@fit.vutbr.cz
 *
 * @brief     The implementation of the class responsible for working with file headers.
 *
 * @version   kspaceFirstOrder3D 32.16
 *
 * @date      24 August    2017, 09:51 (created) \n
 *            24 August    2017, 11:07 (revised)
 *
 * @copyright Copyright (C) 2017 Jiri Jaros and Bradley Treeby.
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

#include <iostream>
#include <stdexcept>
#include <ctime>

// Linux build
#ifdef __linux__
  #include <unistd.h>
#endif

//Windows 64 build
#ifdef _WIN64
  #include<stdio.h>
  #include<Winsock2.h>
  #pragma comment(lib, "Ws2_32.lib")
#endif

#include <Hdf5/Hdf5FileHeader.h>

#include <Parameters/Parameters.h>
#include <Utils/ErrorMessages.h>


//--------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------- Constants -----------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

const string THDF5_FileHeader::HDF5_FileTypesNames[]  = {"input","output", "checkpoint", "unknown"};

const string THDF5_FileHeader::HDF5_MajorFileVersionsNames[] = {"1"};

const string THDF5_FileHeader::HDF5_MinorFileVersionsNames[] = {"0","1"};





//----------------------------------------------------------------------------//
//------------------------    THDF5_File_Header ------------------------------//
//-----------------------------    Public     --------------------------------//
//----------------------------------------------------------------------------//

/**
 * Constructor.
 */
THDF5_FileHeader::THDF5_FileHeader()
{
  HDF5_FileHeaderValues.clear();
  PopulateHeaderFileMap();
}// end of constructor
//------------------------------------------------------------------------------


/**
 * Copy constructor.
 * @param [in] other
 */
THDF5_FileHeader::THDF5_FileHeader(const THDF5_FileHeader & other)
{
  HDF5_FileHeaderValues = other.HDF5_FileHeaderValues;
  HDF5_FileHeaderNames  = other.HDF5_FileHeaderNames;
}// end of copy constructor
//------------------------------------------------------------------------------


/**
 * Destructor.
 *
 */
THDF5_FileHeader::~THDF5_FileHeader()
{
  HDF5_FileHeaderValues.clear();
  HDF5_FileHeaderNames.clear();
}// end of destructor
//------------------------------------------------------------------------------



/**
 * Read header from the input file.
 * @param [in] InputFile - Input file to read from
 * @throw ios:failure when en error happens
 */
void THDF5_FileHeader::ReadHeaderFromInputFile(THDF5_File & InputFile)
{
  // Get file root handle
  hid_t FileRootHandle = InputFile.GetRootGroup();
  // read file type
  HDF5_FileHeaderValues[hdf5_fhi_file_type] =
          InputFile.ReadStringAttribute(FileRootHandle,"/",
                                        HDF5_FileHeaderNames[hdf5_fhi_file_type].c_str());

  if (GetFileType() == hdf5_ft_input)
  {
    HDF5_FileHeaderValues[hdf5_fhi_created_by]
            = InputFile.ReadStringAttribute(FileRootHandle, "/",
                                            HDF5_FileHeaderNames[hdf5_fhi_created_by].c_str());
    HDF5_FileHeaderValues[hdf5_fhi_creation_date]
            = InputFile.ReadStringAttribute(FileRootHandle, "/",
                                            HDF5_FileHeaderNames[hdf5_fhi_creation_date].c_str());
    HDF5_FileHeaderValues[hdf5_fhi_file_description]
            = InputFile.ReadStringAttribute(FileRootHandle, "/",
                                            HDF5_FileHeaderNames[hdf5_fhi_file_description].c_str());
    HDF5_FileHeaderValues[hdf5_fhi_major_version]
            = InputFile.ReadStringAttribute(FileRootHandle, "/",
                                            HDF5_FileHeaderNames[hdf5_fhi_major_version].c_str());
    HDF5_FileHeaderValues[hdf5_fhi_minor_version]
            = InputFile.ReadStringAttribute(FileRootHandle, "/",
                                            HDF5_FileHeaderNames[hdf5_fhi_minor_version].c_str());
  }
  else
  {
    throw ios::failure(kErrFmtBadInputFileType);
  }
}// end of ReadHeaderFromInputFile
//------------------------------------------------------------------------------


/**
 * Read header from output file (necessary for checkpoint-restart).
 * Read only execution times (the others are read from the input file, or
 * calculated based on the very last leg of the simulation).
 * This function is called only if checkpoint-restart is enabled.
 * @param [in] OutputFile
 * @throw ios:failure when en error happens
 */
void THDF5_FileHeader::ReadHeaderFromOutputFile(THDF5_File & OutputFile)
{
  // Get file root handle
  hid_t FileRootHandle = OutputFile.GetRootGroup();

  HDF5_FileHeaderValues[hdf5_fhi_file_type]
          = OutputFile.ReadStringAttribute(FileRootHandle, "/",
                                           HDF5_FileHeaderNames[hdf5_fhi_file_type].c_str());

  if (GetFileType() == hdf5_ft_output)
  {
    HDF5_FileHeaderValues[hdf5_fhi_total_execution_time]
            = OutputFile.ReadStringAttribute(FileRootHandle, "/",
                                             HDF5_FileHeaderNames[hdf5_fhi_total_execution_time].c_str());
    HDF5_FileHeaderValues[hdf5_fhi_data_load_time]
            = OutputFile.ReadStringAttribute(FileRootHandle, "/",
                                             HDF5_FileHeaderNames[hdf5_fhi_data_load_time].c_str());
    HDF5_FileHeaderValues[hdf5_fhi_preprocessing_time]
            = OutputFile.ReadStringAttribute(FileRootHandle, "/",
                                             HDF5_FileHeaderNames[hdf5_fhi_preprocessing_time].c_str());
    HDF5_FileHeaderValues[hdf5_fhi_simulation_time]
            = OutputFile.ReadStringAttribute(FileRootHandle, "/",
                                             HDF5_FileHeaderNames[hdf5_fhi_simulation_time].c_str());
    HDF5_FileHeaderValues[hdf5_fhi_postprocessing_time]
            = OutputFile.ReadStringAttribute(FileRootHandle, "/",
                                             HDF5_FileHeaderNames[hdf5_fhi_postprocessing_time].c_str());
  }
  else
  {
    throw ios::failure(kErrFmtBadOutputFIleType);
  }
}// end of ReadHeaderFromOutputFile
//------------------------------------------------------------------------------

/**
 * Read the file header form the checkpoint file. We need the header to verify
 * the file version and type
 * @param [in] CheckpointFile
 * @throw ios:failure when en error happens
 */
void THDF5_FileHeader::ReadHeaderFromCheckpointFile(THDF5_File & CheckpointFile)
{
  // Get file root handle
  hid_t FileRootHandle = CheckpointFile.GetRootGroup();
  // read file type
  HDF5_FileHeaderValues[hdf5_fhi_file_type] =
          CheckpointFile.ReadStringAttribute(FileRootHandle,"/",
                                             HDF5_FileHeaderNames[hdf5_fhi_file_type].c_str());

  if (GetFileType() == hdf5_ft_checkpoint)
  {
    HDF5_FileHeaderValues[hdf5_fhi_created_by]
            = CheckpointFile.ReadStringAttribute(FileRootHandle, "/",
                                                 HDF5_FileHeaderNames[hdf5_fhi_created_by].c_str());
    HDF5_FileHeaderValues[hdf5_fhi_creation_date]
            = CheckpointFile.ReadStringAttribute(FileRootHandle, "/",
                                                 HDF5_FileHeaderNames[hdf5_fhi_creation_date].c_str());
    HDF5_FileHeaderValues[hdf5_fhi_file_description]
            = CheckpointFile.ReadStringAttribute(FileRootHandle, "/",
                                                 HDF5_FileHeaderNames[hdf5_fhi_file_description].c_str());
    HDF5_FileHeaderValues[hdf5_fhi_major_version]
            = CheckpointFile.ReadStringAttribute(FileRootHandle, "/",
                                                 HDF5_FileHeaderNames[hdf5_fhi_major_version].c_str());
    HDF5_FileHeaderValues[hdf5_fhi_minor_version]
            = CheckpointFile.ReadStringAttribute(FileRootHandle, "/",
                                                 HDF5_FileHeaderNames[hdf5_fhi_minor_version].c_str());
  }
  else
  {
    throw ios::failure(kErrFmtBadCheckpointFileType);
  }
}// end of ReadHeaderFromCheckpointFile
//------------------------------------------------------------------------------

/**
 * Write header into the output file.
 * @param [in] OutputFile
 */
void THDF5_FileHeader::WriteHeaderToOutputFile(THDF5_File & OutputFile)
{
  // Get file root handle
  hid_t FileRootHandle = OutputFile.GetRootGroup();

  for (map<THDF5_FileHeaderItems, string>::iterator it = HDF5_FileHeaderNames.begin();
       it != HDF5_FileHeaderNames.end();
       it++)
  {
    OutputFile.WriteStringAttribute(FileRootHandle,
                                    "/",
                                    it->second.c_str(),
                                    HDF5_FileHeaderValues[it->first]);
  }
}// end of WriteHeaderToOutputFile
//------------------------------------------------------------------------------


/**
 * Write header to the output file (only a subset of all possible fields are written).
 * @param [in] CheckpointFile
 */
void THDF5_FileHeader::WriteHeaderToCheckpointFile(THDF5_File & CheckpointFile)
{
  // Get file root handle
  hid_t FileRootHandle = CheckpointFile.GetRootGroup();

  // Write header
  CheckpointFile.WriteStringAttribute(FileRootHandle,
                                      "/",
                                      HDF5_FileHeaderNames [hdf5_fhi_file_type].c_str(),
                                      HDF5_FileHeaderValues[hdf5_fhi_file_type].c_str());

  CheckpointFile.WriteStringAttribute(FileRootHandle,
                                      "/",
                                      HDF5_FileHeaderNames [hdf5_fhi_created_by].c_str(),
                                      HDF5_FileHeaderValues[hdf5_fhi_created_by].c_str());

  CheckpointFile.WriteStringAttribute(FileRootHandle,
                                      "/",
                                      HDF5_FileHeaderNames [hdf5_fhi_creation_date].c_str(),
                                      HDF5_FileHeaderValues[hdf5_fhi_creation_date].c_str());

  CheckpointFile.WriteStringAttribute(FileRootHandle,
                                      "/",
                                      HDF5_FileHeaderNames [hdf5_fhi_file_description].c_str(),
                                      HDF5_FileHeaderValues[hdf5_fhi_file_description].c_str());

  CheckpointFile.WriteStringAttribute(FileRootHandle,                                      "/",
                                      HDF5_FileHeaderNames [hdf5_fhi_major_version].c_str(),
                                      HDF5_FileHeaderValues[hdf5_fhi_major_version].c_str());

  CheckpointFile.WriteStringAttribute(FileRootHandle,
                                      "/",
                                      HDF5_FileHeaderNames [hdf5_fhi_minor_version].c_str(),
                                      HDF5_FileHeaderValues[hdf5_fhi_minor_version].c_str());
}// end of WriteHeaderToCheckpointFile
//------------------------------------------------------------------------------

/**
 * Get File type.
 * @return FileType
 */
THDF5_FileHeader::THDF5_FileType  THDF5_FileHeader::GetFileType()
{
  for (int i = hdf5_ft_input; i < hdf5_ft_unknown ; i++)
  {
    if (HDF5_FileHeaderValues[hdf5_fhi_file_type] == HDF5_FileTypesNames[static_cast<THDF5_FileType >(i)])
    {
      return static_cast<THDF5_FileType >(i);
    }
  }

  return THDF5_FileHeader::hdf5_ft_unknown;
}// end of GetFileType
//------------------------------------------------------------------------------

/**
 * Set file type.
 * @param [in] FileType
 */
void THDF5_FileHeader::SetFileType(const THDF5_FileHeader::THDF5_FileType FileType)
{
  HDF5_FileHeaderValues[hdf5_fhi_file_type] = HDF5_FileTypesNames[FileType];
}// end of SetFileType
//------------------------------------------------------------------------------


/**
 * Get File Version an enum
 * @return file version as an enum
 */
THDF5_FileHeader::THDF5_FileVersion THDF5_FileHeader::GetFileVersion()
{
  if ((HDF5_FileHeaderValues[hdf5_fhi_major_version] == HDF5_MajorFileVersionsNames[0]) &&
      (HDF5_FileHeaderValues[hdf5_fhi_minor_version] == HDF5_MinorFileVersionsNames[0]))
  {
    return hdf5_fv_10;
  }

  if ((HDF5_FileHeaderValues[hdf5_fhi_major_version] == HDF5_MajorFileVersionsNames[0]) &&
      (HDF5_FileHeaderValues[hdf5_fhi_minor_version] == HDF5_MinorFileVersionsNames[1]))
  {
    return hdf5_fv_11;
  }

  return hdf5_fv_unknown;
}// end of GetFileVersion
//------------------------------------------------------------------------------

/**
 * Set actual date and time.
 *
 */
void THDF5_FileHeader::SetActualCreationTime()
{
  struct tm *current;
  time_t now;
  time(&now);
  current = localtime(&now);

  char DateString[20];

  sprintf(DateString, "%02i/%02i/%02i, %02i:%02i:%02i",
          current->tm_mday, current->tm_mon+1, current->tm_year-100,
          current->tm_hour, current->tm_min, current->tm_sec);

  HDF5_FileHeaderValues[hdf5_fhi_creation_date] = DateString;
}// end of SetCreationTime
//------------------------------------------------------------------------------

/**
 * Set Host name.
 *
 */
void THDF5_FileHeader::SetHostName()
{
  char   HostName[256];

  //Linux build
  #ifdef __linux__
    gethostname(HostName, 256);
  #endif

  //Windows build
  #ifdef _WIN64
    WSADATA wsaData;
    WSAStartup(MAKEWORD(2, 2), &wsaData);
	  gethostname(HostName, 256);

    WSACleanup();
  #endif

  HDF5_FileHeaderValues[hdf5_fhi_host_name] = HostName;
}// end of SetHostName
//------------------------------------------------------------------------------


/**
 * Set memory consumption.
 * @param [in] TotalMemory
 */
void THDF5_FileHeader::SetMemoryConsumption(size_t TotalMemory)
{
  char Text[20] = "";
  sprintf(Text, "%ld MB",TotalMemory);

  HDF5_FileHeaderValues[hdf5_fhi_total_memory_consumption]     = Text;

  sprintf(Text, "%ld MB",TotalMemory / TParameters::GetInstance()->GetNumberOfThreads());
  HDF5_FileHeaderValues[hdf5_fhi_peak_core_memory_consumption] = Text;
}// end of SetMemoryConsumption
//------------------------------------------------------------------------------


/**
 * Set execution times in file header.
 * @param [in] TotalTime
 * @param [in] LoadTime
 * @param [in] PreProcessingTime
 * @param [in] SimulationTime
 * @param [in] PostprocessingTime
 */
void THDF5_FileHeader::SetExecutionTimes(const double TotalTime,
                                         const double LoadTime,
                                         const double PreProcessingTime,
                                         const double SimulationTime,
                                         const double PostprocessingTime)
{
  char Text [30] = "";

  sprintf(Text,"%8.2fs", TotalTime);
  HDF5_FileHeaderValues[hdf5_fhi_total_execution_time] = Text;

  sprintf(Text,"%8.2fs", LoadTime);
  HDF5_FileHeaderValues[hdf5_fhi_data_load_time] = Text;

  sprintf(Text,"%8.2fs", PreProcessingTime);
  HDF5_FileHeaderValues[hdf5_fhi_preprocessing_time] = Text;


  sprintf(Text,"%8.2fs", SimulationTime);
  HDF5_FileHeaderValues[hdf5_fhi_simulation_time] = Text;

  sprintf(Text,"%8.2fs", PostprocessingTime);
  HDF5_FileHeaderValues[hdf5_fhi_postprocessing_time] = Text;
}// end of SetExecutionTimes
//------------------------------------------------------------------------------

/**
 * Get execution times stored in the output file header
 * @param [out] TotalTime
 * @param [out] LoadTime
 * @param [out] PreProcessingTime
 * @param [out] SimulationTime
 * @param [out] PostprocessingTime
 */
void THDF5_FileHeader::GetExecutionTimes(double& TotalTime,
                                         double& LoadTime,
                                         double& PreProcessingTime,
                                         double& SimulationTime,
                                         double& PostprocessingTime)
{
  TotalTime          = atof(HDF5_FileHeaderValues[hdf5_fhi_total_execution_time].c_str());
  LoadTime           = atof(HDF5_FileHeaderValues[hdf5_fhi_data_load_time].c_str());
  PreProcessingTime  = atof(HDF5_FileHeaderValues[hdf5_fhi_preprocessing_time].c_str());
  SimulationTime     = atof(HDF5_FileHeaderValues[hdf5_fhi_simulation_time].c_str());
  PostprocessingTime = atof(HDF5_FileHeaderValues[hdf5_fhi_postprocessing_time].c_str());
}// end of GetExecutionTimes
//------------------------------------------------------------------------------

/**
 * Set Number of cores.
 *
 */
void THDF5_FileHeader::SetNumberOfCores()
{
  char Text[12] = "";
  sprintf(Text, "%ld",TParameters::GetInstance()->GetNumberOfThreads());

  HDF5_FileHeaderValues[hdf5_fhi_number_of_cores] = Text;
}// end of SetNumberOfCores
//------------------------------------------------------------------------------

//----------------------------------------------------------------------------//
//------------------------    THDF5_File_Header ------------------------------//
//---------------------------    Protected     --------------------------------//
//----------------------------------------------------------------------------//


/**
 * Create map with names for the header.
 *
 */
void THDF5_FileHeader::PopulateHeaderFileMap()
{
  HDF5_FileHeaderNames.clear();

  HDF5_FileHeaderNames[hdf5_fhi_created_by]                   = "created_by";
  HDF5_FileHeaderNames[hdf5_fhi_creation_date]                = "creation_date";
  HDF5_FileHeaderNames[hdf5_fhi_file_description]             = "file_description";
  HDF5_FileHeaderNames[hdf5_fhi_major_version]                = "major_version";
  HDF5_FileHeaderNames[hdf5_fhi_minor_version]                = "minor_version";
  HDF5_FileHeaderNames[hdf5_fhi_file_type]                    = "file_type";

  HDF5_FileHeaderNames[hdf5_fhi_host_name]                    = "host_names";
  HDF5_FileHeaderNames[hdf5_fhi_number_of_cores]              = "number_of_cpu_cores" ;
  HDF5_FileHeaderNames[hdf5_fhi_total_memory_consumption]     = "total_memory_in_use";
  HDF5_FileHeaderNames[hdf5_fhi_peak_core_memory_consumption] = "peak_core_memory_in_use";

  HDF5_FileHeaderNames[hdf5_fhi_total_execution_time]         = "total_execution_time";
  HDF5_FileHeaderNames[hdf5_fhi_data_load_time]               = "data_loading_phase_execution_time";
  HDF5_FileHeaderNames[hdf5_fhi_preprocessing_time]           = "pre-processing_phase_execution_time";
  HDF5_FileHeaderNames[hdf5_fhi_simulation_time]              = "simulation_phase_execution_time";
  HDF5_FileHeaderNames[hdf5_fhi_postprocessing_time]          = "post-processing_phase_execution_time";
}// end of PopulateHeaderFileMap
//------------------------------------------------------------------------------
