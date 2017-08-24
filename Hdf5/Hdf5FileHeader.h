/**
 * @file      Hdf5FileHeader.h
 *
 * @author    Jiri Jaros \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            jarosjir@fit.vutbr.cz
 *
 * @brief     The header file containing the class processing file headers.
 *            Detail about the file header are described below.
 *
 * @version   kspaceFirstOrder3D 2.16
 *
 * @date      24 August    2017, 09:51 (created) \n
 *            24 August    2017, 11:07 (revised)
 *
 * @section   Hdf5FileHeader HDF5 File Header Structure
 *
 *The header includes following information
 *
 * \verbatim
+----------------------------------------------------------------------------------------------------------------------+
|                                           Input File / Checkpoint File Header                                        |
+----------------------------------------------------------------------------------------------------------------------+
| created_by                              Short description of the tool that created this file                         |
| creation_date                           Date when the file was created                                               |
| file_description                        Short description of the content of the file (e.g. simulation name)          |
| file_type                               Type of the file (input)                                                     |
| major_version                           Major version of the file definition (1)                                     |
| minor_version                           Minor version of the file definition (1)                                     |
+----------------------------------------------------------------------------------------------------------------------+
 \endverbatim
 *
 * \verbatim
+----------------------------------------------------------------------------------------------------------------------+
|                                                    Output File Header                                                |
+----------------------------------------------------------------------------------------------------------------------+
| created_by                              Short description of the tool that created this file                         |
| creation_date                           Date when the file was created                                               |
| file_description                        Short description of the content of the file (e.g. simulation name)          |
| file_type                               Type of the file (output)                                                    |
| major_version                           Major version of the file definition (1)                                     |
| minor_version                           Minor version of the file definition (1)                                     |
+----------------------------------------------------------------------------------------------------------------------+
| host_names                              List of hosts (computer names) the simulation was executed on                |
| number_of_cpu_cores                     Number of CPU cores used for the simulation                                  |
| data_loading_phase_execution_time       Time taken to load data from the file                                        |
| pre-processing_phase_execution_time     Time taken to pre-process data                                               |
| simulation_phase_execution_time         Time taken to run the simulation                                             |
| post-processing_phase_execution_time    Time taken to complete the post-processing phase                             |
| total_execution_time                    Total execution time                                                         |
| peak_core_memory_in_use                 Peak memory required per core during the simulation                          |
| total_memory_in_use Total               Peak memory in use                                                           |
+----------------------------------------------------------------------------------------------------------------------+
 \endverbatim
 *
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

#ifndef HDF5_FILE_HEADER_H
#define HDF5_FILE_HEADER_H

#include <map>

#include <Hdf5/Hdf5File.h>
#include <Utils/DimensionSizes.h>
#include <Utils/MatrixNames.h>

/**
 * @class THDF5_FileHeader
 * @brief Class for HDF5 header
 * @details This class manages all information that can be stored in the input
 * output or checkpoint file.
 */
class THDF5_FileHeader
{
  public:

    /**
     * @enum THDF5_FileHeaderItems
     * @brief List of all header items
     */
    enum THDF5_FileHeaderItems{
                                hdf5_fhi_created_by                   =  0,
                                hdf5_fhi_creation_date                =  1,
                                hdf5_fhi_file_description             =  2,
                                hdf5_fhi_major_version                =  3,
                                hdf5_fhi_minor_version                =  4,
                                hdf5_fhi_file_type                    =  5,
                                hdf5_fhi_host_name                    =  6,
                                hdf5_fhi_total_memory_consumption     =  7,
                                hdf5_fhi_peak_core_memory_consumption =  8,
                                hdf5_fhi_total_execution_time         =  9,
                                hdf5_fhi_data_load_time               = 10,
                                hdf5_fhi_preprocessing_time           = 11,
                                hdf5_fhi_simulation_time              = 12,
                                hdf5_fhi_postprocessing_time          = 13,
                                hdf5_fhi_number_of_cores              = 14
                              };

    /**
     * @enum  THDF5_FileType
     * @brief HDF5 file type.
     */
    enum THDF5_FileType       {
                                hdf5_ft_input  = 0,
                                hdf5_ft_output = 1,
                                hdf5_ft_checkpoint = 2,
                                hdf5_ft_unknown = 3
                              };

    /**
     * @enum  THDF5_FileVersion
     * @brief HDF5 file version.
     */
    enum THDF5_FileVersion    {
                                hdf5_fv_10 = 0,
                                hdf5_fv_11 = 1,
                                hdf5_fv_unknown = 2
                              };


    /// Constructor.
    THDF5_FileHeader();
    /// Copy constructor.
    THDF5_FileHeader(const THDF5_FileHeader & other);
    /// Destructor.
    ~THDF5_FileHeader();

    /// Read header from the input file.
    void ReadHeaderFromInputFile(THDF5_File & InputFile);
    /// Read Header from output file (necessary for checkpoint-restart).
    void ReadHeaderFromOutputFile(THDF5_File & OutputFile);
    /// Read Header from checkpoint file (necessary for checkpoint-restart).
    void ReadHeaderFromCheckpointFile(THDF5_File & CheckpointFile);

    /// Write header to the output file
    void WriteHeaderToOutputFile(THDF5_File & OutputFile);
    /// Write header to the output file
    void WriteHeaderToCheckpointFile(THDF5_File & CheckpointFile);

    /**
     * @brief Set code name.
     * @details Set code name to the header.
     * @param [in] CodeName - code version
     */
    void SetCodeName(const string& CodeName)
    {
      HDF5_FileHeaderValues[hdf5_fhi_created_by] = CodeName;
    };

    /// Set creation time.
    void SetActualCreationTime();

    /**
     * @brief   Get string version of current Major version.
     * @details Get string version of current Major version.
     * @return  current version
     */
    static string GetCurrentHDF5_MajorVersion()
    {
      return HDF5_MajorFileVersionsNames[0];
    };

    /**
     * @brief   Get string version of current Minor version.
     * @details Get string version of current Minor version.
     * @return  current minor version
     */
    static string GetCurrentHDF5_MinorVersion()
    {
      return HDF5_MinorFileVersionsNames[1];
    };

    /// Set major file version.
    void SetMajorFileVersion()
    {
      HDF5_FileHeaderValues[hdf5_fhi_major_version] = GetCurrentHDF5_MajorVersion();
    };
    /// Set minor file version.
    void SetMinorFileVersion()
    {
      HDF5_FileHeaderValues[hdf5_fhi_minor_version] = GetCurrentHDF5_MinorVersion();
    };

    /// Set major file version in a string.
    THDF5_FileVersion GetFileVersion();

    /**
     * @brief   Check major file version.
     * @details Check major file version.
     * @return true if ok
     */
    bool CheckMajorFileVersion()
    {
      return (HDF5_FileHeaderValues[hdf5_fhi_major_version] == GetCurrentHDF5_MajorVersion());
    };
    /**
     * @brief   Check minor file version.
     * @details Check minor file version.
     * @return true if ok
     */
    bool CheckMinorFileVersion()
    {
      return (HDF5_FileHeaderValues[hdf5_fhi_minor_version] <= GetCurrentHDF5_MinorVersion());
    };

    /// Get File type.
    THDF5_FileHeader::THDF5_FileType GetFileType();
    /// Set file type.
    void SetFileType(const THDF5_FileHeader::THDF5_FileType FileType);

    /// Set host name.
    void SetHostName();
    /// Set memory consumption.
    void SetMemoryConsumption(const size_t TotalMemory);
    /// Set execution times.
    void SetExecutionTimes(const double TotalTime,
                           const double LoadTime,
                           const double PreProcessingTime,
                           const double SimulationTime,
                           const double PostprocessingTime);

    /// Get execution times stored in the output file header.
    void GetExecutionTimes(double& TotalTime,
                           double& LoadTime,
                           double& PreProcessingTime,
                           double& SimulationTime,
                           double& PostprocessingTime);
    /// Set number of cores.
    void SetNumberOfCores();

private:
    /// Populate the map with the header items.
    void PopulateHeaderFileMap();

    /// map for the header values.
    map<THDF5_FileHeaderItems, string> HDF5_FileHeaderValues;
    /// map for the header names.
    map<THDF5_FileHeaderItems, string> HDF5_FileHeaderNames;

    ///String representation of different file types.
    static const string HDF5_FileTypesNames[];
    /// String representations of Major file versions.
    static const string HDF5_MajorFileVersionsNames[];
    /// String representations of Major file versions.
    static const string HDF5_MinorFileVersionsNames[];

};// THDF5_FileHeader
//------------------------------------------------------------------------------


#endif /* HDF5_FILE_HEADER_H */

