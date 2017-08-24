/**
 * @file        BaseMatrix.h
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The header file of the common ancestor of all matrix classes.
 *              A pure abstract class.
 *
 * @version     kspaceFirstOrder3D 2.16
 *
 * @date        11 July      2012, 11:34 (created) \n
 *              24 August    2017, 12:20 (revised)
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


#ifndef BASEMATRIX_H
#define	BASEMATRIX_H


#include <Utils/DimensionSizes.h>
#include <Hdf5/Hdf5File.h>

/**
 * @class TBaseMatrix
 * @brief Abstract base class, the common ancestor defining the common interface
 *      and allowing derived classes to be allocated, freed and loaded from the file
 *      using the Matrix container.
 *
 * @details Abstract base class, the common ancestor defining the common interface
 *      and allowing derived classes to be allocated, freed and loaded from the file
 *      using the Matrix container.
 */
class TBaseMatrix
{
  public:
    /// Default constructor
    TBaseMatrix() {};

    /// Get dimension sizes of the matrix
    virtual struct DimensionSizes GetDimensionSizes() const  = 0;

    /// Get total element count of the matrix
    virtual size_t GetTotalElementCount()              const = 0;

    /// Get total allocated element count (might differ from the total element count used for the simulation because of e.g. padding).
    virtual size_t GetTotalAllocatedElementCount()     const  = 0;

    /**
     * @brief   Read matrix from the HDF5 file
     * @details Read matrix from the HDF5 file
     * @param [in] HDF5_File  - Handle to the HDF5 file
     * @param [in] MatrixName - HDF5 dataset name to read from
     */
    virtual void ReadDataFromHDF5File(Hdf5File & HDF5_File,
                                      const char * MatrixName) = 0;

    /**
     * @brief   Write data into the HDF5 file
     * @details Write data into the HDF5 file
     * @param [in] HDF5_File        - Handle to the HDF5 file
     * @param [in] MatrixName       - HDF5 dataset name to write to
     * @param [in] CompressionLevel - Compression level for the HDF5 dataset
     */
    virtual void WriteDataToHDF5File(Hdf5File & HDF5_File,
                                     const char * MatrixName,
                                     const size_t CompressionLevel) = 0;

    /// Destructor
    virtual ~TBaseMatrix() {};

};// end of TBaseMatrix

#endif	/* BASEMATRIX_H */