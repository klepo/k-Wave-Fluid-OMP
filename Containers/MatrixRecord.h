/**
 * @file        MatrixRecord.h
 * @author      Jiri Jaros \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The header file containing metadata about matrices stored in the matrix container.
 *
 * @version     kspaceFirstOrder3D 2.16
 *
 * @date        27 August    2017, 08:54 (created) \n
 *              27 August    2017, 08:54 (revised)
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


#ifndef MATRIX_RECORD_H
#define MATRIX_RECORD_H

#include <string>

#include <MatrixClasses/BaseMatrix.h>
#include <Utils/MatrixNames.h>

/**
 * @struct TMatrixRecord
 * @brief  A structure storing details about the matrix. The matrix container
 * stores this structures.
 * @details A structure storing details about the matrix. The matrix container
 * stores the list of these records with the data.
 */
struct TMatrixRecord
{
  /**
   * @enum TMatrixDataType
   * @brief All possible types of the matrix.
   */
  enum TMatrixDataType {mdtReal, mdtComplex, mdtIndex, mdtFFTW, mdtUxyz};

  /// Pointer to the matrix object.
  BaseMatrix   * MatrixPtr;
  /// Matrix data type.
  TMatrixDataType MatrixDataType;
  /// Matrix dimension sizes.
  DimensionSizes dimensionSizes;
  /// Is the matrix content loaded from the HDF5 file.
  bool            LoadData;
  /// Is the matrix necessary to be preserver when checkpoint is enabled.
  bool            Checkpoint;
  /// HDF5 matrix name.
  std::string          HDF5MatrixName;

  /// Default constructor.
  TMatrixRecord() : MatrixPtr(NULL), MatrixDataType(mdtReal),
          dimensionSizes(), LoadData(false), Checkpoint(false),
          HDF5MatrixName("")
  {};

  /// Copy constructor.
  TMatrixRecord(const TMatrixRecord& src);

  /// operator =
  TMatrixRecord& operator = (const TMatrixRecord& src);

  /// Set all values of the record.
  void SetAllValues(BaseMatrix *          MatrixPtr,
                    const TMatrixDataType  MatrixDataType,
                    const DimensionSizes  dimensionSizes,
                    const bool             LoadData,
                    const bool             Checkpoint,
                    const std::string      HDF5MatrixName);

  // Destructor.
  virtual ~TMatrixRecord() {};
};// end of TMatrixRecord
//------------------------------------------------------------------------------

#endif	/* MATRIX_RECORD_H */
