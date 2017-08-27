/**
 * @file        MatrixRecord.cpp
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The implementation file containing metadata about matrices stored in the matrix container.
 *
 * @version     kspaceFirstOrder3D 2.16
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


#include <Containers/MatrixRecord.h>



//============================================================================//
//                              TMatrixRecord                                 //
//============================================================================//

//----------------------------------------------------------------------------//
//--------------------------- Public methods ---------------------------------//
//----------------------------------------------------------------------------//


/**
 * Copy constructor of TMatrixRecord.
 * @param [in] src
 */
TMatrixRecord::TMatrixRecord(const TMatrixRecord& src) :
        MatrixPtr(src.MatrixPtr),
        MatrixDataType(src.MatrixDataType),
        dimensionSizes(src.dimensionSizes),
        LoadData(src.LoadData),
        Checkpoint(src.Checkpoint),
        HDF5MatrixName(src.HDF5MatrixName)
{

}// end of TMatrixRecord
//------------------------------------------------------------------------------


/**
 * operator = of TMatrixRecord
 * @param  [in] src
 * @return this
 */
TMatrixRecord& TMatrixRecord::operator = (const TMatrixRecord& src)
{
  if (this != &src)
  {
    MatrixPtr       = src.MatrixPtr;
    MatrixDataType  = src.MatrixDataType;
    dimensionSizes  = src.dimensionSizes;
    LoadData        = src.LoadData;
    Checkpoint      = src.Checkpoint;
    HDF5MatrixName  = src.HDF5MatrixName;
  }

  return *this;
}// end of operator =
//------------------------------------------------------------------------------



/**
 * Set all values for the record
 * @param [in] MatrixPtr        - Pointer to the MatrixClass object
 * @param [in] MatrixDataType   - Matrix data type
 * @param [in] DimensionSizes   - Dimension sizes
 * @param [in] LoadData         - Load data from file?
 * @param [in] Checkpoint       - Checkpoint this matrix?
 * @param [in] HDF5MatrixName   - HDF5 matrix name
 */
void TMatrixRecord::SetAllValues(BaseMatrix *          MatrixPtr,
                                 const TMatrixDataType  MatrixDataType,
                                 const DimensionSizes   dimensionSizes,
                                 const bool             LoadData,
                                 const bool             Checkpoint,
                                 const std::string           HDF5MatrixName)
{
  this->MatrixPtr       = MatrixPtr;
  this->MatrixDataType  = MatrixDataType;
  this->dimensionSizes  = dimensionSizes;
  this->LoadData        = LoadData;
  this->Checkpoint      = Checkpoint;
  this->HDF5MatrixName  = HDF5MatrixName;
}// end of SetAllValues
//------------------------------------------------------------------------------


//----------------------------------------------------------------------------//
//------------------------- Protected methods --------------------------------//
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//-------------------------- Private methods ---------------------------------//
//----------------------------------------------------------------------------//
