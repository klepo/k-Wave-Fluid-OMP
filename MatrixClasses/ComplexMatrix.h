/**
 * @file        ComplexMatrix.h
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The header file with the class for complex matrices.
 *
 * @version     kspaceFirstOrder3D 2.16
 *
 * @date        11 July      2011, 14:02 (created) \n
 *              03 September 2017, 13:19 (revised)
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

#ifndef COMPLEX_MATRIX_H
#define COMPLEX_MATRIX_H

#include <complex>

#include <MatrixClasses/BaseFloatMatrix.h>
#include <MatrixClasses/RealMatrix.h>

#include <Utils/DimensionSizes.h>

/// Datatype for complex single precision numbers.
using FloatComplex = std::complex<float>;


/**
 * @class   ComplexMatrix
 * @brief   The class for complex matrices.
 * @details The class for complex matrices.
 */
class ComplexMatrix : public BaseFloatMatrix
{
  public:

   /// Default constructor not allowed.
    ComplexMatrix() = delete;
    /**
     * @brief Constructor.
     * @param [in] dimensionSizes - Dimension sizes of the matrix.
     */
    ComplexMatrix(const DimensionSizes& dimensionSizes);
    /// Copy constructor not allowed.
    ComplexMatrix(const ComplexMatrix&) = delete;
    /// Destructor.
    virtual ~ComplexMatrix();

    /// Operator= is not allowed.
    ComplexMatrix& operator= (const ComplexMatrix&);

    /**
     * @brief   Read matrix from HDF5 file.
     * @details Read matrix from HDF5 file.
     * @param [in] file       - Handle to the HDF5 file.
     * @param [in] matrixName - HDF5 dataset name to read from.
     * @throw ios::failure    - If error occurred.
     */
    virtual void readData(Hdf5File&   file,
                          MatrixName& matrixName);

    /**
     * @brief   Write data into HDF5 file.
     * @details Write data into HDF5 file.
     * @param [in] file             - Handle to the HDF5 file
     * @param [in] matrixName       - HDF5 dataset name to write to.
     * @param [in] compressionLevel - Compression level for the HDF5 dataset.
     * @throw ios::failure          - If an error occurred.
     */
    virtual void writeData(Hdf5File&    file,
                           MatrixName&  matrixName,
                           const size_t compressionLevel);

    /**
     * @brief Get raw complex data out of the class (for direct kernel access).
     * @return Mutable matrix data
     */
    virtual FloatComplex* getComplexData()
    {
      return reinterpret_cast<FloatComplex*> (mData);
    };

    /**
     * @brief  Get raw complex data out of the class (for direct kernel access).
     * @return Imutable matrix data
     */
    virtual const FloatComplex* getComplexData() const
    {
      return reinterpret_cast<FloatComplex*> (mData);
    };

    /**
     * @brief  Operator [].
     * @param [in] index - 1D index into the matrix.
     * @return An element of the matrix.
     */
    inline FloatComplex& operator[](const size_t& index)
    {
      return reinterpret_cast<FloatComplex*>(mData)[index];
    };
    /**
     * @brief   Operator [], constant version.
     * @param [in] index - 1D index into the matrix.
     * @return An element of the matrix.
     */
    inline const FloatComplex& operator[](const size_t& index) const
    {
      return reinterpret_cast<FloatComplex*> (mData)[index];
    };

    /**
     * @brief   Get element from 3D matrix.
     * @details Get element from 3D matrix.
     * @param [in] x - x dimension
     * @param [in] y - y dimension
     * @param [in] z - z dimension
     * @return a complex element of the class
     */
    inline  FloatComplex& GetElementFrom3D(const size_t x,
                                           const size_t y,
                                           const size_t z)
    {
      return reinterpret_cast<FloatComplex*> (mData)[z * (mSlabSize>>1) + y * (mRowSize>>1) + x];
    };

    /**
     * @brief   Get element from 3D matrix, constant version.
     * @details Get element from 3D matrix, constant version.
     * @param [in] x - x dimension
     * @param [in] y - y dimension
     * @param [in] z - z dimension
     * @return a complex element of the class
     */
    inline const FloatComplex& GetElementFrom3D(const size_t x,
                                                const size_t y,
                                                const size_t z) const
    {
      return reinterpret_cast<FloatComplex*> (mData)[z * (mSlabSize >> 1) + y * (mRowSize >> 1) + x];
    };

  protected:

  private:
   /**
     * @brief Initialize dimension sizes
     * @param [in] dimensionSizes - Dimension sizes of the matrix.
     */
    void initDimensions(const DimensionSizes& dimensionSizes);

};// end of ComplexMatrix
//----------------------------------------------------------------------------------------------------------------------

#endif	/* COMPLEX_MATRIX_H */
