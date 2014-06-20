/**
 * @file        LongMatrix.h
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au
 * 
 * @brief       The header file containing the class for 64b integer matrices
 * 
 * @version     kspaceFirstOrder3D 2.14
 * @date        26 July     2011, 15:16 (created) \n
 *              20 June     2014, 15:37 (revised)
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

#ifndef LONGMATRIXDATA_H
#define	LONGMATRIXDATA_H


#include <MatrixClasses/BaseLongMatrix.h>

#include <Utils/DimensionSizes.h>

/**
 * @class TLongMatrix
 * @brief The class for 64b integers. It is used for sensor_mask_index or sensor_corners_mask 
 * to get the address of sampled voxels. 
 */
class TLongMatrix : public TBaseLongMatrix
{
  public:
    
    /// Constructor 
    TLongMatrix(struct TDimensionSizes DimensionSizes);
    
    /// Destructor    
    virtual ~TLongMatrix() { FreeMemory(); };
     
    /// Read data from the HDF5 file
    virtual void ReadDataFromHDF5File(THDF5_File & HDF5_File, 
                                      const char * MatrixName);
    /// Write data into the HDF5 file
    virtual void WriteDataToHDF5File(THDF5_File & HDF5_File, 
                                     const char * MatrixName, 
                                     const int CompressionLevel);
    
    /**
     * Operator []
     * @param index - 1D index into the matrix
     * @return  Value of the index
     */
    long& operator [](const size_t& index) 
    {
      return pMatrixData[index]; 
    };
    
    /**
     * Operator [] const
     * @param index - 1D index into the matrix
     * @return  Value of the index
     */
    const long & operator [](const size_t& index) const 
    {
      return pMatrixData[index]; 
    };
        
    /**
     * Get the top left corner of the index-th cuboid 
     * @param [in] index - Id of the corner
     * @return the top left corner
     */
    TDimensionSizes GetTopLeftCorner(const size_t index) const 
    {
      size_t X =  pMatrixData[6*index];
      size_t Y =  pMatrixData[6*index+1];
      size_t Z =  pMatrixData[6*index+2];
                 
      return TDimensionSizes(X, Y, Z);
    }; 
    
    /**
     * Get the bottom right corner of the index-th cuboid 
     * @param [in] index -Id of the corner
     * @return the bottom right corner
     */
    TDimensionSizes GetBottomRightCorner(const size_t index) const 
    {
      size_t X =  pMatrixData[6*index + 3];
      size_t Y =  pMatrixData[6*index + 4];
      size_t Z =  pMatrixData[6*index + 5];
                 
      return TDimensionSizes(X, Y, Z);
    };                 
    
    ///  Recompute indices MATALAB->C++     
    void RecomputeIndicesToCPP();
    
    ///  Recompute indices C++ -> MATLAB
    void RecomputeIndicesToMatlab();
    
    /// Get the total number of elements to be sampled within all cuboids
    size_t GetTotalNumberOfElementsInAllCuboids() const;
    
    
    
  protected:
    /// Default constructor not allowed for public
    TLongMatrix()  : TBaseLongMatrix() {};
    
    /// Copy constructor not allowed for public
    TLongMatrix(const TLongMatrix& src);

    /// Operator =  not allowed for public
    TLongMatrix& operator = (const TLongMatrix& src);
    
  private:

    
};// end of TLongMatrixData
//------------------------------------------------------------------------------
#endif	/* LONGMATRIXDATA_H */

