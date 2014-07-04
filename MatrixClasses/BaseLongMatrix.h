/**
 * @file        BaseLongMatrix.h
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au   \n
 * 
 * @brief       The header file containing the base class for 
 *              64b-wide integers (long for Linux/ size_t for Windows)
 
 * 
 * @version     kspaceFirstOrder3D 2.14
 * @date        26 July 2011, 2:17        (created) \n
 *              04 July 2014, 13:09       (revised) \n
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


#ifndef BASELONGMATRIXDATA_H
#define	BASELONGMATRIXDATA_H


#include <MatrixClasses/BaseMatrix.h>

#include <Utils/DimensionSizes.h>


using namespace std;

/**
 * @class TBaseLongMatrix
 * @brief Abstract base class for long based matrices defining basic interface.
 *        Higher dimensional matrices stored as 1D arrays, row-major order..
 */
class TBaseLongMatrix : public TBaseMatrix{
public:
    /// Default constructor
    TBaseLongMatrix(): TBaseMatrix(),  pTotalElementCount(0), pTotalAllocatedElementCount(0),
                           pDimensionSizes(), pDataRowSize(0), p2DDataSliceSize (0),  pMatrixData (NULL) 
    {};
    
    
    /// Get dimension sizes of the matrix
    inline struct TDimensionSizes GetDimensionSizes() const {
        return pDimensionSizes;
    }
    
    /// Get total element count of the matrix 
    virtual size_t GetTotalElementCount() const {
        return pTotalElementCount;
    };        
            
    /// Get total allocated element count (might differ from total element count used for the simulation because of padding).
    virtual size_t GetTotalAllocatedElementCount() const {
        return pTotalAllocatedElementCount;
    };        
    
    /// Destructor
    virtual ~TBaseLongMatrix(){};
        
    
    /// Zero all elements of the matrix (NUMA first touch).
    virtual void ZeroMatrix();
    
    /// Get raw data out of the class (for direct kernel access).
    virtual long* GetRawData() {return pMatrixData;}
    
    /// Get raw data out of the class (for direct kernel access).
    virtual const long* GetRawData() const 
    { 
      return pMatrixData;
    }

protected:
    /// Total number of elements
    size_t pTotalElementCount;
    /// Total number of allocated elements (the array size).
    size_t pTotalAllocatedElementCount;
            
    /// Dimension sizes 
    struct TDimensionSizes pDimensionSizes;     
 
    /// Size of 1D row in X dimension
    size_t pDataRowSize;
    /// Size of 2D slab (X,Y)
    size_t p2DDataSliceSize;
    
    /// Raw matrix data
    long* pMatrixData;
    

    /// Memory allocation
    virtual void AllocateMemory();

    /// Memory deallocation        
    virtual void FreeMemory() ;
 
     /// Copy constructor is not directly allowed 
    TBaseLongMatrix(const TBaseLongMatrix& src);
    /// operator =  is not directly allowed
    TBaseLongMatrix & operator =(const TBaseLongMatrix& src);
    

private:

};// end of TBaseLongMatrix
//------------------------------------------------------------------------------

#endif	/* TBASEINTDATA_H */

