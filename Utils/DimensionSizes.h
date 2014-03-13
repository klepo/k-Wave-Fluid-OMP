/**
 * @file        DimensionSizes.h
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au   \n
 * 
 * @brief       The header file containing the structure with 3D dimension sizes
 * 
 * @version     kspaceFirstOrder3D 2.14
 * @date        09 August     2011, 12:34     (created) \n
 *              12 March      2014, 14:10     (revised)
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


#ifndef DIMENSIONSIZES_H
#define	DIMENSIONSIZES_H

#include <stdlib.h>
#include <stdio.h>

using namespace std;

#ifdef __AVX__
/**
 * @var DATA_ALIGNMENT
 * @brief memory alignment for AVX(32B)
 */
const int DATA_ALIGNMENT  = 32;
#else

/**
 * @var DATA_ALIGNMENT
 * @brief memory alignment for SSE, SSE2, SSE3, SSE4 (16B)
 */const int DATA_ALIGNMENT  = 16;
#endif

/**
 * @struct TDimensionSizes
 * @brief  Structure with  4D dimension sizes (3 in space and 1 in time), 
 * The structure can be used for 3D (the time is then set to 1). \n
 * The structure contains only POD, so no C++ stuff is necessary
 */
struct TDimensionSizes
{
    /// X dimension size
    size_t X;
    /// Y dimension size
    size_t Y;
    /// Z dimension size            
    size_t Z;
    /// Number of time steps (for time series datasets)
    size_t T;

    /// default constructor
    TDimensionSizes() : X(0), Y(0), Z(0), T(0) {};
    
    /**
     * @brief Constructor 
     * @param [in] x 
     * @param [in] y
     * @param [in] z
     * @param [in] t
     */
    TDimensionSizes(size_t x, size_t y, size_t z, size_t t = 0) : X(x), Y(y), Z(z), T(t)   {};
    
    /**
     * @brief Get element count, in 3D only spatial domain, in 4D with time
     * @return spatial element count
     */
    size_t GetElementCount() const
    {
      if (Is3D()) return X * Y * Z;
      else return X * Y * Z * T;
    };
        
    /// Is it a 3D object 
    bool Is3D() const
    {
      return (T == 0);
    }
    
    /// is it a 3D object with time?    
    bool Is4D() const
    {
      return (T > 0);
    }
    /**
     * @brief Operator ==
     * @param [in] other     - the second operand to compare with
     * @return true if ==
     */    
    bool operator==(const TDimensionSizes &other) const
    {
      return ((X == other.X) && (Y == other.Y) && (Z == other.Z) && (T == other.T));
    }
        
    /**
     * @brief Operator !=
     * @param [in] other     - the second operand to compare with
     * @return true if !=
     */
    bool operator!=(const TDimensionSizes &other) const
    {
      return ((X != other.X) || (Y != other.Y) || (Z != other.Z) || (T != other.T));
    }

    /**
     * Operator - 
     * Get the size of the cube by subtracting two corners
     * @param op1 - usually bottom right corner
     * @param op2 - usually top left corner
     * @return  the size of the inner cuboid
     */
     inline friend TDimensionSizes operator-(const TDimensionSizes &op1, const TDimensionSizes &op2)
     {
       // +1 because of planes (10.10.1 - 60.40.1)
       if (op1.Is3D() && op2.Is3D()){
         return TDimensionSizes(op1.X - op2.X + 1, op1.Y - op2.Y + 1, op1.Z - op2.Z + 1);
       }
       else
       {
         return TDimensionSizes(op1.X - op2.X + 1, op1.Y - op2.Y + 1, op1.Z - op2.Z + 1, op1.T - op2.T + 1);
       }        
     }
     
}; // end of TDimensionSizes
//------------------------------------------------------------------------------


#endif	/* DIMENSIONSIZES_H */

