/**
 * @file        DimensionSizes.h
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au   \n
 * 
 * @brief       The header file containing the structure with 3D dimension sizes
 * 
 * @version     kspaceFirstOrder3D 2.13
 * @date        9 August 2011, 12:34     (created) \n
 *              14 September 2012, 14:20 (revised)
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

/**
 * @var SSE_ALIGNMENT
 * @brief memory alignment for SSE2 (16B)
 */
const int SSE_ALIGNMENT = 16;

/**
 * @struct TDimensionSizes
 * @brief  Structure with  3D dimension sizes  
 */
struct TDimensionSizes {
    /// X dimension size
    size_t X;
    /// Y dimension size
    size_t Y; 
    /// Z dimension size            
    size_t Z; 
    
    /// default constructor
    TDimensionSizes(): X(0), Y(0), Z(0) {};   
    
    /**
     * @brief Constructor 
     * @param [in] x
     * @param [in] y
     * @param [in] z
     */
    TDimensionSizes(size_t x, size_t y, size_t z): X(x), Y(y), Z(z) {};    
   
    /**
     * @brief Copy constructor
     * @param src
     */
    TDimensionSizes(const TDimensionSizes & src) :
            X(src.X), Y(src.Y), Z(src.Z) {};
            
    /**
     * @brief Get element count
     * @return element count
     */
    size_t GetElementCount () const {return X * Y * Z; };
    
    
    /**
     * Operator = 
     * @param src
     * @return new value of the element
     */
     TDimensionSizes& operator = (const TDimensionSizes & src){
         if (this != &src){
                X  = src.X;
                Y  = src.Y;
                Z  = src.Z;                
         }    
        return *this;
      };
    
     
     
    /**
     * @brief Operator ==
     * @param [in] other     - the second operand to compare with
     * @return          true if ==
     */
    bool operator== (const TDimensionSizes &other) const {
        return ((X == other.X) && (Y == other.Y) &&  (Z == other.Z));
        
    }
    
    /**
     * @brief Operator !=
     * @param [in] other     - the second operand to compare with
     * @return          true if !=
     */
    bool operator!= (const TDimensionSizes &other) const {
        return ((X != other.X) || (Y != other.Y) ||  (Z != other.Z));
  }
    
    
}; // end of TDimensionSizes
//------------------------------------------------------------------------------


#endif	/* DIMENSIONSIZES_H */

