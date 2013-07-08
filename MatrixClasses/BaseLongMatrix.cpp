/**
 * @file        BaseLongMatrix.cpp
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au   \n
 * 
 * @brief       The implementation file containing the base class for 
 *              64b-wide integers (long for Linux/ size_t for Windows)
 * 
 * @version     kspaceFirstOrder3D 2.14
 * @date        26 July 2011, 14:17 (created) \n
 *              17 September 2012, 15:35 (revised)
 * 
 * @section License
 * This file is part of the C++ extension of the k-Wave Toolbox (http://www.k-wave.org). .\n
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



#include <string.h>
#include <malloc.h>

#include <MatrixClasses/BaseLongMatrix.h>

#include <Utils/DimensionSizes.h>
#include <Utils/ErrorMessages.h>


//----------------------------------------------------------------------------//
//                              Constants                                     //
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
//                              Definitions                                   //
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
//                              Implementation                                //
//                              public methods                                //
//----------------------------------------------------------------------------//

/**
 *  * Zero all allocated elements.
 * 
 */
void TBaseLongMatrix::ZeroMatrix(){

    memset(pMatrixData,0,pTotalAllocatedElementCount*sizeof(long));
    
}// end of ZeroMatrix
//------------------------------------------------------------------------------



//----------------------------------------------------------------------------//
//                              Implementation                                //
//                             protected methods                              //
//----------------------------------------------------------------------------//


/**
 * Memory allocation based on the total number of elements. \n
 * Memory is aligned by the SSE_ALIGNMENT and all elements are zeroed.
 */
void TBaseLongMatrix::AllocateMemory(){
    /* No memory allocated before this function*/

    pMatrixData = (long *) memalign(SSE_ALIGNMENT,pTotalAllocatedElementCount * sizeof (long));
    
    if (!pMatrixData) {
        fprintf(stderr,Matrix_ERR_FMT_NotEnoughMemory, "TBaseLongMatrix");
        throw bad_alloc();   
    }    

    ZeroMatrix();
    
}// end of AllocateMemory
//------------------------------------------------------------------------------

/**
 * Free memory
 */
 void TBaseLongMatrix::FreeMemory(){
     
    if (pMatrixData) free(pMatrixData);
    pMatrixData = NULL;

}// end of MemoryDealocation      
//------------------------------------------------------------------------------




//----------------------------------------------------------------------------//
//                              Implementation                                //
//                              private methods                               //
//----------------------------------------------------------------------------//
