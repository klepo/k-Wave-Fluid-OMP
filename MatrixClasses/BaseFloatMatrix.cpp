/**
 * @file        BaseFloatMatrix.cpp
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au   \n
 * 
 * @brief       The implementation file containing the base class for 
 *              single precisions floating point numbers (floats)
 * 
 * @version     kspaceFirstOrder3D 2.14
 * @date        11 July 2011, 12:13      (created) \n
 *              17 September 2012, 15:35 (revised)
 * 
 * @section license
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



#include <string.h>
#include <malloc.h>
#include <assert.h>

#include <MatrixClasses/BaseFloatMatrix.h>

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
 * 
 * Copy data from another matrix with same size.
 * 
 * @param [in] src - source matrix
 * 
 */
void TBaseFloatMatrix::CopyData(TBaseFloatMatrix & src){
        
    memcpy(pMatrixData,src.pMatrixData,sizeof(float)*pTotalAllocatedElementCount);
        
}// end of CopyDataSameSize
//------------------------------------------------------------------------------


/**
 * Zero all allocated elements in parallel. \n
 * Also work as the first touch strategy on NUMA machines
 * 
 */
void TBaseFloatMatrix::ZeroMatrix(){

    #pragma omp parallel for schedule (static)
    for (size_t i=0; i < pTotalAllocatedElementCount; i++){
        pMatrixData[i] = 0.0f;
    }
    
}// end of ZeroMatrix
//------------------------------------------------------------------------------


/**
 * Divide a scalar by the elements of matrix.
 * @param [in] scalar - scalar to be divided 
 * 
 */
void TBaseFloatMatrix::ScalarDividedBy(const float  scalar){

    #pragma omp parallel for schedule (static)
    for (size_t i=0; i < pTotalAllocatedElementCount; i++){
        pMatrixData[i] = scalar / pMatrixData[i];
        
    }
}// end of ScalarDividedBy
//------------------------------------------------------------------------------



//----------------------------------------------------------------------------//
//                              Implementation                                //
//                             protected methods                              //
//----------------------------------------------------------------------------//

/**
 * Memory allocation based on the total number of elements. \n
 * Memory is aligned by the SSE_ALIGNMENT and all elements are zeroed.
 */
void TBaseFloatMatrix::AllocateMemory(){
    /* No memory allocated before this function*/
    assert(pMatrixData == NULL);
    
    pMatrixData = ( float *) memalign(DATA_ALIGNMENT,pTotalAllocatedElementCount * sizeof (float));
    
    if (!pMatrixData) {
        fprintf(stderr,Matrix_ERR_FMT_NotEnoughMemory, "TBaseFloatMatrix");
        throw bad_alloc();
    }
    
    ZeroMatrix();
    
}// end of AllocateMemory
//------------------------------------------------------------------------------

/**
 * Free memory
 */
 void TBaseFloatMatrix::FreeMemory(){
     
    if (pMatrixData) free(pMatrixData);
         
    pMatrixData = NULL;

}// end of MemoryDealocation      
//------------------------------------------------------------------------------


//----------------------------------------------------------------------------//
//                              Implementation                                //
//                              private methods                               //
//----------------------------------------------------------------------------//
