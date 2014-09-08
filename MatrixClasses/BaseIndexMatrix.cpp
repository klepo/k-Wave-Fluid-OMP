/**
 * @file        BaseIndexMatrix.cpp
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The implementation file containing the base class for index
 *              matrices (based on the size_t datatype)
 *
 * @version     kspaceFirstOrder3D 2.15
 *
 * @date        26 July      2011, 14:17 (created) \n
 *              01 September 2014, 15:55 (revised)
 *
 * @section License
 * This file is part of the C++ extension of the k-Wave Toolbox (http://www.k-wave.org). .\n
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



#include <string.h>
#include <immintrin.h>
#include <MatrixClasses/BaseIndexMatrix.h>

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
 *  Zero all allocated elements.
 *
 */
void TBaseIndexMatrix::ZeroMatrix()
{
  ///@TODO: This breaks the first touch policy! - however we don't know the distribution
  memset(pMatrixData, 0, pTotalAllocatedElementCount*sizeof(size_t));
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
void TBaseIndexMatrix::AllocateMemory()
{
  /* No memory allocated before this function*/

  pMatrixData = (size_t *) _mm_malloc(pTotalAllocatedElementCount * sizeof (size_t), DATA_ALIGNMENT);

  if (!pMatrixData)
  {
    fprintf(stderr,Matrix_ERR_FMT_NotEnoughMemory, "TBaseIndexMatrix");
    throw bad_alloc();
  }

  ZeroMatrix();
}// end of AllocateMemory
//------------------------------------------------------------------------------

/**
 * Free memory
 */
void TBaseIndexMatrix::FreeMemory()
{
  if (pMatrixData) _mm_free(pMatrixData);
  pMatrixData = NULL;
}// end of MemoryDealocation
//------------------------------------------------------------------------------




//----------------------------------------------------------------------------//
//                              Implementation                                //
//                              private methods                               //
//----------------------------------------------------------------------------//
