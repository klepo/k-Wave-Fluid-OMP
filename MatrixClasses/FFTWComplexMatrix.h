/**
 * @file        FFTWComplexMatrix.h
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au
 *
 * @brief       The header file containing the class that implements
 *              3D FFT using the FFTW interface
 *
 * @version     kspaceFirstOrder3D 2.14
 * @date        09 August    2011, 13:10 (created) \n
 *              07 July      2014, 18:01 (revised)
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

#ifndef FFTWCOMPLEXMATRIX_H
#define	FFTWCOMPLEXMATRIX_H


#include <fftw3.h>

#include <MatrixClasses/ComplexMatrix.h>

/**
 * @class TFFTWComplexMatrix
 * @brief Class implementing 3D Real-To-Complex and Complex-To-Real transforms
 *      using FFTW interface
 *
 */
class TFFTWComplexMatrix : public TComplexMatrix
{
  public:
    /// Constructor
    TFFTWComplexMatrix(struct TDimensionSizes DimensionSizes);
    /// Destructor
    virtual ~TFFTWComplexMatrix();

    /// Create FFTW plan for Real-to-Complex
    void CreateFFTPlan3D_R2C(TRealMatrix& InMatrix);
    /// Create FFTW plan for Complex-to-Real
    void CreateFFTPlan3D_C2R(TRealMatrix& OutMatrix);


    /// Compute 3D out-of-place Real-to-Complex FFT
    void Compute_FFT_3D_R2C(TRealMatrix& InMatrix);
    /// Compute 3D out-of-place Complex-to-Real FFT
    void Compute_iFFT_3D_C2R(TRealMatrix& OutMatrix);

    /// Export wisdom to the file
    static void ExportWisdom();
    /// Import wisdom from the file
    static void ImportWisdom();
    /// Destroy wisdom in the file (delete it)
    static void DeleteStoredWisdom();

protected:
    /// Default constructor not allowed for public
    TFFTWComplexMatrix() :
                 TComplexMatrix(),
                 fftw_plan_R2C(NULL),
                 fftw_plan_C2R(NULL)
    {};

    /// Copy constructor not allowed for public
    TFFTWComplexMatrix(const TFFTWComplexMatrix& src);

    /// Operator = not allowed for public
    TFFTWComplexMatrix & operator = (const TFFTWComplexMatrix& src);

    /// Get Wisdom file name, base on the checkpoint filename
    static string GetWisdomFileName();

    /// Allocate memory for the FFTW matrix
    virtual void AllocateMemory();
    /// Free memory of the FFTW matrix
    virtual void FreeMemory();

    /// FFTW plan flag
    static const unsigned TFFTWComplexMatrix_FFT_FLAG  = FFTW_MEASURE;

    static const string FFTW_Wisdom_FileName_Extension;// = "FFTW_Wisdom";

    /// FFTW plan for the Real-to-Complex transform
    fftwf_plan fftw_plan_R2C;
    /// FFTW plan for the Complex-to-Real transform
    fftwf_plan fftw_plan_C2R;

private:

};// TFFTWComplexMatrix

#endif	/* FFTWCOMPLEXMATRIX_H */

