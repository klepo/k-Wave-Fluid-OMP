/**
 * @file        FftwComplexMatrix.h
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The header file containing the class that implements
 *              3D FFT using the FFTW interface
 *
 * @version     kspaceFirstOrder3D 2.16
 *
 * @date        09 August    2011, 13:10 (created) \n
 *              30 August    2017, 16:04 (revised)
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

#ifndef FFTW_COMPLEX_MATRIX_H
#define FFTW_COMPLEX_MATRIX_H


#include <fftw3.h>

#include <MatrixClasses/ComplexMatrix.h>

/**
 * @class FftwComplexMatrix
 * @brief   Class implementing 3D and 1D Real-To-Complex and Complex-To-Real transforms using FFTW interface.
 * @details Class implementing 3D and 1D Real-To-Complex and Complex-To-Real transforms using FFTW interface.
 *
 */
class FftwComplexMatrix : public ComplexMatrix
{
  public:
    /// Default constructor not allowed for public
    FftwComplexMatrix() = delete;
    /**
     * @brief Constructor, inherited from ComplexMatrix.
     * @param [in] dimensionSizes - Dimension sizes of the matrix.
     */
    FftwComplexMatrix(const DimensionSizes& dimensionSizes);
    /// Copy constructor not allowed for public.
    FftwComplexMatrix(const FftwComplexMatrix& src) = delete;

    /// Destructor.
    virtual ~FftwComplexMatrix();

    /// Operator = not allowed for public.
    FftwComplexMatrix & operator = (const FftwComplexMatrix& src);


    /**
     * @brief Create FFTW plan for 3D Real-to-Complex.
     * @param [in] inMatrix      - Input matrix serving as scratch place for planning.
     * @throw std::runtime_error - If the plan can't be created.
     *
     * @warning Unless FFTW_ESTIMATE flag is specified, the content of the inMatrix is destroyed!
     */
    void createR2CFftPlan3D(RealMatrix& inMatrix);
    /**
     * @brief Create FFTW plan for 3D Complex-to-Real.
     * @param [in] outMatrix     - Output matrix serving as scratch place for planning.
     * @throw std::runtime_error - If the plan can't be created.
     *
     * @warning Unless FFTW_ESTIMATE flag is specified, the content of the outMatrix is destroyed!
     */
    void createC2RFftPlan3D(RealMatrix& outMatrix);

    /**
     * @brief Create an FFTW plan for 1D Real-to-Complex in the x dimension.
     *
     * There are two versions of this routine for GCC+FFTW and ICPC + MKL, otherwise it will not build! \n
     * The FFTW version processes the whole matrix at one while the MKL slab by slab.
     *
     * @param [in,out] inMatrix  - RealMatrix of which to create the plan.
     * @throw std::runtime_error - If the plan can't be created.
     *
     * @warning Unless FFTW_ESTIMATE flag is specified, the content of the inMatrix is destroyed!
     */
    void createR2CFftPlan1DX(RealMatrix& inMatrix);
    /**
     * @brief Create an FFTW plan for 1D Real-to-Complex in the y dimension.
     *
     * There are two versions of this routine for GCC+FFTW and ICPC + MKL, otherwise it will not build! \n
     * The FFTW version processes the whole matrix at one while the MKL slab by slab
     *
     * @param   [in,out] inMatrix  - RealMatrix of which to create the plan.
     * @throw   std::runtime_error - If the plan can't be created.
     *
     * @warning Unless FFTW_ESTIMATE flag is specified, the content of the inMatrix is destroyed!
     * @warning The FFTW matrix must be able to store 2 * (nx * (ny /2 + 1) * nz) elements - possibly more than
     *          reduced dims!
     */
    void createR2CFftPlan1DY(RealMatrix& inMatrix);
    /**
     * @brief Create an FFTW plan for 1D Real-to-Complex in the z dimension.
     *
     * There are two versions of this routine for GCC+FFTW and ICPC + MKL, otherwise it will not build! \n
     * The FFTW version processes the whole matrix at one while the MKL slab by slab
     *
     * @param   [in,out] inMatrix  - RealMatrix of which to create the plan.
     * @throw   std::runtime_error - If the plan can't be created.
     *
     * @warning Unless FFTW_ESTIMATE flag is specified, the content of the inMatrix is destroyed!
     * @warning The FFTW matrix must be able to store 2 * (nx * by * (nz / 2 + 1)) elements - possibly more than
     *          reduced dims!
     */
    void createR2CFftPlan1DZ(RealMatrix& inMatrix);

    /**
     * @brief Create FFTW plan for Complex-to-Real in the x dimension.
     *
     * There are two versions of this routine for GCC+FFTW and ICPC + MKL, otherwise it will not build! \n
     * The FFTW version processes the whole matrix at one while the MKL slab by slab.
     *
     * @param [in, out] outMatrix - RealMatrix of which to create the plan.
     * @throw std::runtime_error  - If the plan can't be created.
     *
     * @warning Unless FFTW_ESTIMATE flag is specified, the content of the outMatrix is destroyed!
     */
    void createC2RFftPlan1DX(RealMatrix& outMatrix);
    /**
     * @brief Create FFTW plan for Complex-to-Real in the y dimension.
     *
     * There are two versions of this routine for GCC+FFTW and ICPC + MKL, otherwise it will not build! \n
     * The FFTW version processes the whole matrix at one while the MKL slab by slab.
     *
     * @param [in, out] outMatrix - RealMatrix of which to create the plan.
     * @throw std::runtime_error  - If the plan can't be created.
     *
     * @warning Unless FFTW_ESTIMATE flag is specified, the content of the outMatrix is destroyed!
     */
    void createC2RFftPlan1DY(RealMatrix& outMatrix);
    /**
     * @brief Create FFTW plan for Complex-to-Real in the z dimension.
     *
     * There are two versions of this routine for GCC+FFTW and ICPC + MKL, otherwise it will not build! \n
     * The FFTW version processes the whole matrix at one while the MKL slab by slab.
     *
     * @param [in, out] outMatrix - RealMatrix of which to create the plan.
     * @throw std::runtime_error  - If the plan can't be created.
     *
     * @warning Unless FFTW_ESTIMATE flag is specified, the content of the outMatrix is destroyed!
     */
    void createC2RFftPlan1DZ(RealMatrix& outMatrix);



    /**
     * @brief Compute forward out-of-place 3D Real-to-Complex FFT.
     *
     * @param [in] inMatrix      - Input data for the forward FFT.
     * @throw std::runtime_error - If the plan is not valid.
     */
    void computeR2CFft3D(RealMatrix& inMatrix);
    /**
     * @brief Compute forward out-of-place 3D Complex-to-Real FFT.
     *
     * @param [out] outMatrix    - Output of the inverse FFT.
     * @throw std::runtime_error - If the plan is not valid.
     */
    void computeC2RFft3D(RealMatrix& outMatrix);

    /**
     * @brief Compute 1D out-of-place Real-to-Complex FFT in the x dimension.
     *
     * There are two versions of this routine for GCC+FFTW and ICPC + MKL, otherwise it will not build!
     * The FFTW version processes the whole matrix at one while the MKL slab by slab.
     *
     * @param [in] inMatrix      - Input matrix
     * @throw std::runtime_error - If the plan is not valid.
     */
    void computeR2CFft1DX(RealMatrix& inMatrix);
    /**
     * @brief Compute 1D out-of-place Real-to-Complex FFT in the y dimension.
     *
     * There are two versions of this routine for GCC+FFTW and ICPC + MKL, otherwise it will not build!
     * The FFTW version processes the whole matrix at one while the MKL slab by slab.
     *
     * @param [in] inMatrix      - Input matrix
     * @throw std::runtime_error - If the plan is not valid.
     */
    void computeR2CFft1DY(RealMatrix& inMatrix);
    /**
     * @brief Compute 1D out-of-place Real-to-Complex FFT in the z dimension.
     *
     * There are two versions of this routine for GCC+FFTW and ICPC + MKL, otherwise it will not build!
     * The FFTW version processes the whole matrix at one while the MKL slab by slab.
     *
     * @param [in] inMatrix      - Input matrix
     * @throw std::runtime_error - If the plan is not valid.
     */
    void computeR2CFft1DZ(RealMatrix& inMatrix);

    /**
     * @brief Compute 1D out-of-place Complex-to-Real FFT in the x dimension.
     *
     * There are two versions of this routine for GCC+FFTW and ICPC + MKL, otherwise it will not build!
     * The FFTW version processes the whole matrix at one while the MKL slab by slab.
     *
     * @param [out] outMatrix    - Output matrix
     * @throw std::runtime_error - If the plan is not valid.
     */
    void computeC2RFft1DX(RealMatrix& outMatrix);
    /**
     * @brief Compute 1D out-of-place Complex-to-Real FFT in the y dimension.
     *
     * There are two versions of this routine for GCC+FFTW and ICPC + MKL, otherwise it will not build!
     * The FFTW version processes the whole matrix at one while the MKL slab by slab.
     *
     * @param [out] outMatrix    - Output matrix
     * @throw std::runtime_error - If the plan is not valid.
     */
    void computeC2RFft1DY(RealMatrix& outMatrix);
    /**
     * @brief Compute 1D out-of-place Complex-to-Real FFT in the z dimension.
     *
     * There are two versions of this routine for GCC+FFTW and ICPC + MKL, otherwise it will not build!
     * The FFTW version processes the whole matrix at one while the MKL slab by slab.
     *
     * @param [out] outMatrix      Output matrix
     * @throw std::runtime_error - If the plan is not valid.
     */
    void computeC2RFft1DZ(RealMatrix& outMatrix);

    /**
     * @brief   Export wisdom to the file.
     * @warning This routine is supported only while compiling with the GNU C++ compiler.
     */
    static void exportWisdom();
    /**
     * @brief   Import wisdom from the file.
     * @warning This routine is supported only while compiling with the GNU C++ compiler.
     */
    static void importWisdom();
    /// Destroy wisdom in the file (delete it).
    static void deleteStoredWisdom();

protected:

    /**
     * Get Wisdom file name (derive it form the checkpoint filename).
     * @return the filename for wisdom.
     */
    static std::string getWisdomFileName();

    /// Allocate memory for the FFTW matrix.
    virtual void allocateMemory();
    /// Free memory of the FFTW matrix.
    virtual void freeMemory();

    /// FFTW plan flag.
    static const unsigned kFftMeasureFlag  = FFTW_MEASURE;
    /// The file extension for FFTW wisdom.
    static const std::string kFftWisdomFileExtension;// = "FFTW_Wisdom";

    /// FFTW plan for the 3D Real-to-Complex transform.
    fftwf_plan mR2CFftPlan3D;
    /// FFTW plan for the 3D Complex-to-Real transform.
    fftwf_plan mC2RFftPlan3D;

    /// FFTW plan for the 1D Real-to-Complex transform in the x dimension.
    fftwf_plan mR2CFftPlan1DX;
    /// FFTW plan for the 3D Real-to-Complex transform in the y dimension.
    fftwf_plan mR2CFftPlan1DY;
    /// FFTW plan for the 3D Real-to-Complex transform in the z dimension.
    fftwf_plan mR2CFftPlan1DZ;

    /// FFTW plan for the 3D Complex-to-Real transform in the x dimension.
    fftwf_plan mC2RFftPlan1DX;
    /// FFTW plan for the 3D Complex-to-Real transform in the y dimension.
    fftwf_plan mC2RFftPlan1DY;
    /// FFTW plan for the 3Z Complex-to-Real transform in the z dimension.
    fftwf_plan mC2RFftPlan1DZ;

private:

};// end of FftwComplexMatrix
//----------------------------------------------------------------------------------------------------------------------
#endif	/* FFTW_COMPLEX_MATRIX_H */
