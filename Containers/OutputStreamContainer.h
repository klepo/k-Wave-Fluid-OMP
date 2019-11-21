/**
 * @file      OutputStreamContainer.h
 *
 * @author    Jiri Jaros \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            jarosjir@fit.vutbr.cz
 *
 * @brief     The header file defining the output stream container.
 *
 * @version   kspaceFirstOrder 2.17
 *
 * @date      27 August    2017, 08:58 (created) \n
 *            20 February  2019, 14:45 (revised)
 *
 * @copyright Copyright (C) 2019 Jiri Jaros and Bradley Treeby.
 *
 * This file is part of the C++ extension of the [k-Wave Toolbox](http://www.k-wave.org).
 *
 * This file is part of the k-Wave. k-Wave is free software: you can redistribute it and/or modify it under the terms
 * of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with k-Wave.
 * If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).
 */

#ifndef OUTPUT_STREAM_CONTAINER_H
#define OUTPUT_STREAM_CONTAINER_H

#include <map>

#include <Containers/MatrixContainer.h>
#include <OutputStreams/BaseOutputStream.h>

#include <Utils/MatrixNames.h>
#include <Utils/DimensionSizes.h>


/**
 * @class   OutputStreamContainer
 * @brief   A container for output streams.
 * @details The output stream container maintains matrices used to sample data.
 * These may or may not require some scratch place or reuse temp matrices.
 */
class OutputStreamContainer
{
  public:
    /**
      * @enum    OutputStreamIdx
      * @brief   Output streams identifiers in k-Wave.
      * @details Output streams identifiers in k-Wave.
      */
    enum class OutputStreamIdx
    {
      /// Pressure time series.
      kPressureRaw,
      /// Compressed pressure time series.
      kPressureC,
      /// RMS of pressure over sensor mask.
      kPressureRms,
      /// Max of pressure over sensor mask.
      kPressureMax,
      /// Min of pressure over sensor mask.
      kPressureMin,
      /// Max of pressure over all domain.
      kPressureMaxAll,
      /// Min of pressure over all domain.
      kPressureMinAll,

      /// Velocity x time series.
      kVelocityXRaw,
      /// Velocity y time series.
      kVelocityYRaw,
      /// Velocity z time series.
      kVelocityZRaw,
      /// Compressed velocity x time series.
      kVelocityXC,
      /// Compressed velocity y time series.
      kVelocityYC,
      /// Compressed velocity z time series.
      kVelocityZC,
      /// Non staggered velocity x time series.
      kVelocityXNonStaggeredRaw,
      /// Non staggered velocity y time series.
      kVelocityYNonStaggeredRaw,
      /// Non staggered velocity z time series.
      kVelocityZNonStaggeredRaw,
      /// Compressed non staggered velocity x time series.
      kVelocityXNonStaggeredC,
      /// Compressed non staggered velocity y time series.
      kVelocityYNonStaggeredC,
      /// Compressed non staggered velocity z time series.
      kVelocityZNonStaggeredC,

      /// RMS of velocity x over sensor mask.
      kVelocityXRms,
      /// RMS of velocity y over sensor mask.
      kVelocityYRms,
      /// RMS of velocity z over sensor mask.
      kVelocityZRms,
      /// Max of velocity x over sensor mask.
      kVelocityXMax,
      /// Max of velocity y over sensor mask.
      kVelocityYMax,
      /// Max of velocity z over sensor mask.
      kVelocityZMax,
      /// Min of velocity x over sensor mask.
      kVelocityXMin,
      /// Min of velocity y over sensor mask.
      kVelocityYMin,
      /// Min of velocity z over sensor mask.
      kVelocityZMin,

      /// Max of velocity x over all domain.
      kVelocityXMaxAll,
      /// Max of velocity y over all domain.
      kVelocityYMaxAll,
      /// Max of velocity z over all domain.
      kVelocityZMaxAll,
      /// Min of velocity x over all domain.
      kVelocityXMinAll,
      /// Min of velocity y over all domain.
      kVelocityYMinAll,
      /// Min of velocity z over all domain.
      kVelocityZMinAll,

      /// Average intensity x over sensor mask
      kIntensityXAvg,
      /// Average intensity y over sensor mask
      kIntensityYAvg,
      /// Average intensity z over sensor mask
      kIntensityZAvg,
      /// Average intensity x over sensor mask using compression
      kIntensityXAvgC,
      /// Average intensity y over sensor mask using compression
      kIntensityYAvgC,
      /// Average intensity z over sensor mask using compression
      kIntensityZAvgC,

      /// Q term (volume rate of heat deposition)
      kQTerm,
      /// Q term (volume rate of heat deposition) using compression
      kQTermC,
    };// end of OutputStreamIdx


    /// Constructor.
    OutputStreamContainer();
    /// Copy constructor not allowed.
    OutputStreamContainer(const OutputStreamContainer&) = delete;
    /// Destructor.
    ~OutputStreamContainer();

    /// Operator = not allowed.
    OutputStreamContainer& operator=(OutputStreamContainer&) = delete;

    /**
     * @brief  Get size of the container.
     * @return the size of the container
     */
    size_t size() const
    {
      return mContainer.size();
    };

    /**
     * @brief   Is the container empty?
     * @return  true - If the container is empty.
     */
    bool empty() const
    {
      return mContainer.empty();
    };

    /**
     * @brief operator []
     * @param [in] outputStreamIdx - Id of the output stream.
     * @return An element of the container.
     */
    BaseOutputStream& operator[](const OutputStreamIdx outputStreamIdx)
    {
      return (* (mContainer[outputStreamIdx]));
    };

    /**
     * @brief Add all streams in simulation in the container, set all streams records here!
     *
     * Please note, the matrix container has to be populated before calling this routine.
     *
     * @param [in] matrixContainer - matrix container to link the steams with sampled matrices and sensor masks.
     */
    void init(MatrixContainer& matrixContainer);

    /// Create all streams - opens the datasets.
    void createStreams();
    /// Reopen streams after checkpoint file (datasets).
    void reopenStreams();

    /// Sample all streams.
    void sampleStreams();
    /// Post-process all streams and flush them to the file.
    void postProcessStreams();
    /// Post-process 2 all streams and flush them to the file.
    void postProcessStreams2();
    /// Checkpoint streams.
    void checkpointStreams();

    /// Close all streams.
    void closeStreams();

    /// Free all streams - destroy them.
    void freeStreams();

  protected:
    /**
     * @brief Create a new output stream.
     * @param [in] matrixContainer  - name of the HDF5 dataset or group
     * @param [in] sampledMatrixIdx - code id of the matrix
     * @param [in] fileObjectName   - name of the HDF5 dataset or group
     * @param [in] reduceOp         - reduction operator
     * @param [in] bufferToReuse    - buffer to reuse
     * @return New output stream with defined links.
     */
    BaseOutputStream* createOutputStream(MatrixContainer&                       matrixContainer,
                                         const MatrixContainer::MatrixIdx       sampledMatrixIdx,
                                         const MatrixName&                      fileObjectName,
                                         const BaseOutputStream::ReduceOperator reduceOp,
                                         float*                                 bufferToReuse = nullptr,
                                         bool                                   doNotSaveFlag = false);

  private:
    /// Map with output streams.
    std::map<OutputStreamIdx, BaseOutputStream*> mContainer;

}; // end of OutputStreamContainer
//----------------------------------------------------------------------------------------------------------------------

#endif	/* OUTPUT_STREAM_CONTAINER_H */
