/**
 * @file        OutputStreamContainer.h
 * @author      Jiri Jaros \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The header file defining the output stream container.
 *
 * @version     kspaceFirstOrder3D 2.16
 *
 * @date        27 August    2017, 08:58 (created) \n
 *              27 August    2017, 09:58 (revised)
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

#ifndef OUTPUT_STREAM_CONTAINER_H
#define OUTPUT_STREAM_CONTAINER_H

#include <map>

#include <Containers/MatrixContainer.h>
#include <OutputStreams/BaseOutputStream.h>

#include <Utils/MatrixNames.h>
#include <Utils/DimensionSizes.h>

/**
 * @class TOutputStreamContainer
 * @brief A container for output streams.
 * @details The output stream container maintains matrices used for sampling data.
 * These may or may not require some scratch place or reuse temp matrices.
 */
class TOutputStreamContainer
{
  public:
    /// Constructor.
    TOutputStreamContainer() {};
    /// Destructor.
    virtual ~TOutputStreamContainer();

    /**
     * @brief Get size of the container.
     * @details Get size of the container.
     */
    size_t size() const
    {
      return OutputStreamContainer.size();
    };

    /**
     * @brief  Is the container empty?
     * @details  Is the container empty?
     */
    bool empty() const
    {
      return OutputStreamContainer.empty();
    };

    /**
     * @brief Operator [].
     * @details Operator [].
     * @param [in] MatrixID
     * @return Ouptut stream
     */
    BaseOutputStream & operator [] (const TMatrixID MatrixID)
    {
      return (* (OutputStreamContainer[MatrixID]));
    };

    /// Create all streams in container (no file manipulation).
    void AddStreamsIntoContainer(TMatrixContainer & MatrixContainer);

    /// Create all streams - opens the datasets.
    void CreateStreams();
    /// Reopen streams after checkpoint file (datasets).
    void ReopenStreams();

    /// Sample all streams.
    void SampleStreams();
    /// Post-process all streams and flush them to the file.
    void PostProcessStreams();
    /// Checkpoint streams.
    void CheckpointStreams();

    /// Close all streams.
    void CloseStreams();

    /// Free all streams - destroy them.
    void FreeAllStreams();

  protected:
    /// Create a new output stream.
    BaseOutputStream * CreateNewOutputStream(TMatrixContainer & MatrixContainer,
                                                  const TMatrixID    SampledMatrixID,
                                                  const char *       HDF5_DatasetName,
                                                  const BaseOutputStream::ReduceOperator ReductionOp,
                                                  float *            BufferToReuse = NULL);

    /// Copy constructor not allowed for public.
    TOutputStreamContainer(const TOutputStreamContainer &);
    /// Operator = not allowed for public.
    TOutputStreamContainer & operator = (TOutputStreamContainer &);

  private:
    /// Output stream map.
    typedef map < TMatrixID, BaseOutputStream * > TOutputStreamMap;
    /// Map with output streams.
    TOutputStreamMap OutputStreamContainer;

}; // end of TOutputStreamContainer
//------------------------------------------------------------------------------

#endif	/* OUTPUT_STREAM_CONTAINER_H */