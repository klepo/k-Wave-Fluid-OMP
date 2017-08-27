/**
 * @file        OutputStreamContainer.cpp
 * @author      Jiri Jaros              \n
 *              Faculty of Information Technology\n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The implementation file for the output stream container.
 *
 * @version     kspaceFirstOrder3D 2.16
 * @date        27 August    2017, 08:59 (created) \n
 *              27 August    2017, 08:59 (revised)
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


#include <Parameters/Parameters.h>
#include <Containers/OutputStreamContainer.h>

#include <OutputStreams/BaseOutputStream.h>
#include <OutputStreams/IndexOutputStream.h>
#include <OutputStreams/CuboidOutputStream.h>
#include <OutputStreams/WholeDomainOutputStream.h>


//============================================================================//
//                        TOutputStreamContainer                              //
//============================================================================//

//----------------------------------------------------------------------------//
//--------------------------- Public methods ---------------------------------//
//----------------------------------------------------------------------------//


/**
 * Destructor
 */
TOutputStreamContainer::~TOutputStreamContainer()
{
  OutputStreamContainer.clear();
}// end of Destructor
//------------------------------------------------------------------------------

/**
 * Add all streams in simulation in the container, set all streams records here!
 * Please note, the Matrixcontainer has to be populated before calling this routine.
 *
 * @param [in] MatrixContainer - matrix container to link the steams with
 *                               sampled matrices and sensor masks
 */
void TOutputStreamContainer::AddStreamsIntoContainer(TMatrixContainer & MatrixContainer)
{

  Parameters& Params = Parameters::getInstance();

  using RO = BaseOutputStream::ReduceOperator;

  float * TempBufferX = MatrixContainer.GetMatrix<RealMatrix>(Temp_1_RS3D).getData();
  float * TempBufferY = MatrixContainer.GetMatrix<RealMatrix>(Temp_2_RS3D).getData();
  float * TempBufferZ = MatrixContainer.GetMatrix<RealMatrix>(Temp_3_RS3D).getData();

  //--------------------- Pressure ------------------/
  if (Params.getStorePressureRawFlag())
  {
    OutputStreamContainer[p_sensor_raw] = CreateNewOutputStream(MatrixContainer,
                                                                p,
                                                                kPressureRawName.c_str(),
                                                                RO::kNone,
                                                                TempBufferX);
  }// IsStore_p_raw

  if (Params.getStorePressureRmsFlag())
  {
    OutputStreamContainer[p_sensor_rms] = CreateNewOutputStream(MatrixContainer,
                                                                p,
                                                                kPressureRmsName.c_str(),
                                                                RO::kRms);
  }

  if (Params.getStorePressureMaxFlag())
  {
    OutputStreamContainer[p_sensor_max] = CreateNewOutputStream(MatrixContainer,
                                                                p,
                                                                kPressureMaxName.c_str(),
                                                                RO::kMax);
  }

  if (Params.getStorePressureMinFlag())
  {
    OutputStreamContainer[p_sensor_min] = CreateNewOutputStream(MatrixContainer,
                                                                p,
                                                                kPressureMinName.c_str(),
                                                                RO::kMin);
  }

  if (Params.getStorePressureMaxAllFlag())
  {
    OutputStreamContainer[p_sensor_max_all] =
            new WholeDomainOutputStream(Params.getOutputFile(),
                                             kPressureMaxAllName.c_str(),
                                             MatrixContainer.GetMatrix<RealMatrix>(p),
                                             RO::kMax);
  }

  if (Params.getStorePressureMinAllFlag())
  {
    OutputStreamContainer[p_sensor_min_all] =
            new WholeDomainOutputStream(Params.getOutputFile(),
                                             kPressureMinAllName.c_str(),
                                             MatrixContainer.GetMatrix<RealMatrix>(p),
                                             RO::kMin);
  }

  //--------------------- Velocity ------------------/
  if (Params.getStoreVelocityRawFlag())
  {
    OutputStreamContainer[ux_sensor_raw] = CreateNewOutputStream(MatrixContainer,
                                                                 ux_sgx,
                                                                 kUxName.c_str(),
                                                                 RO::kNone,
                                                                 TempBufferX);
    OutputStreamContainer[uy_sensor_raw] = CreateNewOutputStream(MatrixContainer,
                                                                 uy_sgy,
                                                                 kUyName.c_str(),
                                                                 RO::kNone,
                                                                 TempBufferY);
    OutputStreamContainer[uz_sensor_raw] = CreateNewOutputStream(MatrixContainer,
                                                                 uz_sgz,
                                                                 kUzName.c_str(),
                                                                 RO::kNone,
                                                                 TempBufferZ);
  }

  if (Params.getStoreVelocityNonStaggeredRawFlag())
  {
    OutputStreamContainer[ux_shifted_sensor_raw] = CreateNewOutputStream(MatrixContainer,
                                                                         ux_shifted,
                                                                         kUxNonStaggeredName.c_str(),
                                                                         RO::kNone,
                                                                         TempBufferX);
    OutputStreamContainer[uy_shifted_sensor_raw] = CreateNewOutputStream(MatrixContainer,
                                                                         uy_shifted,
                                                                         kUyNonStaggeredName.c_str(),
                                                                         RO::kNone,
                                                                         TempBufferY);
    OutputStreamContainer[uz_shifted_sensor_raw] = CreateNewOutputStream(MatrixContainer,
                                                                         uz_shifted,
                                                                         kUzNonStaggeredName.c_str(),
                                                                         RO::kNone,
                                                                         TempBufferZ);
  }

  if (Params.getStoreVelocityRmsFlag())
  {
    OutputStreamContainer[ux_sensor_rms] = CreateNewOutputStream(MatrixContainer,
                                                                 ux_sgx,
                                                                 kUxRmsName.c_str(),
                                                                 RO::kRms);
    OutputStreamContainer[uy_sensor_rms] = CreateNewOutputStream(MatrixContainer,
                                                                 uy_sgy,
                                                                 kUyRmsName.c_str(),
                                                                 RO::kRms);
    OutputStreamContainer[uz_sensor_rms] = CreateNewOutputStream(MatrixContainer,
                                                                 uz_sgz,
                                                                 kUzRmsName.c_str(),
                                                                 RO::kRms);
  }

  if (Params.getStoreVelocityMaxFlag())
  {
    OutputStreamContainer[ux_sensor_max] = CreateNewOutputStream(MatrixContainer,
                                                                 ux_sgx,
                                                                 kUxMaxName.c_str(),
                                                                 RO::kMax);
    OutputStreamContainer[uy_sensor_max] = CreateNewOutputStream(MatrixContainer,
                                                                 uy_sgy,
                                                                 kUyMaxName.c_str(),
                                                                 RO::kMax);
    OutputStreamContainer[uz_sensor_max] = CreateNewOutputStream(MatrixContainer,
                                                                 uz_sgz,
                                                                 kUzMaxName.c_str(),
                                                                 RO::kMax);
  }

  if (Params.getStoreVelocityMinFlag())
  {
    OutputStreamContainer[ux_sensor_min] = CreateNewOutputStream(MatrixContainer,
                                                                 ux_sgx,
                                                                 kUxMinName.c_str(),
                                                                 RO::kMin);
    OutputStreamContainer[uy_sensor_min] = CreateNewOutputStream(MatrixContainer,
                                                                 uy_sgy,
                                                                 kUyMinName.c_str(),
                                                                 RO::kMin);
    OutputStreamContainer[uz_sensor_min] = CreateNewOutputStream(MatrixContainer,
                                                                 uz_sgz,
                                                                 kUzMinName.c_str(),
                                                                 RO::kMin);
  }

  if (Params.getStoreVelocityMaxAllFlag())
  {
    OutputStreamContainer[ux_sensor_max_all] =
            new WholeDomainOutputStream(Params.getOutputFile(),
                                             kUxMaxAllName.c_str(),
                                             MatrixContainer.GetMatrix<RealMatrix>(ux_sgx),
                                             RO::kMax);
    OutputStreamContainer[uy_sensor_max_all] =
            new WholeDomainOutputStream(Params.getOutputFile(),
                                             kUyMaxAllName.c_str(),
                                             MatrixContainer.GetMatrix<RealMatrix>(uy_sgy),
                                             RO::kMax);
    OutputStreamContainer[uz_sensor_max_all] =
            new WholeDomainOutputStream(Params.getOutputFile(),
                                             kUzMaxAllName.c_str(),
                                             MatrixContainer.GetMatrix<RealMatrix>(uz_sgz),
                                             RO::kMax);
  }

  if (Params.getStoreVelocityMinAllFlag())
  {
    OutputStreamContainer[ux_sensor_min_all] =
            new WholeDomainOutputStream(Params.getOutputFile(),
                                             kUxMinAllName.c_str(),
                                             MatrixContainer.GetMatrix<RealMatrix>(ux_sgx),
                                             RO::kMin);
    OutputStreamContainer[uy_sensor_min_all] =
            new WholeDomainOutputStream(Params.getOutputFile(),
                                             kUyMinAllName.c_str(),
                                             MatrixContainer.GetMatrix<RealMatrix>(uy_sgy),
                                             RO::kMin);
    OutputStreamContainer[uz_sensor_min_all] =
            new WholeDomainOutputStream(Params.getOutputFile(),
                                             kUzMinAllName.c_str(),
                                             MatrixContainer.GetMatrix<RealMatrix>(uz_sgz),
                                             RO::kMin);
  }
}// end of AddStreamsdIntoContainer
//------------------------------------------------------------------------------

/**
 * Create all streams.
 */
void TOutputStreamContainer::CreateStreams()
{
  for (TOutputStreamMap::iterator it = OutputStreamContainer.begin(); it != OutputStreamContainer.end(); it++)
  {
    if (it->second)
    {
      (it->second)->create();
    }
  }
}// end of CreateStreams
//------------------------------------------------------------------------------

/**
 * Reopen all streams after restarting form checkpoint.
 */
void TOutputStreamContainer::ReopenStreams()
{
  for (TOutputStreamMap::iterator it = OutputStreamContainer.begin(); it != OutputStreamContainer.end(); it++)
  {
    if (it->second)
    {
      (it->second)->reopen();
    }
  }
}// end of ReopenStreams
//------------------------------------------------------------------------------


/**
 * Sample all streams.
 */
void TOutputStreamContainer::SampleStreams()
{
  for (TOutputStreamMap::iterator it = OutputStreamContainer.begin(); it != OutputStreamContainer.end(); it++)
  {
    if (it->second)
    {
      (it->second)->sample();
    }
  }
}// end of SampleStreams
//------------------------------------------------------------------------------


/**
 * Checkpoint streams without post-processing (flush to the file).
 */
void TOutputStreamContainer::CheckpointStreams()
{
  for (TOutputStreamMap::iterator it = OutputStreamContainer.begin(); it != OutputStreamContainer.end(); it++)
  {
    if (it->second)
    {
      (it->second)->checkpoint();
    }
  }
}// end of CheckpointStreams
//------------------------------------------------------------------------------

/**
 * /// Post-process all streams and flush them to the file.
 */
void TOutputStreamContainer::PostProcessStreams()
{
  for (TOutputStreamMap::iterator it = OutputStreamContainer.begin(); it != OutputStreamContainer.end(); it++)
  {
    if (it->second)
    {
      (it->second)->postProcess();
    }
  }
}// end of CheckpointStreams
//------------------------------------------------------------------------------


/**
 * Close all streams (apply post-processing if necessary, flush data and close).
 */
void TOutputStreamContainer::CloseStreams()
{
  for (TOutputStreamMap::iterator it = OutputStreamContainer.begin(); it != OutputStreamContainer.end(); it++)
  {
    if (it->second)
    {
      (it->second)->close();
    }
  }
}// end of CloseStreams
//------------------------------------------------------------------------------

/**
 *  Free all streams- destroy them.
 */
void TOutputStreamContainer::FreeAllStreams()
{
  for (TOutputStreamMap::iterator it = OutputStreamContainer.begin(); it != OutputStreamContainer.end(); it++)
  {
    if (it->second)
    {
      delete it->second;
    }
  }
  OutputStreamContainer.clear();
}// end of FreeAllStreams
//------------------------------------------------------------------------------


//----------------------------------------------------------------------------//
//--------------------------- Protected methods ------------------------------//
//----------------------------------------------------------------------------//


/**
 * Create a new output stream.
 * @param [in] MatrixContainer  - name of the HDF5 dataset or group
 * @param [in] SampledMatrixID  - code id of the matrix
 * @param [in] HDF5_DatasetName - name of the HDF5 dataset or group
 * @param [in] ReductionOp      - reduction operator
 * @param [in] BufferToReuse   - buffer to reuse
 * @return new output stream with defined links.
 */
BaseOutputStream * TOutputStreamContainer::CreateNewOutputStream(TMatrixContainer & MatrixContainer,
                                                                      const TMatrixID    SampledMatrixID,
                                                                      const char *       HDF5_DatasetName,
                                                                      const BaseOutputStream::ReduceOperator  ReductionOp,
                                                                      float *            BufferToReuse)
{
  Parameters& Params = Parameters::getInstance();

  BaseOutputStream * Stream = NULL;

  if (Params.getSensorMaskType() == Parameters::SensorMaskType::kIndex)
  {
    Stream = new IndexOutputStream(Params.getOutputFile(),
                                        HDF5_DatasetName,
                                        MatrixContainer.GetMatrix<RealMatrix>(SampledMatrixID),
                                        MatrixContainer.GetMatrix<IndexMatrix>(sensor_mask_index),
                                        ReductionOp,
                                        BufferToReuse);
  }
  else
  {
    Stream = new CuboidOutputStream(Params.getOutputFile(),
                                         HDF5_DatasetName,
                                         MatrixContainer.GetMatrix<RealMatrix>(SampledMatrixID),
                                         MatrixContainer.GetMatrix<IndexMatrix>(sensor_mask_corners),
                                         ReductionOp,
                                         BufferToReuse);
  }

  return Stream;
}// end of CreateNewOutputStream
//------------------------------------------------------------------------------

//----------------------------------------------------------------------------//
//--------------------------- Private methods --------------------------------//
//----------------------------------------------------------------------------//
