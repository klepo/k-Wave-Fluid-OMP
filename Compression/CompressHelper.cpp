/**
 * @file      CompressHelper.cpp
 *
 * @author    Petr Kleparnik \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            ikleparnik@fit.vutbr.cz
 *
 * @brief     The implementation file containing CompressHelper class definition.
 *
 * @version   kspaceFirstOrder 2.17
 *
 * @date      8  September 2016, 12:00 (created) \n
 *            2  April     2019, 15:11 (revised)
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

#include <Compression/CompressHelper.h>

// Initialization of the singleton instance flag
bool CompressHelper::sCompressHelperInstanceFlag = false;

// Initialization of the instance
CompressHelper* CompressHelper::sCompressHelperInstance = nullptr;

/**
 * @brief Initialization
 * @param[in] period Period
 * @param[in] mos Multiple of overlap size
 * @param[in] harmonics Number of harmonics
 * @param[in] normalize Normalizes basis functions for compression (optional)
 */
void CompressHelper::init(float period, hsize_t mos, hsize_t harmonics, bool normalize)
{
  mOSize = hsize_t(period * mos);
  mBSize = mOSize * 2 + 1;
  mPeriod = period;
  mMos = mos;
  mHarmonics = harmonics;
  mStride = harmonics * 2;
  mB = new float[mBSize]();
  mE = new floatC[harmonics * mBSize]();
  mBE = new floatC[harmonics * mBSize]();
  mBE_1 = new floatC[harmonics * mBSize]();

  generateFunctions(mBSize, mOSize, period, harmonics, mB, mE, mBE, mBE_1, normalize);
}

/**
 * @brief Destructor
 */
CompressHelper::~CompressHelper()
{
  if (mB)
  {
    delete[] mB;
  }
  mB = nullptr;

  if (mE)
  {
    delete[] mE;
  }
  mE = nullptr;

  if (mBE)
  {
    delete[] mBE;
  }
  mBE = nullptr;

  if (mBE_1)
  {
    delete[] mBE_1;
  }
  mBE_1 = nullptr;

  sCompressHelperInstanceFlag = false;
  if (sCompressHelperInstance)
  {
    delete sCompressHelperInstance;
  }
  sCompressHelperInstance = nullptr;
}

/**
 * @brief Get instance of singleton class
 * @return Instance of singleton class
 */
CompressHelper &CompressHelper::getInstance()
{
  if(!sCompressHelperInstanceFlag)
  {
    sCompressHelperInstance = new CompressHelper();
    sCompressHelperInstanceFlag = true;
    return *sCompressHelperInstance;
  }
  else
  {
    return *sCompressHelperInstance;
  }
}

/**
 * @brief Computes period from data
 * @param[in] dataSrc Source data
 * @param[in] length Source length
 * @return Period
 */
float CompressHelper::findPeriod(const float* dataSrc, hsize_t length)
{
  float* dataTmp = new float[length]();
  float* peaksTmp = new float[length]();
  float* locsTmp = new float[length]();
  hsize_t peaksCount;
  float period;

  //xcorr(dataSrc, dataSrc, dataTmp, length, length);
  findPeaks(dataSrc, locsTmp, peaksTmp, length, peaksCount);

  float* newLocsTmp = new float[peaksCount]();
  float* newLocs = new float[peaksCount]();

  // Find max peak
  float m = std::numeric_limits<float>::min();
  for (hsize_t i = 0; i < peaksCount; i++)
  {
    if (peaksTmp[i] > m)
    {
      m = peaksTmp[i];
    }
  }

  // Filter peaks under 0.5 * max
  hsize_t j = 0;
  for (hsize_t i = 0; i < peaksCount; i++)
  {
    if (peaksTmp[i] > 0.5f * m)
    {
      newLocsTmp[j] = locsTmp[i];
      j++;
    }
  }

  float* locs = new float[j - 1]();
  diff(newLocsTmp, locs, j);
  period = median(locs, j - 1);

  if (dataTmp)
  {
    delete[] dataTmp;
    dataTmp = nullptr;
  }
  if (peaksTmp)
  {
    delete[] peaksTmp;
    peaksTmp = nullptr;
  }
  if (locsTmp)
  {
    delete[] locsTmp;
    locsTmp = nullptr;
  }
  if (newLocsTmp)
  {
    delete[] newLocsTmp;
    newLocsTmp = nullptr;
  }
  if (newLocs)
  {
    delete[] newLocs;
    newLocs = nullptr;
  }
  if (locs)
  {
    delete[] locs;
    locs = nullptr;
  }
  return period;
}

/**
 * @brief Computes sample value for given step using coefficients (decompression)
 * @param[in] cC Current coefficients
 * @param[in] lC Last coefficients
 * @param[in] stepLocal Local step between coefficients
 * @return Step value
 */
float CompressHelper::computeTimeStep(const float* cC, const float* lC, hsize_t stepLocal) const
{
  float stepValue = 0.0f;
  hsize_t sH = stepLocal;
  for (hsize_t h = 0; h < mHarmonics; h++)
  {
    stepValue += real(conj(reinterpret_cast<const floatC*>(cC)[h]) * getBE()[sH]) +
                 real(conj(reinterpret_cast<const floatC*>(lC)[h]) * getBE_1()[sH]);
    sH += mBSize;
  }
  return stepValue;
}

/**
 * @brief Returns complex exponencial window basis
 * @return Complex exponencial window basis
 */
const floatC* CompressHelper::getBE() const
{
  return mBE;
}

/**
 * @brief Returns inverted complex exponencial window basis
 * @return Inverted complex exponencial window basis
 */
const floatC* CompressHelper::getBE_1() const
{
  return mBE_1;
}

/**
 * @brief Returns overlap size
 * @return Overlap size
 */
hsize_t CompressHelper::getOSize() const
{
  return mOSize;
}

/**
 * @brief Return base size
 * @return Base size
 */
hsize_t CompressHelper::getBSize() const
{
  return mBSize;
}

/**
 * @brief Returns period
 * @return Period
 */
float CompressHelper::getPeriod() const
{
  return mPeriod;
}

/**
 * @brief Returns multiple of overlap size
 * @return Multiple of overlap size
 */
hsize_t CompressHelper::getMos() const
{
  return mMos;
}

/**
 * @brief Return number of harmonics
 * @return Number of harmonics
 */
hsize_t CompressHelper::getHarmonics() const
{
  return mHarmonics;
}

/**
 * @brief Returns coefficients stride for one step (harmonics * 2;)
 * @return Coefficients stride for one step
 */
hsize_t CompressHelper::getStride() const
{
  return mStride;
}

/**
 * @brief Constructor
 */
CompressHelper::CompressHelper()
{

}

/**
 * @brief Cross correlation
 * @param[in] dataSrc1 Source data 1
 * @param[in] dataSrc2 Source data 2
 * @param[out] dataDst Destination
 * @param[in] lengthSrc1 Source length 1
 * @param[in] lengthSrc2 Source length 2
 */
void CompressHelper::xcorr(const float* dataSrc1, const float* dataSrc2, float* dataDst, hsize_t lengthSrc1, hsize_t lengthSrc2)
{
  hsize_t i, j;
  hssize_t i1;
  float tmp;

  for (i = 0; i < lengthSrc1 + lengthSrc2 - 1; i++)
  {
    dataDst[i] = 0.0f;
    i1 = hssize_t(i);
    tmp = 0.0f;
    for (j = 0; j < lengthSrc2; j++)
    {
      if (i1 >= 0 && i1 < hssize_t(lengthSrc1))
      {
        tmp = tmp + (dataSrc1[i1] * dataSrc2[lengthSrc2 - 1 - j]);
      }
      i1 = i1 - 1;
      dataDst[i] = tmp;
    }
  }
}

/**
 * @brief Convolution
 * @param[in] dataSrc1 Source data 1
 * @param[in] dataSrc2 Source data 2
 * @param[out] dataDst Destination
 * @param[in] lengthSrc1 Source length 1
 * @param[in] lengthSrc2 Source length 2
 */
void CompressHelper::conv(const float* dataSrc1, const float* dataSrc2, float* dataDst, hsize_t lengthSrc1, hsize_t lengthSrc2)
{
  hsize_t i, j;
  hssize_t i1;
  float tmp;

  for (i = 0; i < lengthSrc1 + lengthSrc2 - 1; i++)
  {
    dataDst[i] = 0.0f;
    i1 = hssize_t(i);
    tmp = 0.0f;
    for (j = 0; j < lengthSrc2; j++)
    {
      if (i1 >= 0 && i1 < hssize_t(lengthSrc1))
      {
        tmp = tmp + (dataSrc1[i1] * dataSrc2[j]);
      }
      i1 = i1 - 1;
      dataDst[i] = tmp;
    }
  }
}

/**
 * @brief Finds peaks in source signal
 * @param[in] dataSrc Source data
 * @param[out] locsDst Destination for locations
 * @param[out] peaksDst Destination for peaks
 * @param[in] lengthSrc Source length
 * @param[out] lengthDst Destination length
 */
void CompressHelper::findPeaks(const float* dataSrc, float* locsDst, float* peaksDst, hsize_t lengthSrc, hsize_t &lengthDst)
{
  for (hsize_t i = 0; i < lengthSrc; i++)
  {
    locsDst[i] = 0.0f;
    peaksDst[i] = 0.0f;
  }

  lengthDst = 0;

  for (hsize_t i = 0; i < lengthSrc; i++)
  {
    // Peaks between
    if (i > 0 && i < lengthSrc - 1 && lengthSrc > 2 && dataSrc[i] > dataSrc[i - 1] && dataSrc[i] >= dataSrc[i + 1])
    {
      float d1 = dataSrc[i] - dataSrc[i - 1];
      float d2 = dataSrc[i] - dataSrc[i + 1];
      locsDst[lengthDst] = i + d1 / (d1 + d2) - 0.5f;
      peaksDst[lengthDst] = dataSrc[i];
      lengthDst++;
    }
  }
}

/**
 * @brief Differences
 * @param[in] dataSrc Source data
 * @param[out] dataDst Destination
 * @param[in] length Source length
 */
void CompressHelper::diff(const float* dataSrc, float* dataDst, hsize_t length)
{
  for (hsize_t i = 0; i < length - 1; i++)
  {
    dataDst[i] = dataSrc[i + 1] - dataSrc[i];
  }
}

/**
 * @brief Differences
 * @param[in] dataSrc Source data
 * @param[out] dataDst Destination
 * @param[in] length Source length
 */
void CompressHelper::diff(const hsize_t* dataSrc, hsize_t* dataDst, hsize_t length)
{
  for (hsize_t i = 0; i < length - 1; i++)
  {
    dataDst[i] = dataSrc[i + 1] - dataSrc[i];
  }
}

/**
 * @brief Mean
 * @param[in] dataSrc Source data
 * @param[in] length Source length
 * @return Mean
 */
float CompressHelper::mean(const float* dataSrc, hsize_t length)
{
  float sum = 0.0f;
  for (hsize_t i = 0; i < length; i++)
  {
    sum += dataSrc[i];
  }
  return sum / length;
}

/**
 * @brief Mean
 * @param[in] dataSrc Source data
 * @param[in] length Source length
 * @return Mean
 */
hsize_t CompressHelper::mean(const hsize_t* dataSrc, hsize_t length)
{
  float sum = 0.0f;
  for (hsize_t i = 0; i < length; i++)
  {
    sum += dataSrc[i];
  }
  return hsize_t(round(sum / length));
}

/**
 * @brief Median
 * @param[in] dataSrc Source data
 * @param[in] length Source length
 * @return Median
 */
float CompressHelper::median(const float* dataSrc, hsize_t length)
{
  std::vector<float> dataSrcVector(dataSrc, dataSrc + length);
  std::sort(dataSrcVector.begin(), dataSrcVector.end());
  return dataSrcVector[size_t(length / 2)];
}

/**
 * @brief Median
 * @param[in] dataSrc Source data
 * @param[in] length Source length
 * @return Median
 */
hsize_t CompressHelper::median(const hsize_t* dataSrc, hsize_t length)
{
  std::vector<hsize_t> dataSrcVector(dataSrc, dataSrc + length);
  std::sort(dataSrcVector.begin(), dataSrcVector.end());
  return dataSrcVector[size_t(length / 2)];
}

/**
 * @brief Generates basis function values
 * @param[in] bSize Base size
 * @param[in] oSize Overlap size
 * @param[in] period Period
 * @param[in] harmonics Harmonics
 * @param[out] b Window basis
 * @param[out] e Exponencial basis
 * @param[out] bE Complex exponencial window basis
 * @param[out] bE_1 Inverted complex exponencial window basis
 * @param[in] normalize Normalization flag (optional)
 */
void CompressHelper::generateFunctions(hsize_t bSize, hsize_t oSize, float period, hsize_t harmonics, float* b, floatC* e, floatC* bE, floatC* bE_1, bool normalize) const
{
  // Generate basis function (window)
  triangular(oSize, b);  // Triangular window
  //hann(oSize, b);        // Hann window

  // Generate complex exponential window basis functions
  for (hsize_t ih = 0; ih < harmonics; ih++)
  {
    generateE(period, ih, ih + 1, bSize, e);
    generateBE(ih, bSize, oSize, b, e, bE, bE_1, normalize);
  }
}

/**
 * @brief Generates triangular window
 * @param[in] oSize Overlap size
 * @param[out] w Window
 */
void CompressHelper::triangular(hsize_t oSize, float* w) const
{
  for (hsize_t x = 0; x < oSize; x++)
  {
    w[x] = float(x) / oSize;
  }
  for (hsize_t x = oSize; x < 2 * oSize + 1; x++)
  {
    w[x] = 2.0f - float(x) / oSize;
  }
}

/**
 * @brief Generates Hann window
 * @param[in] oSize Overlap size
 * @param[out] w Window
 */
void CompressHelper::hann(hsize_t oSize, float* w) const
{
  for (hsize_t x = 0; x < 2 * oSize + 1; x++)
  {
    w[x] = float(pow(sin(M_PI * x / (2.0f * oSize)), 2));
  }
}

/**
 * @brief Generates complex exponencial basis
 * @param[in] period Period
 * @param[in] ih Harmonics index
 * @param[in] h Harmonic multiple
 * @param[in] bSize Base size
 * @param[out] e Exponencial basis
 */
void CompressHelper::generateE(float period, hsize_t ih, hsize_t h, hsize_t bSize, floatC* e) const
{
  floatC i(0.0f, -1.0f);
  for (hsize_t x = 0; x < bSize; x++)
  {
    hsize_t hx = ih * bSize + x;
    e[hx] = std::exp(i * (2.0f * float(M_PI) / (period / float(h))) * float(x));
  }
}

/**
 * @brief Generates complex exponencial window basis
 * @param[in] ih Harmonics index
 * @param[in] bSize Base size
 * @param[in] oSize Overlap size
 * @param[in] b Window basis
 * @param[in] e Exponencial basis
 * @param[out] bE Complex exponencial window basis
 * @param[out] bE_1 Inverted complex exponencial window basis
 * @param[in] normalize Normalization flag (optional)
 */
void CompressHelper::generateBE(hsize_t ih, hsize_t bSize, hsize_t oSize, const float* b, const floatC* e, floatC* bE, floatC* bE_1, bool normalize) const
{
  for (hsize_t x = 0; x < hsize_t(bSize); x++)
  {
    hsize_t hx = ih * bSize + x;
    bE[hx] = b[x] * e[hx];
    bE_1[hx] = b[(x + oSize) % (bSize - 1)] * e[ih * bSize + ((x + oSize) % (bSize - 1))];
    if (normalize)
    {
      bE[hx] *= (2.0f / float(oSize));
      bE_1[hx] *= (2.0f / float(oSize));
    }
  }
}
