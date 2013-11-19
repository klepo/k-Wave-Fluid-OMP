/**
 * @file        TimeMeasure.h
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia   \n
 *              jiri.jaros@anu.edu.au
 * 
 * @brief       The header file containing the class measuring elapsed time
 * 
 * @version     kspaceFirstOrder3D 2.14
 * @date        15 August 2012, 9:35          (created) \n
 *              17 September 2012, 15:35      (revised)
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
 *  */


#ifndef TIMEMEASURE_H
#define	TIMEMEASURE_H

#include <exception>

#ifdef _OPENMP  
  #include <omp.h>
#else
#include <sys/time.h>
#endif


/**
 * @class  TTimeMesssure
 * @brief  Class measuring elapsed time
 */
class TTimeMesssure{
public :
    /// Start timestamp of the interval
    double StartTime;
    /// Stop timestamp of the interval
    double StopTime;
    
    /// Default constructor
    TTimeMesssure() : StartTime(0.0), StopTime(0.0) {};
    
    /**
     * @brief Copy constructor
     * @param [in] src - the other class to copy from
     */
    TTimeMesssure(const TTimeMesssure & src) :
        StartTime(src.StartTime),
        StopTime (src.StopTime) { };
    
    /**
     * @brief operator =
     * @param [in] src - source
     * @return 
     */
    TTimeMesssure& operator = (const TTimeMesssure & src){
        if (this != &src){
            StartTime = src.StartTime;
            StopTime  = src.StopTime;            
        }       
        return *this;
    };
    
    
    /// Destructor
    virtual ~TTimeMesssure() {};
    
    ///Get start timestamp 
    void Start() {
       #ifdef _OPENMP
         StartTime = omp_get_wtime(); 
         
       #else
          timeval ActTime;
          gettimeofday(&ActTime, NULL);
          StartTime = ActTime.tv_sec + ActTime.tv_usec * 1.0e-6;
       #endif
    };
    ///Get stop timestamp
    void Stop()  {
        #ifdef _OPENMP
          StopTime = omp_get_wtime();    
        #else
          timeval ActTime;
          gettimeofday(&ActTime, NULL);
          StopTime = ActTime.tv_sec + ActTime.tv_usec * 1.0e-6;
       #endif
    };
    
    /**
     * Get elapsed time
     * @return elapsed time between start timestamp and stop timestamp
     */
    double GetElapsedTime() const {return StopTime - StartTime;};
};// end of TTimeMesssure
//------------------------------------------------------------------------------

#endif	/* TIMEMEASURE_H */

