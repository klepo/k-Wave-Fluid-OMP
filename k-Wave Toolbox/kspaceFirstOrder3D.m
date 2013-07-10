function [sensor_data, mem_usage] = kspaceFirstOrder3D(kgrid, medium, source, sensor, varargin)
%KSPACEFIRSTORDER3D     3D time-domain simulation of wave propagation.
%
% DESCRIPTION:
%       kspaceFirstOrder3D simulates the time-domain propagation of
%       compressional waves through a three-dimensional homogeneous or
%       heterogeneous acoustic medium given four input structures: kgrid,
%       medium, source, and sensor. The computation is based on a
%       first-order k-space model which accounts for power law absorption
%       and a heterogeneous sound speed and density. If medium.BonA is
%       specified, cumulative nonlinear effects are also modelled. At each
%       time-step (defined by kgrid.t_array), the acoustic field parameters
%       at the positions defined by sensor.mask are recorded and stored. If
%       kgrid.t_array is set to 'auto', this array is automatically
%       generated using makeTime. An anisotropic absorbing boundary layer
%       called a perfectly matched layer (PML) is implemented to prevent
%       waves that leave one side of the domain being reintroduced from the
%       opposite side (a consequence of using the FFT to compute the
%       spatial derivatives in the wave equation). This allows infinite
%       domain simulations to be computed using small computational grids. 
%
%       For a homogeneous medium the formulation is exact and the
%       time-steps are only limited by the effectiveness of the perfectly
%       matched layer. For a heterogeneous medium, the solution represents
%       a leap-frog pseudospectral method with a k-space correction that
%       improves the accuracy of computing the temporal derivatives. This
%       allows larger time-steps to be taken for the same level of accuracy
%       compared to conventional pseudospectral time-domain methods. The
%       computational grids are staggered both spatially and temporally. 
%
%       An initial pressure distribution can be specified by assigning a
%       matrix (the same size as the computational grid) of arbitrary
%       numeric values to source.p0. A time varying pressure source can
%       similarly be specified by assigning a binary matrix (i.e., a matrix
%       of 1's and 0's with the same dimensions as the computational grid)
%       to source.p_mask where the 1's represent the grid points that form
%       part of the source. The time varying input signals are then
%       assigned to source.p. This must be the same length as kgrid.t_array
%       and can be a single time series (in which case it is applied to all
%       source elements), or a matrix of time series following the source
%       elements using MATLAB's standard column-wise linear matrix index
%       ordering. A time varying velocity source can be specified in an
%       analogous fashion, where the source location is specified by
%       source.u_mask, and the time varying input velocity is assigned to
%       source.ux, source.uy, and source.uz.
%
%       The field values are returned as arrays of time series at the
%       sensor locations defined by sensor.mask. This can be given either
%       as a binary matrix (i.e., a matrix of 1's and 0's with the same
%       dimensions as the computational grid) representing the grid points
%       within the computational grid that will collect the data, or as a
%       series of arbitrary Cartesian coordinates within the grid at which
%       the pressure values are calculated at each time step via
%       interpolation. The Cartesian points must be given as a 3 by N
%       matrix corresponding to the x, y, and z positions, respectively. If
%       no output is required, the sensor input can be replaced with an
%       empty array []. Both the source and sensor inputs can also be
%       replaced by an object of the kWaveTransducer class created using
%       makeTransducer. If sensor.mask is given as a set of Cartesian
%       coordinates, the computed sensor_data is returned in the same
%       order. If sensor.mask is given as a binary matrix, sensor_data is
%       returned using MATLAB's standard column-wise linear matrix index
%       ordering. In both cases, the recorded data is indexed as
%       sensor_data(sensor_point_index, time_index). For a binary sensor
%       mask, the field values at a particular time can be restored to the
%       sensor positions within the computation grid using
%       unmaskSensorData.  
%
%       By default, the recorded acoustic pressure field is passed directly
%       to the output sensor_data. However, other acoustic parameters can
%       also be recorded by setting sensor.record to a cell array of the
%       form {'p', 'u', 'p_max', ...}. For example, both the particle
%       velocity and the acoustic pressure can be return by setting
%       sensor.record = {'p', 'u'}. If sensor.record is given, the output
%       sensor_data is returned as a structure with the different outputs
%       appended as structure fields. For example, if sensor.record = {'p',
%       'p_final', 'p_max', 'u'}, the output would contain fields
%       sensor_data.p, sensor_data.p_final, sensor_data.p_max,
%       sensor_data.ux, sensor_data.uy, and sensor_data.uz. Most of the
%       output parameters are recorded at the given sensor positions and
%       are indexed as sensor_data.field(sensor_point_index, time_index).
%       The exceptions are the averaged quantities ('p_max', 'p_rms',
%       'u_max', 'p_rms', 'I_avg'), the 'all' quantities ('p_max_all',
%       'p_min_all', 'u_max_all', 'u_min_all'), and the final quantities
%       ('p_final', 'u_final'). The averaged quantities are indexed as
%       sensor_data.p_max(sensor_point_index), while the final and 'all'
%       quantities are returned over the entire grid and are indexed as
%       sensor_data.p_final(nx, ny, nz).
%
%       kspaceFirstOrder3D may also be used for time reversal image
%       reconstruction by assigning the time varying pressure recorded over
%       an arbitrary sensor surface to the input field
%       sensor.time_reversal_boundary_data. This data is then enforced in
%       time reversed order as a time varying Dirichlet boundary condition
%       over the sensor surface given by sensor.mask. The boundary data
%       must be indexed as
%       sensor.time_reversal_boundary_data(sensor_point_index, time_index).
%       If sensor.mask is given as a set of Cartesian coordinates, the
%       boundary data must be given in the same order. An equivalent binary
%       sensor mask (computed using nearest neighbour interpolation) is
%       then used to place the pressure values into the computational grid
%       at each time step. If sensor.mask is given as a binary matrix of
%       sensor points, the boundary data must be ordered using MATLAB's
%       standard column-wise linear matrix indexing. If no additional
%       inputs are required, the source input can be replaced with an empty
%       array []. 
%
%       Acoustic attenuation compensation can also be included during time
%       reversal image reconstruction by assigning the absorption
%       parameters medium.alpha_coeff and medium.alpha_power and reversing
%       the sign of the absorption term by setting medium.alpha_sign = [-1,
%       1]. This forces the propagating waves to grow according to the
%       absorption parameters instead of decay. The reconstruction should
%       then be regularised by assigning a filter to medium.alpha_filter
%       (this can be created using getAlphaFilter).
%
%       Note: To run a simple reconstruction example using time reversal
%       (that commits the 'inverse crime' of using the same numerical
%       parameters and model for data simulation and image reconstruction),
%       the sensor_data returned from a k-Wave simulation can be passed
%       directly to sensor.time_reversal_boundary_data with the input
%       fields source.p0 and source.p removed or set to zero. 
%
% USAGE:
%       sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor)
%       sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...) 
%
% INPUTS:
% The minimum fields that must be assigned to run an initial value problem
% (for example, a photoacoustic forward simulation) are marked with a *. 
%
%       kgrid*              - k-Wave grid structure returned by makeGrid
%                             containing Cartesian and k-space grid fields  
%       kgrid.t_array*      - evenly spaced array of time values [s] (set
%                             to 'auto' by makeGrid) 
%
%
%       medium.sound_speed* - sound speed distribution within the acoustic
%                             medium [m/s] 
%       medium.sound_speed_ref - reference sound speed used within the
%                             k-space operator (phase correction term)
%                             [m/s]
%       medium.density*     - density distribution within the acoustic
%                             medium [kg/m^3] 
%       medium.BonA         - parameter of nonlinearity
%       medium.alpha_power  - power law absorption exponent
%       medium.alpha_coeff  - power law absorption coefficient 
%                             [dB/(MHz^y cm)] 
%       medium.alpha_mode   - optional input to force either the absorption
%                             or dispersion terms in the equation of state
%                             to be excluded; valid inputs are
%                             'no_absorption' or 'no_dispersion' 
%       medium.alpha_filter - frequency domain filter applied to the
%                             absorption and dispersion terms in the
%                             equation of state 
%       medium.alpha_sign   - two element array used to control the sign of
%                             absorption and dispersion terms in the
%                             equation of state  
%
%
%       source.p0*          - initial pressure within the acoustic medium
%       source.p            - time varying pressure at each of the source
%                             positions given by source.p_mask 
%       source.p_mask       - binary matrix specifying the positions of the
%                             time varying pressure source distribution
%       source.p_mode       - optional input to control whether the input
%                             pressure is injected as a mass source or
%                             enforced as a dirichlet boundary condition;
%                             valid inputs are 'additive' (the default) or
%                             'dirichlet'    
%       source.ux           - time varying particle velocity in the
%                             x-direction at each of the source positions
%                             given by source.u_mask 
%       source.uy           - time varying particle velocity in the
%                             y-direction at each of the source positions
%                             given by source.u_mask 
%       source.uz           - time varying particle velocity in the
%                             z-direction at each of the source positions
%                             given by source.u_mask  
%       source.u_mask       - binary matrix specifying the positions of the
%                             time varying particle velocity distribution 
%       source.u_mode       - optional input to control whether the input
%                             velocity is applied as a force source or
%                             enforced as a dirichlet boundary condition;
%                             valid inputs are 'additive' (the default) or
%                             'dirichlet'
%
%
%       sensor.mask*        - binary matrix or a set of Cartesian points
%                             where the pressure is recorded at each
%                             time-step  
%       sensor.record       - cell array of the acoustic parameters to
%                             record in the form sensor.record = {'p', 'u',
%                             ...}; valid inputs are:  
%                               'p' (acoustic pressure)
%                               'p_max' (maximum pressure)
%                               'p_min' (minimum pressure)
%                               'p_rms' (RMS pressure)
%                               'p_final' (final pressure field at all grid points)
%                               'p_max_all' (maximum pressure at all grid points)
%                               'p_min_all' (minimum pressure at all grid points)
%                               'u' (particle velocity)
%                               'u_max' (maximum particle velocity)
%                               'u_min' (minimum particle velocity)
%                               'u_rms' (RMS particle velocity)
%                               'u_final' (final particle velocity field at all grid points)
%                               'u_max_all' (maximum particle velocity at all grid points)
%                               'u_min_all' (minimum particle velocity at all grid points)
%                               'I' (time varying acoustic intensity)
%                               'I_avg' (average acoustic intensity) 
%       sensor.record_start_index - time index at which the sensor should
%                             start recording the data specified by
%                             sensor.record (default = 1) 
%       sensor.time_reversal_boundary_data - time varying pressure
%                             enforced as a Dirichlet boundary condition
%                             over sensor.mask  
%       sensor.frequency_response - two element array specifying the center
%                             frequency and percentage bandwidth of a
%                             frequency domain Gaussian filter applied to
%                             the sensor_data
%
% Note: For heterogeneous medium parameters, medium.sound_speed and
% medium.density must be given in matrix form with the same dimensions as
% kgrid. For homogeneous medium parameters, these can be given as single
% numeric values. If the medium is homogeneous and velocity inputs or
% outputs are not required, it is not necessary to specify medium.density.
%
% OPTIONAL INPUTS:
%       Optional 'string', value pairs that may be used to modify the
%       default computational settings.
%
%       'CartInterp'- Interpolation mode used to extract the pressure when
%                     a Cartesian sensor mask is given. If set to 'nearest'
%                     and more than one Cartesian point maps to the same
%                     grid point, duplicated data points are discarded and
%                     sensor_data will be returned with less points than
%                     that specified by sensor.mask (default = 'nearest').
%       'CreateLog' - Boolean controlling whether the command line output
%                     is saved using the diary function with a date and
%                     time stamped filename (default = false). 
%       'DataCast'  - String input of the data type that variables are cast
%                     to before computation. For example, setting to
%                     'single' will speed up the computation time (due to
%                     the improved efficiency of fftn and ifftn for this
%                     data type) at the expense of a loss in precision.
%                     This variable is also useful for utilising GPU
%                     parallelisation through libraries such as GPUmat or
%                     AccelerEyesJacket by setting 'DataCast' to
%                     'GPUsingle' or 'gsingle' (default = 'off').
%       'DataRecast'- Boolean controlling whether the output data is cast
%                     back to double precision. If set to false,
%                     sensor_data will be returned in the data format set
%                     using the 'DataCast' option.
%       'DisplayMask' - Binary matrix overlayed onto the animated
%                     simulation display. Elements set to 1 within the
%                     display mask are set to black within the display
%                     (default = sensor.mask).
%       'LogScale'  - Boolean controlling whether the pressure field is log
%                     compressed before display (default = false). The data
%                     is compressed by scaling both the positive and
%                     negative values between 0 and 1 (truncating the data
%                     to the given plot scale), adding a scalar value
%                     (compression factor) and then using the corresponding
%                     portion of a log10 plot for the compression (the
%                     negative parts are remapped to be negative thus the
%                     default color scale will appear unchanged). The
%                     amount of compression can be controlled by adjusting
%                     the compression factor which can be given in place of
%                     the Boolean input. The closer the compression factor
%                     is to zero, the steeper the corresponding part of the
%                     log10 plot used, and the greater the compression (the
%                     default compression factor is 0.02).
%       'MovieArgs' - Settings for movie2avi. Parameters must be given as
%                     {param, value, ...} pairs within a cell array
%                     (default = {}).
%       'MovieName' - Name of the movie produced when 'RecordMovie' is set
%                     to true (default = 'date-time-kspaceFirstOrder2D').
%       'PlotFreq'  - The number of iterations which must pass before the
%                     simulation plot is updated (default = 10).
%       'PlotLayout'- Boolean controlling whether three plots are produced
%                     of the initial simulation layout (initial pressure,
%                     sound speed, density) (default = false).
%       'PlotPML'   - Boolean controlling whether the perfectly matched
%                     layer is shown in the simulation plots. If set to
%                     false, the PML is not displayed (default = true).
%       'PlotScale' - [min, max] values used to control the scaling for
%                     imagesc (visualisation). If set to 'auto', a
%                     symmetric plot scale is chosen automatically for each
%                     plot frame.
%       'PlotSim'   - Boolean controlling whether the simulation iterations
%                     are progressively plotted (default = true).
%       'PMLAlpha'  - Absorption within the perfectly matched layer in
%                     Nepers per grid point (default = 2).
%       'PMLInside' - Boolean controlling whether the perfectly matched
%                     layer is inside or outside the grid. If set to false,
%                     the input grids are enlarged by PMLSize before
%                     running the simulation (default = true). 
%       'PMLSize'   - Size of the perfectly matched layer in grid points.
%                     By default, the PML is added evenly to all sides of
%                     the grid, however, both PMLSize and PMLAlpha can be
%                     given as three element arrays to specify the x, y,
%                     and z properties, respectively. To remove the PML,
%                     set the appropriate PMLAlpha to zero rather than
%                     forcing the PML to be of zero size (default = 10).
%       'RecordMovie' - Boolean controlling whether the displayed image
%                     frames are captured and stored as a movie using
%                     movie2avi (default = false).  
%       'Smooth'    - Boolean controlling whether source.p0,
%                     medium.sound_speed, and medium.density are smoothed
%                     using smooth before computation. 'Smooth' can either
%                     be given as a single Boolean value or as a 3 element
%                     array to control the smoothing of source.p0,
%                     medium.sound_speed, and medium.density,
%                     independently.  
%       'SaveToDisk'- String containing a filename (including pathname if
%                     required). If set, after the precomputation phase,
%                     the input variables used in the time loop are saved
%                     the specified location in HDF5 format. The simulation
%                     then exits. The saved variables can be used to run
%                     simulations using the C++ code.
%       'StreamToDisk' - Boolean controlling whether sensor_data is
%                     periodically saved to disk to avoid storing the
%                     complete matrix in memory. StreamToDisk may also be
%                     given as an integer which specifies the number of
%                     times steps that are taken before the data is saved
%                     to disk (default = 200).
%
% OUTPUTS:
% If sensor.record is not defined by the user:
%       sensor_data - time varying pressure recorded at the sensor
%                     positions given by sensor.mask
%
% If sensor.record is defined by the user:
%       sensor_data.p         - time varying pressure recorded at the
%                               sensor positions given by sensor.mask
%                               (returned if 'p' is set)  
%       sensor_data.p_max     - maximum pressure recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'p_max' is set)  
%       sensor_data.p_min     - minimum pressure recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'p_min' is set)  
%       sensor_data.p_rms     - rms of the time varying pressure recorded
%                               at the sensor positions given by
%                               sensor.mask (returned if 'p_rms' is set)  
%       sensor_data.p_final   - final pressure field at all grid points
%                               within the domain (returned if 'p_final' is
%                               set)
%       sensor_data.p_max_all - maximum pressure recorded at all grid
%                               points within the domain (returned if
%                               'p_max_all' is set) 
%       sensor_data.p_min_all - minimum pressure recorded at all grid
%                               points within the domain (returned if
%                               'p_min_all' is set)  
%       sensor_data.ux        - time varying particle velocity in the
%                               x-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'u' is set) 
%       sensor_data.uy        - time varying particle velocity in the
%                               y-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'u' is set)   
%       sensor_data.uz        - time varying particle velocity in the
%                               z-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'u' is set)   
%       sensor_data.ux_max    - maximum particle velocity in the
%                               x-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'u_max' is set)   
%       sensor_data.uy_max    - maximum particle velocity in the
%                               y-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'u_max' is set)   
%       sensor_data.uz_max    - maximum particle velocity in the
%                               z-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'u_max' is set)  
%       sensor_data.ux_min    - minimum particle velocity in the
%                               x-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'u_min' is set) 
%       sensor_data.uy_min    - minimum particle velocity in the
%                               y-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'u_min' is set)   
%       sensor_data.uz_min    - minimum particle velocity in the
%                               z-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'u_min' is set) 
%       sensor_data.ux_rms    - rms of the time varying particle velocity
%                               in the x-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'u_rms' is set)   
%       sensor_data.uy_rms    - rms of the time varying particle velocity
%                               in the y-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'u_rms' is set)   
%       sensor_data.uz_rms    - rms of the time varying particle velocity
%                               in the z-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'u_rms' is set)   
%       sensor_data.ux_final  - final particle velocity field in the
%                               x-direction at all grid points within the
%                               domain (returned if 'u_final' is set) 
%       sensor_data.uy_final  - final particle velocity field in the
%                               y-direction at all grid points within the
%                               domain (returned if 'u_final' is set) 
%       sensor_data.uz_final  - final particle velocity field in the
%                               z-direction at all grid points within the
%                               domain (returned if 'u_final' is set) 
%       sensor_data.ux_max_all- maximum particle velocity in the
%                               x-direction recorded at all grid points
%                               within the domain (returned if 'u_max_all'
%                               is set)   
%       sensor_data.uy_max_all- maximum particle velocity in the
%                               y-direction recorded at all grid points
%                               within the domain (returned if 'u_max_all'
%                               is set)   
%       sensor_data.uz_max_all- maximum particle velocity in the
%                               z-direction recorded at all grid points
%                               within the domain (returned if 'u_max_all'
%                               is set)  
%       sensor_data.ux_min_all- minimum particle velocity in the
%                               x-direction recorded at all grid points
%                               within the domain (returned if 'u_min_all'
%                               is set)   
%       sensor_data.uy_min_all- minimum particle velocity in the
%                               y-direction recorded at all grid points
%                               within the domain (returned if 'u_min_all'
%                               is set)   
%       sensor_data.uz_min_all- minimum particle velocity in the
%                               z-direction recorded at all grid points
%                               within the domain (returned if 'u_min_all'
%                               is set) 
%       sensor_data.Ix        - time varying acoustic intensity in the
%                               x-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'I' is set) 
%       sensor_data.Iy        - time varying acoustic intensity in the
%                               y-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'I' is set) 
%       sensor_data.Iz        - time varying acoustic intensity in the
%                               z-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'I' is set) 
%       sensor_data.Ix_avg    - average acoustic intensity in the
%                               x-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'I_avg' is set)   
%       sensor_data.Iy_avg    - average acoustic intensity in the
%                               y-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'I_avg' is set) 
%       sensor_data.Iz_avg    - average acoustic intensity in the
%                               z-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'I_avg' is set) 
%
% ABOUT:
%       author      - Bradley Treeby and Ben Cox
%       date        - 7th April 2009
%       last update - 4th June 2013
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2013 Bradley Treeby and Ben Cox
%
% See also fftn, ifftn, imagesc, kspaceFirstOrder1D, kspaceFirstOrder2D,
% makeGrid, makeTime, makeTransducer, smooth, unmaskSensorData 

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>.

% suppress mlint warnings that arise from using subscripts
%#ok<*NASGU>
%#ok<*COLND>
%#ok<*NODEF>
%#ok<*INUSL>

% start the timer and store the start time
start_time = clock;
tic;

% update command line status
disp('Running k-Wave simulation...');
disp(['  start time: ' datestr(start_time)]);

% =========================================================================
% DEFINE LITERALS
% =========================================================================

% minimum number of input variables
NUM_REQ_INPUT_VARIABLES = 4;
 
% optional input defaults (these are assigned and then cleared after
% inputChecking)
CARTESIAN_INTERP_DEF = 'linear';
CREATE_LOG_DEF = false;
DATA_CAST_DEF = 'off';
DATA_RECAST_DEF = false;
DISPLAY_MASK_DEF = 'default';
LOG_SCALE_DEF = false;
LOG_SCALE_COMPRESSION_FACTOR_DEF = 0.02;
MOVIE_ARGS_DEF = {};
MOVIE_NAME_DEF = [getDateString '-kspaceFirstOrder3D'];
PLOT_FREQ_DEF = 10;
PLOT_LAYOUT_DEF = false;
PLOT_SCALE_DEF = [-1 1];
PLOT_SIM_DEF = true;
PLOT_PML_DEF = true;
PML_ALPHA_DEF = 2;
PML_INSIDE_DEF = true;
PML_SIZE_DEF = 10;
RECORD_MOVIE_DEF = false;
SAVE_TO_DISK_DEF = false;
SAVE_TO_DISK_FILENAME_DEF = 'kwave_input_data.mat';
SAVE_TO_DISK_EXIT_DEF = true;
SCALE_SOURCE_TERMS_DEF = true;
SMOOTH_P0_DEF = true;
SMOOTH_C0_DEF = false;
SMOOTH_RHO0_DEF = false;
SOURCE_P_MODE_DEF = 'additive';
SOURCE_U_MODE_DEF = 'additive';
STREAM_TO_DISK_DEF = false;
STREAM_TO_DISK_STEPS_DEF = 200;
USE_KSPACE_DEF = true;
USE_SG_DEF = true;

% set default movie compression (these are not cleared after inputChecking)
MOVIE_COMP_WIN = 'Cinepak';
MOVIE_COMP_MAC = 'None';
MOVIE_COMP_LNX = 'None';
MOVIE_COMP_64B = 'None';

% set additional literals (these are not cleared after inputChecking)
COLOR_MAP = getColorMap;
DT_WARNING_CFL = 0.45;  
ESTIMATE_SIM_TIME_STEPS = 50;
LOG_NAME = ['k-Wave-Log-' getDateString];
MFILE = mfilename;
PLOT_SCALE_WARNING = 20;
STREAM_TO_DISK_FILENAME = 'temp_sensor_data.bin';

% =========================================================================
% CHECK INPUT STRUCTURES AND OPTIONAL INPUTS
% =========================================================================

% run subscript to check inputs
kspaceFirstOrder_inputChecking;

% gpu memory counter for GPUmat toolbox
if strncmp(data_cast, 'kWaveGPU', 8);
    total_gpu_mem = GPUmem;
end

% =========================================================================
% UPDATE COMMAND LINE STATUS
% =========================================================================

disp(['  dt: ' scaleSI(dt) 's, t_end: ' scaleSI(t_array(end)) 's, time steps: ' num2str(length(t_array))]);
[x_sc, scale, prefix] = scaleSI(min([kgrid.x_size, kgrid.y_size, kgrid.z_size])); %#ok<ASGLU>
disp(['  input grid size: ' num2str(kgrid.Nx) ' by ' num2str(kgrid.Ny) ' by ' num2str(kgrid.Nz) ' grid points (' num2str(kgrid.x_size*scale) ' by ' num2str(kgrid.y_size*scale) ' by ' num2str(kgrid.z_size*scale) prefix 'm)']); 
if (kgrid.kx_max == kgrid.kz_max) && (kgrid.kx_max == kgrid.ky_max)
    disp(['  maximum supported frequency: ' scaleSI( kgrid.k_max * min(c(:)) / (2*pi) ) 'Hz']);
else
    disp(['  maximum supported frequency: ' scaleSI( kgrid.kx_max * min(c(:)) / (2*pi) ) 'Hz by ' scaleSI( kgrid.ky_max * min(c(:)) / (2*pi) ) 'Hz by ' scaleSI( kgrid.kz_max * min(c(:)) / (2*pi) ) 'Hz']);
end

% =========================================================================
% SMOOTH AND ENLARGE INPUT GRIDS
% =========================================================================

% smooth the initial pressure distribution p0 if required, and then restore
% the maximum magnitude (NOTE: if p0 has any values at the edge of the
% domain, the smoothing may cause part of p0 to wrap to the other side of
% the domain) 
if isfield(source, 'p0') && smooth_p0
    disp('  smoothing p0 distribution...');      
    source.p0 = smooth(kgrid, source.p0, true);
end

% expand the computational grid if the PML is set to be outside the input
% grid defined by the user (kgrid.t_array is re-appended to the expanded
% grid)
if ~PML_inside

    % expand the computational grid, retaining the values for
    % kgrid.t_array
    disp('  expanding computational grid...');
    t_array_temp = kgrid.t_array;
    kgrid = makeGrid(kgrid.Nx + 2*PML_x_size, kgrid.dx, kgrid.Ny + 2*PML_y_size, kgrid.dy, kgrid.Nz + 2*PML_z_size, kgrid.dz);
    kgrid.t_array = t_array_temp;
    clear t_array_temp;
               
    % assign Nt and dt to kgrid if given as a structure
    if isstruct(kgrid)
        kgrid.Nt = length(kgrid.t_array);
        kgrid.dt = kgrid.t_array(2) - kgrid.t_array(1);
    end      
    
    % expand the grid matrices allowing a different PML size in each
    % Cartesian direction
    expand_size = [PML_x_size, PML_y_size, PML_z_size]; %#ok<NASGU>
    kspaceFirstOrder_expandGridMatrices;
    clear expand_size;
    
    % update command line status
    disp(['  computational grid size: ' num2str(kgrid.Nx) ' by ' num2str(kgrid.Ny) ' by ' num2str(kgrid.Nz) ' grid points']);

end

% define index variables to remove the PML from the display if the optional
% input 'PlotPML' is set to false
if ~plot_PML
    x1 = (PML_x_size + 1);
    x2 = kgrid.Nx - PML_x_size;
    y1 = (PML_y_size + 1);
    y2 = kgrid.Ny - PML_y_size;    
    z1 = (PML_z_size + 1);
    z2 = kgrid.Nz - PML_z_size;
else
    x1 = 1;
    x2 = kgrid.Nx;
    y1 = 1;
    y2 = kgrid.Ny;
    z1 = 1;
    z2 = kgrid.Nz;
end    

% define index variables to allow original grid size to be maintained for
% the _final and _all output variabkes if 'PMLInside' is set to false
if ~PML_inside
    x1_inside = (PML_x_size + 1);
    x2_inside = kgrid.Nx - PML_x_size;
    y1_inside = (PML_y_size + 1);
    y2_inside = kgrid.Ny - PML_y_size;    
    z1_inside = (PML_z_size + 1);
    z2_inside = kgrid.Nz - PML_z_size;
else
    x1_inside = 1;
    x2_inside = kgrid.Nx;
    y1_inside = 1;
    y2_inside = kgrid.Ny;
    z1_inside = 1;
    z2_inside = kgrid.Nz;
end

% smooth the sound speed distribution if required
if smooth_c && numDim(c) == 3
    disp('  smoothing sound speed distribution...');      
    c = smooth(kgrid, c);
end
    
% smooth the ambient density distribution if required
if smooth_rho0 && numDim(rho0) == 3
    disp('  smoothing density distribution...');      
    rho0 = smooth(kgrid, rho0);
end

% =========================================================================
% PREPARE STAGGERED COMPUTATIONAL GRIDS AND OPERATORS
% =========================================================================

% interpolate the values of the density at the staggered grid locations
% where sgx = (x + dx/2, y, z), sgy = (x, y + dy/2, z), sgz = (x, y, z +
% dz/2)
if numDim(rho0) == 3 && use_sg
    
    % rho0 is heterogeneous and staggered grids are used
    rho0_sgx = interpn(kgrid.x, kgrid.y, kgrid.z, rho0, kgrid.x + kgrid.dx/2, kgrid.y, kgrid.z, '*linear');
    rho0_sgy = interpn(kgrid.x, kgrid.y, kgrid.z, rho0, kgrid.x, kgrid.y + kgrid.dy/2, kgrid.z, '*linear');
    rho0_sgz = interpn(kgrid.x, kgrid.y, kgrid.z, rho0, kgrid.x, kgrid.y, kgrid.z + kgrid.dz/2, '*linear');
    
    % set values outside of the interpolation range to original values
    rho0_sgx(isnan(rho0_sgx)) = rho0(isnan(rho0_sgx));
    rho0_sgy(isnan(rho0_sgy)) = rho0(isnan(rho0_sgy));    
    rho0_sgz(isnan(rho0_sgz)) = rho0(isnan(rho0_sgz));
    
else
    % rho0 is homogeneous or staggered grids are not used
    rho0_sgx = rho0;
    rho0_sgy = rho0;
    rho0_sgz = rho0;
end

% get the PML operators based on the reference sound speed and PML settings
pml_x = getPML(kgrid.Nx, kgrid.dx, kgrid.dt, c_ref, PML_x_size, PML_x_alpha, false, 1);
pml_x_sgx = getPML(kgrid.Nx, kgrid.dx, kgrid.dt, c_ref, PML_x_size, PML_x_alpha, true && use_sg, 1);
pml_y = getPML(kgrid.Ny, kgrid.dy, kgrid.dt, c_ref, PML_y_size, PML_y_alpha, false, 2);
pml_y_sgy = getPML(kgrid.Ny, kgrid.dy, kgrid.dt, c_ref, PML_y_size, PML_y_alpha, true && use_sg, 2);
pml_z = getPML(kgrid.Nz, kgrid.dz, kgrid.dt, c_ref, PML_z_size, PML_z_alpha, false, 3);
pml_z_sgz = getPML(kgrid.Nz, kgrid.dz, kgrid.dt, c_ref, PML_z_size, PML_z_alpha, true && use_sg, 3);

% define the k-space derivative operators, multiply by the staggered
% grid shift operators, and then re-order using ifftshift (the option
% use_sg exists for debugging) 
if use_sg
    ddx_k_shift_pos = ifftshift( 1i*kgrid.kx_vec .* exp(1i*kgrid.kx_vec*kgrid.dx/2) );
    ddx_k_shift_neg = ifftshift( 1i*kgrid.kx_vec .* exp(-1i*kgrid.kx_vec*kgrid.dx/2) );
    ddy_k_shift_pos = ifftshift( 1i*kgrid.ky_vec .* exp(1i*kgrid.ky_vec*kgrid.dy/2) );
    ddy_k_shift_neg = ifftshift( 1i*kgrid.ky_vec .* exp(-1i*kgrid.ky_vec*kgrid.dy/2) );
    ddz_k_shift_pos = ifftshift( 1i*kgrid.kz_vec .* exp(1i*kgrid.kz_vec*kgrid.dz/2) );
    ddz_k_shift_neg = ifftshift( 1i*kgrid.kz_vec .* exp(-1i*kgrid.kz_vec*kgrid.dz/2) );
else
    ddx_k_shift_pos = ifftshift( 1i*kgrid.kx_vec );
    ddx_k_shift_neg = ifftshift( 1i*kgrid.kx_vec );
    ddy_k_shift_pos = ifftshift( 1i*kgrid.ky_vec );
    ddy_k_shift_neg = ifftshift( 1i*kgrid.ky_vec );
    ddz_k_shift_pos = ifftshift( 1i*kgrid.kz_vec );
    ddz_k_shift_neg = ifftshift( 1i*kgrid.kz_vec );         
end

% create k-space operator (the option use_kspace exists for debugging)
if use_kspace
    kappa = ifftshift( sinc(c_ref*dt*kgrid.k/2) );
else
    kappa = 1;
end

% force the derivative operators and shift oeprators to be in the correct
% direction for use with BSXFUN
ddy_k_shift_pos = ddy_k_shift_pos.'; 
ddy_k_shift_neg = ddy_k_shift_neg.';
ddz_k_shift_pos = permute(ddz_k_shift_pos, [2 3 1]);
ddz_k_shift_neg = permute(ddz_k_shift_neg, [2 3 1]);

% cleanup unused variables
clear ax* ay* az* x0_min x0_max y0_min y0_max z0_min z0_max;

% =========================================================================
% GET SENSOR MASK INDEX
% ========================================================================

% run subscript to get the indices of the active sensor elements
if use_sensor
    kspaceFirstOrder_getSensorMaskIndex;
end

% =========================================================================
% SCALE THE SOURCE TERMS
% =========================================================================

% run subscript to scale the source terms to the correct units
if scale_source_terms
    kspaceFirstOrder_scaleSourceTerms;
end

% =========================================================================
% SAVE DATA TO DISK FOR RUNNING SIMULATION EXTERNAL TO MATLAB
% =========================================================================

% save to disk option for saving the input matrices to disk for running
% simulations using k-Wave++
if save_to_disk
    % run subscript to save files to disk
    kspaceFirstOrder_saveToDisk;
    
    % run subscript to resize the transducer object if the grid has been
    % expanded 
    kspaceFirstOrder_retractTransducerGridSize;
    
    % exit matlab computation if required
    if save_to_disk_exit
        return
    end
end

% =========================================================================
% PREPARE DATA MASKS AND STORAGE VARIABLES
% =========================================================================

% run subscript to create acoustic absorption variables
kspaceFirstOrder_createAbsorptionVariables;

% run subscript to create storage variables
kspaceFirstOrder_createStorageVariables;

% =========================================================================
% DATA CASTING
% =========================================================================

% preallocate the loop variables using the castZeros anonymous function
% (this creates a matrix of zeros in the data type specified by data_cast)
p = castZeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
rhox = castZeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
rhoy = castZeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
rhoz = castZeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
ux_sgx = castZeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
uy_sgy = castZeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
uz_sgz = castZeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
p_k = castZeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);

% run subscript to cast the remaining loop variables to the data type
% specified by data_cast 
if ~strcmp(data_cast, 'off')
    kspaceFirstOrder_dataCast;
end

% =========================================================================
% CREATE INDEX VARIABLES
% =========================================================================

% setup the time index variable
if ~time_rev
    index_start = 1;
    index_step = 1;
    index_end = length(t_array); 
else
    % reverse the order of the input data
    sensor.time_reversal_boundary_data = fliplr(sensor.time_reversal_boundary_data);
    index_start = 1;
    index_step = 1;
    
    % stop one time point before the end so the last points are not
    % propagated
    index_end = length(t_array) - 1;
end

% =========================================================================
% PREPARE VISUALISATIONS
% =========================================================================

% pre-compute suitable axes scaling factor
if plot_layout || plot_sim
    [x_sc, scale, prefix] = scaleSI(max([kgrid.x_vec; kgrid.y_vec; kgrid.z_vec]));  %#ok<ASGLU>
end

% run subscript to plot the simulation layout if 'PlotLayout' is set to true
if plot_layout
    kspaceFirstOrder_plotLayout;
end

% initialise the figures used for animation if 'PlotSim' is set to 'true'
if plot_sim
    img = figure;
    if ~time_rev
        pbar = waitbar(0, 'Computing Pressure Field', 'Visible', 'off');
    else
        pbar = waitbar(0, 'Computing Time Reversed Field', 'Visible', 'off');
    end
    
    % shift the waitbar so it doesn't overlap the figure window
    posn_pbar = get(pbar, 'OuterPosition');
    posn_img = get(img, 'OuterPosition');
    posn_pbar(2) = max(min(posn_pbar(2) - posn_pbar(4), posn_img(2) - posn_pbar(4) - 10), 0);
    set(pbar, 'OuterPosition', posn_pbar, 'Visible', 'on');
end 

% initialise movie parameters if 'RecordMovie' is set to 'true'
if record_movie
    kspaceFirstOrder_initialiseMovieParameters;
end

% =========================================================================
% LOOP THROUGH TIME STEPS
% =========================================================================

% update command line status
disp(['  precomputation completed in ' scaleTime(toc)]);
disp('  starting time loop...');

% restart timing variables
loop_start_time = clock;
tic;

% start time loop
for t_index = index_start:index_step:index_end

    % enforce time reversal bounday condition
    if time_rev

        % load pressure value and enforce as a Dirichlet boundary condition
        p(sensor_mask_index) = sensor.time_reversal_boundary_data(:, t_index);

        % update p_k
        p_k = fftn(p);

        % compute rhox and rhoz using an adiabatic equation of state
        rhox_mod = p./(3*c.^2);
        rhoy_mod = p./(3*c.^2);
        rhoz_mod = p./(3*c.^2);
        rhox(sensor_mask_index) = rhox_mod(sensor_mask_index);
        rhoy(sensor_mask_index) = rhoy_mod(sensor_mask_index);
        rhoz(sensor_mask_index) = rhoz_mod(sensor_mask_index);
           
    end    

    % calculate ux, uy and uz at the next time step using dp/dx, dp/dy and
    % dp/dz at the current time step
    ux_sgx = bsxfun(@times, pml_x_sgx, ...
        bsxfun(@times, pml_x_sgx, ux_sgx) ... 
        - dt./rho0_sgx .* real(ifftn( bsxfun(@times, ddx_k_shift_pos, kappa .* p_k) )) ...
        );
    uy_sgy = bsxfun(@times, pml_y_sgy, ...
        bsxfun(@times, pml_y_sgy, uy_sgy) ...
        - dt./rho0_sgy .* real(ifftn( bsxfun(@times, ddy_k_shift_pos, kappa .* p_k) )) ...
        );
    uz_sgz = bsxfun(@times, pml_z_sgz, ...
        bsxfun(@times, pml_z_sgz, uz_sgz) ...
        - dt./rho0_sgz .* real(ifftn( bsxfun(@times, ddz_k_shift_pos, kappa .* p_k) )) ...
        );                 
    
    % override lazy execution if using the Accelereyes GPU toolbox (this
    % has a significant effect on performance)
    if force_geval
        geval(ux_sgx, uy_sgy, uz_sgz);
    end    
    
    % add in the velocity source terms   
    if ux_source >= t_index
        if strcmp(source.u_mode, 'dirichlet')
            % enforce the source values as a dirichlet boundary condition
            ux_sgx(u_source_index) = source.ux(:, t_index);            
        else
            % add the source values to the existing field values
            ux_sgx(u_source_index) = ux_sgx(u_source_index) + source.ux(:, t_index);
        end
    end
    if uy_source >= t_index
        if strcmp(source.u_mode, 'dirichlet')
            % enforce the source values as a dirichlet boundary condition
            uy_sgy(u_source_index) = source.uy(:, t_index);            
        else
            % add the source values to the existing field values
            uy_sgy(u_source_index) = uy_sgy(u_source_index) + source.uy(:, t_index);
        end        
    end
    if uz_source >= t_index
        if strcmp(source.u_mode, 'dirichlet')
            % enforce the source values as a dirichlet boundary condition
            uz_sgz(u_source_index) = source.uz(:, t_index);            
        else
            % add the source values to the existing field values
            uz_sgz(u_source_index) = uz_sgz(u_source_index) + source.uz(:, t_index);
        end         
    end   
    
    % add in transducer source term; there will normally be less source
    % points than time points, so this is only done until the source points
    % run out (the number of source points is stored in transducer_source)
    if transducer_source >= t_index
        % as only flat transducers are currently supported, assume all the
        % energy is transfered to x-direction velocity, multiply source
        % terms by apodization weights
        ux_sgx(u_source_index) = ux_sgx(u_source_index) + transducer_transmit_apodization.*transducer_input_signal(delay_mask);
        
        % update the delay_mask - this maps the source positions belonging
        % to different transducer elements to time points within
        % transducer_input_signal (this is a single time series) based on the
        % beamforming and focussing delays.  
        delay_mask = delay_mask + 1;
    end
    
    % calculate dux/dx, duydy and duz/dz at the next time step
    duxdx = real(ifftn( bsxfun(@times, ddx_k_shift_neg, kappa .* fftn(ux_sgx)) ));
    duydy = real(ifftn( bsxfun(@times, ddy_k_shift_neg, kappa .* fftn(uy_sgy)) ));
    duzdz = real(ifftn( bsxfun(@times, ddz_k_shift_neg, kappa .* fftn(uz_sgz)) ));        

    % calculate rhox, rhoy and rhoz at the next time step
    if ~nonlinear
        % use linearised mass conservation equation
        rhox = bsxfun(@times, pml_x, bsxfun(@times, pml_x, rhox) - dt.*rho0 .* duxdx);
        rhoy = bsxfun(@times, pml_y, bsxfun(@times, pml_y, rhoy) - dt.*rho0 .* duydy);        
        rhoz = bsxfun(@times, pml_z, bsxfun(@times, pml_z, rhoz) - dt.*rho0 .* duzdz);
    else
        % use nonlinear mass conservation equation (implicit calculation)
        rhox = bsxfun(@times, pml_x, ( bsxfun(@times, pml_x, rhox) - dt.*rho0 .* duxdx ) ./ (1 + 2*dt.*duxdx));
        rhoy = bsxfun(@times, pml_y, ( bsxfun(@times, pml_y, rhoy) - dt.*rho0 .* duydy ) ./ (1 + 2*dt.*duydy));
        rhoz = bsxfun(@times, pml_z, ( bsxfun(@times, pml_z, rhoz) - dt.*rho0 .* duzdz ) ./ (1 + 2*dt.*duzdz));
    end 
            
    % override lazy execution if using the Accelereyes GPU toolbox (this
    % has a significant effect on performance)
    if force_geval
        geval(rhox, rhoy, rhoz);
    end
    
    % add in the pre-scaled pressure source term as a mass source   
    if p_source >= t_index
        if strcmp(source.p_mode, 'dirichlet')
            % enforce the source values as a dirichlet boundary condition
            rhox(p_source_index) = source.p(:, t_index);  
            rhoy(p_source_index) = source.p(:, t_index);  
            rhoz(p_source_index) = source.p(:, t_index);
        else
            % add the source values to the existing field values
            rhox(p_source_index) = rhox(p_source_index) + source.p(:, t_index);
            rhoy(p_source_index) = rhoy(p_source_index) + source.p(:, t_index);  
            rhoz(p_source_index) = rhoz(p_source_index) + source.p(:, t_index);
        end
    end
    
    % calculate p at the next time step
    if ~nonlinear
        switch equation_of_state
            case 'lossless';
                % calculate p using a linear adiabatic equation of state
                p = c.^2.*(rhox + rhoy + rhoz);
            case 'absorbing';
                % calculate p using a linear absorbing equation of state                
                p = c.^2.*( ...
                    (rhox + rhoy + rhoz) ...
                    + absorb_tau.*real(ifftn( absorb_nabla1.*fftn(rho0.*(duxdx + duydy + duzdz)) )) ...
                    - absorb_eta.*real(ifftn( absorb_nabla2.*fftn(rhox + rhoy + rhoz) )) ...
                    );  
        end
    else
        switch equation_of_state
            case 'lossless';
                % calculate p using a nonlinear adiabatic equation of state
                p = c.^2.*(rhox + rhoy + rhoz + medium.BonA.*(rhox + rhoy + rhoz).^2./(2*rho0));
            case 'absorbing';
                % calculate p using a nonlinear absorbing equation of state 
                p = c.^2.*(...
                    (rhox + rhoy + rhoz) ...
                    + absorb_tau.*real(ifftn( absorb_nabla1.*fftn(rho0.*(duxdx + duydy + duzdz)) ))...
                    - absorb_eta.*real(ifftn( absorb_nabla2.*fftn(rhox + rhoy + rhoz) ))...
                    + medium.BonA.*(rhox + rhoy + rhoz).^2./(2*rho0) ...
                    );
        end
    end
             
    % enforce initial conditions if source.p0 is defined instead of time
    % varying sources
    if t_index == 1 && isfield(source, 'p0')
    
        % add the initial pressure to rho as a mass source
        p = source.p0;
        rhox = source.p0./(3*c.^2);
        rhoy = source.p0./(3*c.^2);
        rhoz = source.p0./(3*c.^2);

        % compute u(t = t1 + dt/2) based on the assumption u(dt/2) = -u(-dt/2)
        % which forces u(t = t1) = 0
        ux_sgx = dt./rho0_sgx .* real(ifftn( bsxfun(@times, ddx_k_shift_pos, kappa .* fftn(p)) )) / 2;
        uy_sgy = dt./rho0_sgy .* real(ifftn( bsxfun(@times, ddy_k_shift_pos, kappa .* fftn(p)) )) / 2;
        uz_sgz = dt./rho0_sgz .* real(ifftn( bsxfun(@times, ddz_k_shift_pos, kappa .* fftn(p)) )) / 2;    
        
    end
    
    % precompute fft of p here so p can be modified for visualisation
    p_k = fftn(p);        

    % extract required sensor data from the pressure and particle velocity
    % fields if the number of time steps elapsed is greater than
    % sensor.record_start_index (defaults to 1) 
    if use_sensor && ~time_rev && (t_index >= sensor.record_start_index)
    
        % update index for data storage - if streaming to disk, a smaller
        % matrix is used which is continually overwritten, and then saved to
        % disk each time it is filled
        if stream_to_disk
            file_index = t_index - stream_to_disk*(stream_data_index - 1) - sensor.record_start_index + 1;
        else
            file_index = t_index - sensor.record_start_index + 1;
        end
        
        % store the acoustic pressure if using a transducer object
        if transducer_sensor
            
            % check if an elevation focus is set
            if transducer_receive_elevation_focus

                % update the sensor data buffer
                sensor_data_buffer = circshift(sensor_data_buffer, [0 1]);
                sensor_data_buffer(:, 1) = p(sensor_mask_index);

                % if buffer has been filled, store average pressure
                % across each transducer element accounting for the
                % elevation focus
                if file_index >= sensor_data_buffer_size
                    % get the current values
                    current_vals = sum(transducer_receive_mask.*sensor_data_buffer, 2);

                    % reshape and average
                    sensor_data.transducer(:, file_index - sensor_data_buffer_size + 1) = ...
                        sum(reshape(sum(reshape(current_vals, sensor.number_active_elements*sensor.element_width, sensor.element_length), 2), sensor.element_width, sensor.number_active_elements).', 2);                        
                end
                
            else
                
                % store average pressure across each transducer element
                sensor_data.transducer(:, file_index) = sum(reshape(sum(reshape(p(sensor_mask_index), sensor.number_active_elements*sensor.element_width, sensor.element_length), 2), sensor.element_width, sensor.number_active_elements).', 2);
                
            end
        
        % extract data from specified grid points
        elseif binary_sensor_mask
            
            % store the time history of the acoustic pressure
            if record.p
                sensor_data.p(:, file_index) = p(sensor_mask_index);
            end

            % store the maximum acoustic pressure
            if record.p_max
                if file_index == 1
                    sensor_data.p_max = p(sensor_mask_index);
                else
                    sensor_data.p_max = max(sensor_data.p_max, p(sensor_mask_index));
                end
            end
            
            % store the minimum acoustic pressure
            if record.p_min
                if file_index == 1
                    sensor_data.p_min = p(sensor_mask_index);
                else
                    sensor_data.p_min = min(sensor_data.p_min, p(sensor_mask_index));
                end                
            end

            % store the rms acoustic pressure
            if record.p_rms
                sensor_data.p_rms = sqrt((sensor_data.p_rms.^2*(file_index - 1) + p(sensor_mask_index).^2)./file_index);
            end
            
            % store the time history of the particle velocity
            if record.u
                sensor_data.ux(:, file_index) = ux_sgx(sensor_mask_index);
                sensor_data.uy(:, file_index) = uy_sgy(sensor_mask_index);            
                sensor_data.uz(:, file_index) = uz_sgz(sensor_mask_index);  
            end

            % store the maximum particle velocity
            if record.u_max
                if file_index == 1
                    sensor_data.ux_max = ux_sgx(sensor_mask_index);
                    sensor_data.uy_max = uy_sgy(sensor_mask_index);
                    sensor_data.uz_max = uz_sgz(sensor_mask_index);
                else
                    sensor_data.ux_max = max(sensor_data.ux_max, ux_sgx(sensor_mask_index));            
                    sensor_data.uy_max = max(sensor_data.uy_max, uy_sgy(sensor_mask_index));            
                    sensor_data.uz_max = max(sensor_data.uz_max, uz_sgz(sensor_mask_index));            
                end
            end        
            
            % store the minimum particle velocity
            if record.u_min
                if file_index == 1
                    sensor_data.ux_min = ux_sgx(sensor_mask_index);
                    sensor_data.uy_min = uy_sgy(sensor_mask_index);
                    sensor_data.uz_min = uz_sgz(sensor_mask_index);
                else
                    sensor_data.ux_min = min(sensor_data.ux_min, ux_sgx(sensor_mask_index));            
                    sensor_data.uy_min = min(sensor_data.uy_min, uy_sgy(sensor_mask_index));            
                    sensor_data.uz_min = min(sensor_data.uz_min, uz_sgz(sensor_mask_index));            
                end
            end             

            % store the rms particle velocity
            if record.u_rms
                sensor_data.ux_rms(:) = sqrt((sensor_data.ux_rms(:).^2*(file_index - 1) + ux_sgx(sensor_mask_index).^2)./file_index);
                sensor_data.uy_rms(:) = sqrt((sensor_data.uy_rms(:).^2*(file_index - 1) + uy_sgy(sensor_mask_index).^2)./file_index);
                sensor_data.uz_rms(:) = sqrt((sensor_data.uz_rms(:).^2*(file_index - 1) + uz_sgz(sensor_mask_index).^2)./file_index);
            end        
            
            % store the time history of the acoustic intensity
            if record.I
                if file_index > 1
                    sensor_data.Ix(:, file_index-1) = 0.25*(sensor_data.ux_prev_t + ux_sgx(sensor_mask_index) + ux_sgx(sensor_mask_index_sgx)).*sensor_data.p_prev_t;
                    sensor_data.Iy(:, file_index-1) = 0.25*(sensor_data.uy_prev_t + uy_sgy(sensor_mask_index) + uy_sgy(sensor_mask_index_sgy)).*sensor_data.p_prev_t;
                    sensor_data.Iz(:, file_index-1) = 0.25*(sensor_data.uz_prev_t + uz_sgz(sensor_mask_index) + uz_sgz(sensor_mask_index_sgz)).*sensor_data.p_prev_t;
                end
            end
            
            % store the average acoustic intensity
            if record.I_avg
                if file_index > 1
                    if record.I
                        sensor_data.Ix_avg = sensor_data.Ix_avg + sensor_data.Ix(:, file_index-1);
                        sensor_data.Iy_avg = sensor_data.Iy_avg + sensor_data.Iy(:, file_index-1);
                        sensor_data.Iz_avg = sensor_data.Iz_avg + sensor_data.Iz(:, file_index-1);
                    else
                        sensor_data.Ix_avg = sensor_data.Ix_avg + 0.25*(sensor_data.ux_prev_t + ux_sgx(sensor_mask_index) + ux_sgx(sensor_mask_index_sgx)).*sensor_data.p_prev_t;
                        sensor_data.Iy_avg = sensor_data.Iy_avg + 0.25*(sensor_data.uy_prev_t + uy_sgy(sensor_mask_index) + uy_sgy(sensor_mask_index_sgy)).*sensor_data.p_prev_t;
                        sensor_data.Iz_avg = sensor_data.Iz_avg + 0.25*(sensor_data.uz_prev_t + uz_sgz(sensor_mask_index) + uz_sgz(sensor_mask_index_sgz)).*sensor_data.p_prev_t;
                    end
                end
            end
            
            % update the stored value of the particle velocity at the
            % previous time step if storing the intensity but not the time
            % varying particle velocity
            if (record.I || record.I_avg)
                sensor_data.ux_prev_t = ux_sgx(sensor_mask_index) + ux_sgx(sensor_mask_index_sgx);
                sensor_data.uy_prev_t = uy_sgy(sensor_mask_index) + uy_sgy(sensor_mask_index_sgy);
                sensor_data.uz_prev_t = uz_sgz(sensor_mask_index) + uz_sgz(sensor_mask_index_sgz);
                sensor_data.p_prev_t = p(sensor_mask_index);
            end
            
        % extract data from specified Cartesian coordinates using
        % interpolation (tri and bc are the Delaunay triangulation
        % and Barycentric coordinates returned by gridDataFast3D)    
        else
            
            % store the time history of the acoustic pressure
            if record.p
                sensor_data.p(:, file_index) = sum(p(tri) .* bc, 2);
            end

            % store the maximum acoustic pressure
            if record.p_max
                if file_index == 1
                    sensor_data.p_max = sum(p(tri) .* bc, 2);
                else
                    sensor_data.p_max = max(sensor_data.p_max, sum(p(tri) .* bc, 2));
                end
            end        

            % store the minimum acoustic pressure
            if record.p_min
                if file_index == 1
                    sensor_data.p_min = sum(p(tri) .* bc, 2);
                else
                    sensor_data.p_min = min(sensor_data.p_min, sum(p(tri) .* bc, 2));
                end
            end            
            
            % store the rms acoustic pressure
            if record.p_rms
                sensor_data.p_rms(:) = sqrt((sensor_data.p_rms(:).^2*(file_index - 1) + (sum(p(tri) .* bc, 2)).^2)./file_index);
            end

            % store the time history of the particle velocity
            if record.u
                sensor_data.ux(:, file_index) = sum(ux_sgx(tri) .* bc, 2);
                sensor_data.uy(:, file_index) = sum(uy_sgy(tri) .* bc, 2);
                sensor_data.uz(:, file_index) = sum(uz_sgz(tri) .* bc, 2);
            end

            % store the maximum particle velocity
            if record.u_max
                if file_index == 1
                    sensor_data.ux_max = sum(ux_sgx(tri) .* bc, 2);
                    sensor_data.uy_max = sum(uy_sgy(tri) .* bc, 2);
                    sensor_data.uz_max = sum(uz_sgz(tri) .* bc, 2);
                else
                    sensor_data.ux_max = max(sensor_data.ux_max, sum(ux_sgx(tri) .* bc, 2));
                    sensor_data.uy_max = max(sensor_data.uy_max, sum(uy_sgy(tri) .* bc, 2));
                    sensor_data.uz_max = max(sensor_data.uz_max, sum(uz_sgz(tri) .* bc, 2));
                end
            end   
            
            % store the minimum particle velocity
            if record.u_min
                if file_index == 1
                    sensor_data.ux_min = sum(ux_sgx(tri) .* bc, 2);
                    sensor_data.uy_min = sum(uy_sgy(tri) .* bc, 2);
                    sensor_data.uz_min = sum(uz_sgz(tri) .* bc, 2);
                else
                    sensor_data.ux_min = min(sensor_data.ux_min, sum(ux_sgx(tri) .* bc, 2));
                    sensor_data.uy_min = min(sensor_data.uy_min, sum(uy_sgy(tri) .* bc, 2));
                    sensor_data.uz_min = min(sensor_data.uz_min, sum(uz_sgz(tri) .* bc, 2));
                end
            end             

            % store the rms particle velocity
            if record.u_rms
                sensor_data.ux_rms(:) = sqrt((sensor_data.ux_rms(:).^2*(file_index - 1) + (sum(ux_sgx(tri) .* bc, 2)).^2)./file_index);
                sensor_data.uy_rms(:) = sqrt((sensor_data.uy_rms(:).^2*(file_index - 1) + (sum(uy_sgy(tri) .* bc, 2)).^2)./file_index);
                sensor_data.uz_rms(:) = sqrt((sensor_data.uz_rms(:).^2*(file_index - 1) + (sum(uz_sgz(tri) .* bc, 2)).^2)./file_index);
            end        

            % store the time history of the acoustic intensity
            if record.I
                if file_index > 1
                    sensor_data.Ix(:, file_index-1) = 0.5*(sensor_data.ux_prev_t + sum(ux_sgx(tri_sgx) .* bc_sgx, 2) ).* sensor_data.p_prev_t;
                    sensor_data.Iy(:, file_index-1) = 0.5*(sensor_data.uy_prev_t + sum(uy_sgy(tri_sgy) .* bc_sgy, 2) ).* sensor_data.p_prev_t;
                    sensor_data.Iz(:, file_index-1) = 0.5*(sensor_data.uz_prev_t + sum(uz_sgz(tri_sgz) .* bc_sgz, 2) ).* sensor_data.p_prev_t;
                end
            end
            
            % store the average acoustic intensity
            if record.I_avg
                if file_index > 1
                    if record.I
                        sensor_data.Ix_avg = sensor_data.Ix_avg + sensor_data.Ix(:, file_index-1);
                        sensor_data.Iy_avg = sensor_data.Iy_avg + sensor_data.Iy(:, file_index-1);
                        sensor_data.Iz_avg = sensor_data.Iz_avg + sensor_data.Iz(:, file_index-1);
                    else
                        sensor_data.Ix_avg = sensor_data.Ix_avg + 0.5*(sensor_data.ux_prev_t + sum(ux_sgx(tri_sgx) .* bc_sgx, 2) ).* sensor_data.p_prev_t;
                        sensor_data.Iy_avg = sensor_data.Iy_avg + 0.5*(sensor_data.uy_prev_t + sum(uy_sgy(tri_sgy) .* bc_sgy, 2) ).* sensor_data.p_prev_t;
                        sensor_data.Iz_avg = sensor_data.Iz_avg + 0.5*(sensor_data.uz_prev_t + sum(uz_sgz(tri_sgz) .* bc_sgz, 2) ).* sensor_data.p_prev_t;
                    end
                end
            end
            
            % update the stored value of the particle velocity at the
            % previous time step if storing the intensity but not the time
            % varying particle velocity
            if (record.I || record.I_avg)
                sensor_data.ux_prev_t = sum(ux_sgx(tri_sgx) .* bc_sgx, 2);
                sensor_data.uy_prev_t = sum(uy_sgy(tri_sgy) .* bc_sgy, 2);
                sensor_data.uz_prev_t = sum(uz_sgz(tri_sgz) .* bc_sgz, 2);
                sensor_data.p_prev_t = sum(p(tri) .* bc, 2);
            end
            
        end
             
        % store the maximum acoustic pressure over all the grid
        % elements
        if record.p_max_all
            if file_index == 1
                sensor_data.p_max_all = p(x1_inside:x2_inside, y1_inside:y2_inside, z1_inside:z2_inside);
            else
                sensor_data.p_max_all = max(sensor_data.p_max_all, p(x1_inside:x2_inside, y1_inside:y2_inside, z1_inside:z2_inside));
            end
        end

        % store the minimum acoustic pressure over all the grid
        % elements
        if record.p_min_all
            if file_index == 1
                sensor_data.p_min_all = p(x1_inside:x2_inside, y1_inside:y2_inside, z1_inside:z2_inside);
            else
                sensor_data.p_min_all = min(sensor_data.p_min_all, p(x1_inside:x2_inside, y1_inside:y2_inside, z1_inside:z2_inside));
            end                
        end              
        
        % store the maximum particle velocity over all the grid
        % elements
        if record.u_max_all
            if file_index == 1
                sensor_data.ux_max_all = ux_sgx(x1_inside:x2_inside, y1_inside:y2_inside, z1_inside:z2_inside);
                sensor_data.uy_max_all = uy_sgy(x1_inside:x2_inside, y1_inside:y2_inside, z1_inside:z2_inside);
                sensor_data.uz_max_all = uz_sgz(x1_inside:x2_inside, y1_inside:y2_inside, z1_inside:z2_inside);
            else
                sensor_data.ux_max_all = max(sensor_data.ux_max_all, ux_sgx(x1_inside:x2_inside, y1_inside:y2_inside, z1_inside:z2_inside));            
                sensor_data.uy_max_all = max(sensor_data.uy_max_all, uy_sgy(x1_inside:x2_inside, y1_inside:y2_inside, z1_inside:z2_inside));            
                sensor_data.uz_max_all = max(sensor_data.uz_max_all, uz_sgz(x1_inside:x2_inside, y1_inside:y2_inside, z1_inside:z2_inside));            
            end
        end        

        % store the minimum particle velocity over all the grid
        % elements
        if record.u_min_all
            if file_index == 1
                sensor_data.ux_min_all = ux_sgx(x1_inside:x2_inside, y1_inside:y2_inside, z1_inside:z2_inside);
                sensor_data.uy_min_all = uy_sgy(x1_inside:x2_inside, y1_inside:y2_inside, z1_inside:z2_inside);
                sensor_data.uz_min_all = uz_sgz(x1_inside:x2_inside, y1_inside:y2_inside, z1_inside:z2_inside);
            else
                sensor_data.ux_min_all = min(sensor_data.ux_min_all, ux_sgx(x1_inside:x2_inside, y1_inside:y2_inside, z1_inside:z2_inside));            
                sensor_data.uy_min_all = min(sensor_data.uy_min_all, uy_sgy(x1_inside:x2_inside, y1_inside:y2_inside, z1_inside:z2_inside));            
                sensor_data.uz_min_all = min(sensor_data.uz_min_all, uz_sgz(x1_inside:x2_inside, y1_inside:y2_inside, z1_inside:z2_inside));            
            end
        end          
        
        % if the data is being streamed to disk and sensor_data has just been
        % filled, append the values of sensor_data to the values already saved
        % to disk (this option is currently only supported for recording
        % the time varying pressure)
        if stream_to_disk && (file_index == stream_to_disk)

            % open the file to append the new data
            if stream_data_index == 1
                % create or open the file and overwrite any existing data
                try
                    fid = fopen(STREAM_TO_DISK_FILENAME, 'w+');
                catch ME
                    disp('Error in writing file using ''StreamToDisk''');
                    rethrow(ME);
                end
            else
                % open the file and append new data
                fid = fopen(STREAM_TO_DISK_FILENAME, 'a+');
            end

            % write values at end of file using the precision specified by
            % the data_cast option. If using the GPU, the PCT requires a
            % call to gather, while GPUmat and Accelereyes only require a
            % cast back to single or double.
            if strcmp(data_cast, 'gpuArray')
                if strcmp(precision, 'single')
                    fwrite(fid, gather(sensor_data.p), 'single');
                else
                    fwrite(fid, gather(sensor_data.p), 'double');
                end
            elseif strcmp(precision, 'single')
                fwrite(fid, single(sensor_data.p), 'single');
            else
                fwrite(fid, double(sensor_data.p), 'double');       
            end

            % close the file
            fclose(fid);

            % increment the data file index if there is
            % still data left
            if index_end > t_index
                stream_data_index = stream_data_index + 1;
            end

        end
    end

    % estimate the time to run the simulation
    if t_index == ESTIMATE_SIM_TIME_STEPS
        disp(['  estimated simulation time ' scaleTime(etime(clock, loop_start_time)*index_end/t_index) '...']);

        % display current matlab memory usage
        if nargout == 2 && strncmp(computer, 'PCWIN', 5)
            [mem_usage.user, mem_usage.sys] = memory;
            disp(['  memory used: ' num2str(mem_usage.user.MemUsedMATLAB./1024^3) ' GB (of ' num2str(mem_usage.sys.PhysicalMemory.Total./1024^3) ' GB)']); 
        end
                
        % gpu memory counter for GPUmat toolbox
        if strncmp(data_cast, 'kWaveGPU', 8)
            current_gpu_mem = GPUmem;
            disp(['  GPU memory used: ' num2str((total_gpu_mem - current_gpu_mem)/2^30) ' GB (of ' num2str(total_gpu_mem/2^30) ' GB)']);
            mem_usage.gpu.total = total_gpu_mem;
            mem_usage.gpu.used = total_gpu_mem - current_gpu_mem;
        end
        
        % gpu memory counter for Accelereyes toolbox
        if strcmp(data_cast, 'gsingle') || strcmp(data_cast, 'gdouble')
            gpu_info = ginfo(true);
            disp(['  GPU memory used: ' num2str((gpu_info.gpu_total - gpu_info.gpu_free)/2^30) ' GB (of ' num2str(gpu_info.gpu_total/2^30) ' GB)']);
            mem_usage.gpu.total = gpu_info.gpu_total;
            mem_usage.gpu.used = gpu_info.gpu_total - gpu_info.gpu_free;            
        end
        
        % gpu memory counter for Parallel Computing toolbox
        if strcmp(data_cast, 'gpuArray')
            gpu_info = gpuDevice;
            disp(['  GPU memory used: ' num2str((gpu_info.TotalMemory - gpu_info.FreeMemory)/2^30) ' GB (of ' num2str(gpu_info.TotalMemory/2^30) ' GB)']);
            mem_usage.gpu.total = gpu_info.TotalMemory;
            mem_usage.gpu.used = gpu_info.TotalMemory - gpu_info.FreeMemory;            
        end        
    end
    
    % plot data if required
    if plot_sim && (rem(t_index, plot_freq) == 0 || t_index == 1 || t_index == index_end) 
        
        % update progress bar
        waitbar(t_index/length(t_array), pbar);
        drawnow;

        % ensure p is cast as a CPU variable and remove the PML from the
        % plot if required
        if strcmp(data_cast, 'gpuArray')
            p_plot = double(gather(p(x1:x2, y1:y2, z1:z2)));
        else
            p_plot = double(p(x1:x2, y1:y2, z1:z2));       
        end

        % update plot scale if set to automatic or log
        if plot_scale_auto || plot_scale_log
            kspaceFirstOrder_adjustPlotScale;
        end         
        
        % add display mask onto plot
        if strcmp(display_mask, 'default')
            p_plot(sensor.mask(x1:x2, y1:y2, z1:z2) ~= 0) = plot_scale(2);
        elseif ~strcmp(display_mask, 'off')
            p_plot(display_mask(x1:x2, y1:y2, z1:z2) ~= 0) = plot_scale(2);
        end     

        % update plot
        planeplot(scale*kgrid.x_vec(x1:x2), scale*kgrid.y_vec(y1:y2), scale*kgrid.z_vec(z1:z2), p_plot, '', plot_scale, prefix, COLOR_MAP);

        % save movie frames if required
        if record_movie
            
            % set background color to white
            set(gcf, 'Color', [1 1 1]);

            % save the movie frame
            movie_frames(frame_index) = getframe(gcf); %#ok<AGROW>

            % update frame index
            frame_index  = frame_index  + 1;
            
        end
        
        % update variable used for timing variable to exclude the first
        % time step if plotting is enabled
        if t_index == 1
            loop_start_time = clock;
        end        
    end
end

% assign the final time reversal values
if time_rev
    p(sensor_mask_index) = sensor.time_reversal_boundary_data(:, index_end + 1);
end

% update command line status
disp(['  simulation completed in ' scaleTime(toc)]);

% =========================================================================
% CLEAN UP
% =========================================================================

% clean up used figures
if plot_sim
    close(img);
    close(pbar);
    drawnow;
end

% save the movie frames to disk
if record_movie
    kspaceFirstOrder_saveMovieFile;   
end

% save the final particle velocity if required
if record.u_final
    sensor_data.ux_final = ux_sgx(x1_inside:x2_inside, y1_inside:y2_inside, z1_inside:z2_inside);
    sensor_data.uy_final = uy_sgy(x1_inside:x2_inside, y1_inside:y2_inside, z1_inside:z2_inside);
    sensor_data.uz_final = uz_sgz(x1_inside:x2_inside, y1_inside:y2_inside, z1_inside:z2_inside);
end

% save the final intensity values
if (record.I || record.I_avg)
    
    % update the values of the particle velocity at the next staggered time
    % point
    ux_sgx = bsxfun(@times, pml_x_sgx, ...
        bsxfun(@times, pml_x_sgx, ux_sgx) ... 
        - dt./rho0_sgx .* real(ifftn( bsxfun(@times, ddx_k_shift_pos, kappa .* p_k) )) ...
        );
    uy_sgy = bsxfun(@times, pml_y_sgy, ...
        bsxfun(@times, pml_y_sgy, uy_sgy) ...
        - dt./rho0_sgy .* real(ifftn( bsxfun(@times, ddy_k_shift_pos, kappa .* p_k) )) ...
        );
    uz_sgz = bsxfun(@times, pml_z_sgz, ...
        bsxfun(@times, pml_z_sgz, uz_sgz) ...
        - dt./rho0_sgz .* real(ifftn( bsxfun(@times, ddz_k_shift_pos, kappa .* p_k) )) ...
        );  
    
    % store the time history of the acoustic intensity
    if record.I
        if binary_sensor_mask
            sensor_data.Ix(:, file_index) = 0.25*(sensor_data.ux_prev_t + ux_sgx(sensor_mask_index) + ux_sgx(sensor_mask_index_sgx)).*sensor_data.p_prev_t;
            sensor_data.Iy(:, file_index) = 0.25*(sensor_data.uy_prev_t + uy_sgy(sensor_mask_index) + uy_sgy(sensor_mask_index_sgy)).*sensor_data.p_prev_t;
            sensor_data.Iz(:, file_index) = 0.25*(sensor_data.uz_prev_t + uz_sgz(sensor_mask_index) + uz_sgz(sensor_mask_index_sgz)).*sensor_data.p_prev_t;
        else
            sensor_data.Ix(:, file_index) = 0.5*(sensor_data.ux_prev_t + sum(ux_sgx(tri_sgx) .* bc_sgx, 2) ).*sensor_data.p_prev_t;
            sensor_data.Iy(:, file_index) = 0.5*(sensor_data.uy_prev_t + sum(uy_sgy(tri_sgy) .* bc_sgy, 2) ).*sensor_data.p_prev_t;
            sensor_data.Iz(:, file_index) = 0.5*(sensor_data.uz_prev_t + sum(uz_sgz(tri_sgz) .* bc_sgz, 2) ).*sensor_data.p_prev_t;
        end
    end

    % store the average acoustic intensity
    if record.I_avg
        if record.I
            sensor_data.Ix_avg = sensor_data.Ix_avg + sensor_data.Ix(:, file_index);
            sensor_data.Iy_avg = sensor_data.Iy_avg + sensor_data.Iy(:, file_index);
            sensor_data.Iz_avg = sensor_data.Iz_avg + sensor_data.Iz(:, file_index);
        else
            if binary_sensor_mask
                sensor_data.Ix_avg = sensor_data.Ix_avg + 0.25*(sensor_data.ux_prev_t + ux_sgx(sensor_mask_index) + ux_sgx(sensor_mask_index_sgx)).*sensor_data.p_prev_t;
                sensor_data.Iy_avg = sensor_data.Iy_avg + 0.25*(sensor_data.uy_prev_t + uy_sgy(sensor_mask_index) + uy_sgy(sensor_mask_index_sgy)).*sensor_data.p_prev_t;
                sensor_data.Iz_avg = sensor_data.Iz_avg + 0.25*(sensor_data.uz_prev_t + uz_sgz(sensor_mask_index) + uz_sgz(sensor_mask_index_sgz)).*sensor_data.p_prev_t;
            else
                sensor_data.Ix_avg = sensor_data.Ix_avg + 0.5*(sensor_data.ux_prev_t + sum(ux_sgx(tri_sgx) .* bc_sgx, 2) ).*sensor_data.p_prev_t;
                sensor_data.Iy_avg = sensor_data.Iy_avg + 0.5*(sensor_data.uy_prev_t + sum(uy_sgy(tri_sgy) .* bc_sgy, 2) ).*sensor_data.p_prev_t;
                sensor_data.Iz_avg = sensor_data.Iz_avg + 0.5*(sensor_data.uz_prev_t + sum(uz_sgz(tri_sgz) .* bc_sgz, 2) ).*sensor_data.p_prev_t;
            end
        end
    end
    
    % remove temporary variables
    sensor_data = rmfield(sensor_data, {'ux_prev_t', 'uy_prev_t', 'uz_prev_t', 'p_prev_t'});
    
end

% process the intensity data if required
if record.I_avg
    
    % divide the average values by the number of values recorded (the
    % variable num_recorded_time_points is set in createStorageVariables)
    sensor_data.Ix_avg = sensor_data.Ix_avg / num_recorded_time_points;
    sensor_data.Iy_avg = sensor_data.Iy_avg / num_recorded_time_points;
    sensor_data.Iz_avg = sensor_data.Iz_avg / num_recorded_time_points;
    
end

% save final sensor_data variable to disk if required
if stream_to_disk

    % clear some time loop variables before reloading sensor data to
    % free up memory (in case the sensor data is very large)
    clear duxdx duydy duzdz ux_sgx uy_sgy uz_sgz rhox rhoy rhoz;
    
    % double check the data has not just been saved
    if file_index ~= stream_to_disk
        
        % open the file to append the new data
        if stream_data_index == 1
            % open the file and overwrite any existing data
            fid = fopen(STREAM_TO_DISK_FILENAME, 'w+');
        else
            % open the file and append new data
            fid = fopen(STREAM_TO_DISK_FILENAME, 'a+');
        end
        
        % extract required data
        sensor_data.p = sensor_data.p(:, 1:file_index);
        
        % write values at end of file
        if strcmp(data_cast, 'gpuArray')
            if strcmp(precision, 'single')
                fwrite(fid, gather(sensor_data.p), 'single');
            else
                fwrite(fid, gather(sensor_data.p), 'double');
            end
        elseif strcmp(precision, 'single')
            fwrite(fid, single(sensor_data.p), 'single');
        else
            fwrite(fid, double(sensor_data.p), 'double');       
        end
            
        % close the file
        fclose(fid);
        
    end

    % reload complete streamed data and assign to sensor_data
    fid = fopen(STREAM_TO_DISK_FILENAME, 'r');
    if strcmp(precision, 'single')
        sensor_data.p = fread(fid, [sum(sensor.mask(:)), length(t_array)], '*single');            
    else
        sensor_data.p = fread(fid, [sum(sensor.mask(:)), length(t_array)], '*double');
    end
    fclose(fid);    
    
    % permanently delete the temporary storage
    disp('  removing temporary data...');
    old_state = recycle;
    recycle('off');
    delete(STREAM_TO_DISK_FILENAME);
    recycle(old_state);
   
end

% save the final pressure field if in time reversal mode
if time_rev
    record.p_final = true;
end

% save the final acoustic pressure if required
if record.p_final
    sensor_data.p_final = p(x1_inside:x2_inside, y1_inside:y2_inside, z1_inside:z2_inside);
end

% process the sensor data if recorded using a transducer
if transducer_sensor
    % the pressure values across the grid points in each element have
    % already been automatically summed, so now scale the summed pressure
    % data by the number of grid points in each element to convert to an
    % average 
    sensor_data.transducer = sensor_data.transducer./(sensor.element_length * sensor.element_width);
end

% run subscript to cast variables back to double precision if required
if data_recast
    kspaceFirstOrder_dataRecast;
end

% reorder the sensor points if a binary sensor mask was used for Cartesian
% sensor mask nearest neighbour interpolation (this is performed after
% recasting as the GPU toolboxes do not all support this subscript)
if use_sensor && reorder_data
    kspaceFirstOrder_reorderCartData;
end

% filter the recorded time domain pressure signals if transducer filter
% parameters are given 
if use_sensor && ~time_rev && isfield(sensor, 'frequency_response')
    sensor_data.p = gaussianFilter(sensor_data.p, 1/kgrid.dt, sensor.frequency_response(1), sensor.frequency_response(2));
end

% run subscript to resize the transducer object if the grid has been
% expanded 
kspaceFirstOrder_retractTransducerGridSize;

if ~use_sensor
    % if sensor is not used, return empty sensor data
    sensor_data = [];
elseif transducer_sensor
    % if using a transducer as sensor, reassign sensor_data.transducer to
    % sensor_data
    sensor_data = sensor_data.transducer;
elseif time_rev
    % if computing time reversal, reassign sensor_data.p_final to
    % sensor_data
    sensor_data = sensor_data.p_final;    
elseif ~isfield(sensor, 'record')
    % if sensor.record is not given by the user, reassign sensor_data.p to
    % sensor_data
    sensor_data = sensor_data.p;
end

% update command line status
disp(['  total computation time ' scaleTime(etime(clock, start_time))]);

% switch off log
if create_log
    diary off;
end

function planeplot(x_vec, y_vec, z_vec, data, data_title, plot_scale, prefix, color_map)
% Subfunction to produce a plot of a three-dimensional matrix through the
% three central planes
   
subplot(2, 2, 1), imagesc(y_vec, x_vec, squeeze(data(:, :, round(end/2))), plot_scale);
title([data_title 'x-y plane']);
axis image;

subplot(2, 2, 2), imagesc(z_vec, x_vec, squeeze(data(:, round(end/2), :)), plot_scale);
title('x-z plane');
axis image;
xlabel(['(All axes in ' prefix 'm)']);

subplot(2, 2, 3), imagesc(z_vec, y_vec, squeeze(data(round(end/2), :, :)), plot_scale);
title('y-z plane');
axis image;
colormap(color_map); 
drawnow;