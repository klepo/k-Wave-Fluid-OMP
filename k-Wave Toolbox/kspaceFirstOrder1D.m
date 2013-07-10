function sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, varargin)
%KSPACEFIRSTORDER1D     1D time-domain simulation of wave propagation.
%
% DESCRIPTION:
%       kspaceFirstOrder1D simulates the time-domain propagation of
%       compressional waves through a one-dimensional homogeneous or
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
%       source.ux.  
%
%       The field values are returned as arrays of time series at the
%       sensor locations defined by sensor.mask. This can be given either
%       as a binary matrix (i.e., a matrix of 1's and 0's with the same
%       dimensions as the computational grid) representing the grid points
%       within the computational grid that will collect the data, or as a
%       series of arbitrary Cartesian coordinates within the grid at which
%       the pressure values are calculated at each time step via
%       interpolation. The Cartesian points must be given as a 1 by N
%       matrix. If no output is required, the sensor input can be replaced
%       with an empty array []. If sensor.mask is given as a set of
%       Cartesian coordinates, the computed sensor_data is returned in the
%       same order. If sensor.mask is given as a binary matrix, sensor_data
%       is returned using MATLAB's standard column-wise linear matrix index
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
%       sensor_data.p, sensor_data.p_final, sensor_data.p_max, and
%       sensor_data.ux. Most of the output parameters are recorded at the
%       given sensor positions and are indexed as
%       sensor_data.field(sensor_point_index, time_index). The exceptions
%       are the averaged quantities ('p_max', 'p_rms', 'u_max', 'p_rms',
%       'I_avg'), the 'all' quantities ('p_max_all', 'p_min_all',
%       'u_max_all', 'u_min_all'), and the final quantities ('p_final',
%       'u_final'). The averaged quantities are indexed as
%       sensor_data.p_max(sensor_point_index), while the final and 'all'
%       quantities are returned over the entire grid and are indexed as
%       sensor_data.p_final(nx).
%
%       kspaceFirstOrder1D may also be used for time reversal image
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
%       directly to sensor.time_reversal_boundary_data  with the input
%       fields source.p0 and source.p removed or set to zero.
%
% USAGE:
%       sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor)
%       sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, ...) 
%
% INPUTS:
% The minimum fields that must be assigned to run an initial value problem
% (for example, a photoacoustic forward simulation) are marked with a *. 
%
%       kgrid*              - k-space grid structure returned by makeGrid
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
%                     that specified by sensor.mask.
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
%       'PlotLayout'- Boolean controlling whether a four panel plot of the
%                     initial simulation layout is produced (initial
%                     pressure, sensor mask, sound speed, density)
%                     (default = false).
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
%                     To remove the PML, set the appropriate PMLAlpha to
%                     zero rather than forcing the PML to be of zero size
%                     (default = 20). 
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
%       sensor_data.ux_max    - maximum particle velocity in the
%                               x-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'u_max' is set)   
%       sensor_data.ux_min    - minimum particle velocity in the
%                               x-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'u_min' is set)   
%       sensor_data.ux_rms    - rms of the time varying particle velocity
%                               in the x-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'u_rms' is set)     
%       sensor_data.ux_final  - final particle velocity field in the
%                               x-direction at all grid points within the
%                               domain (returned if 'u_final' is set) 
%       sensor_data.ux_max_all- maximum particle velocity in the
%                               x-direction recorded at all grid points
%                               within the domain (returned if 'u_max_all'  
%       sensor_data.ux_min_all- minimum particle velocity in the
%                               x-direction recorded at all grid points
%                               within the domain (returned if 'u_min_all'
%                               is set)   
%       sensor_data.Ix        - time varying acoustic intensity in the
%                               x-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'I' is set) 
%       sensor_data.Ix_avg    - average acoustic intensity in the
%                               x-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'I_avg' is set)   
%
% ABOUT:
%       author      - Bradley Treeby and Ben Cox
%       date        - 22nd April 2009
%       last update - 4th June 2013
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2013 Bradley Treeby and Ben Cox
%
% See also fft, ifft, getframe, kspaceFirstOrder2D, kspaceFirstOrder3D,
% makeGrid, makeTime, movie2avi, smooth, unmaskSensorData 

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

% suppress mlint warnings that arise from using sub-scripts
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

% optional input defaults
CARTESIAN_INTERP_DEF = 'linear';
CREATE_LOG_DEF = false;
DATA_CAST_DEF = 'off';
DATA_RECAST_DEF = false;
DISPLAY_MASK_DEF = 'default';
LOG_SCALE_DEF = false;
LOG_SCALE_COMPRESSION_FACTOR_DEF = 0.02;
MOVIE_ARGS_DEF = {};
MOVIE_NAME_DEF = [getDateString '-kspaceFirstOrder1D'];
PLOT_FREQ_DEF = 10;
PLOT_LAYOUT_DEF = false;
PLOT_SCALE_DEF = [-1.1 1.1];
PLOT_SIM_DEF = true;
PLOT_PML_DEF = true;
PML_ALPHA_DEF = 2;
PML_INSIDE_DEF = true;
PML_SIZE_DEF = 20;
RECORD_MOVIE_DEF = false;
SCALE_SOURCE_TERMS_DEF = true;
SMOOTH_P0_DEF = true;
SMOOTH_C0_DEF = false;
SMOOTH_RHO0_DEF = false;
SOURCE_P_MODE_DEF = 'additive';
SOURCE_U_MODE_DEF = 'additive';
USE_KSPACE_DEF = true;
USE_SG_DEF = true;

% set default movie compression
MOVIE_COMP_WIN = 'Cinepak';
MOVIE_COMP_MAC = 'None';
MOVIE_COMP_LNX = 'None';
MOVIE_COMP_64B = 'None';

% set additional literals
MFILE = mfilename;
DT_WARNING_CFL = 0.5; 
ESTIMATE_SIM_TIME_STEPS = 50;
LOG_NAME = ['k-Wave-Log-' getDateString];
PLOT_SCALE_WARNING = 5;

% =========================================================================
% CHECK INPUT STRUCTURES AND OPTIONAL INPUTS
% =========================================================================

% run subscript to check inputs
kspaceFirstOrder_inputChecking;

% =========================================================================
% UPDATE COMMAND LINE STATUS
% =========================================================================

disp(['  dt: ' scaleSI(dt) 's, t_end: ' scaleSI(t_array(end)) 's, time steps: ' num2str(length(t_array))]);
disp(['  input grid size: ' num2str(kgrid.Nx) ' grid points (' scaleSI(kgrid.x_size) 'm)']);
disp(['  maximum supported frequency: ' scaleSI( kgrid.k_max * min(c(:)) / (2*pi) ) 'Hz']); 

% =========================================================================
% SMOOTH AND ENLARGE INPUT GRIDS
% =========================================================================

% smooth the initial pressure distribution p0 if required, and then restore
% the maximum magnitude (NOTE: if p0 has any values at the edge of the
% domain, the smoothing may cause part of p0 to wrap to the other side of
% the domain) 
if smooth_p0 && ~time_rev
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
    kgrid = makeGrid(kgrid.Nx + 2*PML_x_size, kgrid.dx);
    kgrid.t_array = t_array_temp;
    clear t_array_temp;    
           
    % assign Nt and dt to kgrid if given as a structure
    if isstruct(kgrid)
        kgrid.Nt = length(kgrid.t_array);
        kgrid.dt = kgrid.t_array(2) - kgrid.t_array(1);
    end       
    
    % expand the grid matrices
    expand_size = PML_x_size; %#ok<NASGU>
    kspaceFirstOrder_expandGridMatrices;
    clear expand_size;
        
    % update command line status
    disp(['  computation grid size: ' num2str(kgrid.Nx) ' grid points']);
    
end

% define index values to remove the PML from the display if the optional
% input 'PlotPML' is set to false
if ~plot_PML
    % create indexes to allow the source input to be placed into the larger
    % simulation grid
    x1 = (PML_x_size + 1);
    x2 = kgrid.Nx - PML_x_size;
else
    % create indexes to place the source input exactly into the simulation
    % grid
    x1 = 1;
    x2 = kgrid.Nx;
end    

% smooth the sound speed distribution if required
if smooth_c && numel(c) > 1
    disp('  smoothing sound speed distribution...');     
    c = smooth(kgrid, c);
end
    
% smooth the ambient density distribution if required
if smooth_rho0 && numel(rho0) > 1
    disp('  smoothing density distribution...');     
    rho0 = smooth(kgrid, rho0);
end

% =========================================================================
% PREPARE STAGGERED COMPUTATIONAL GRIDS AND OPERATORS
% =========================================================================

% interpolate the values of the density at the staggered grid locations
% where sgx = (x + dx/2)
if numel(rho0) > 1 && use_sg
    
    % rho0 is heterogeneous and staggered grids are used
    rho0_sgx = interp1(kgrid.x, rho0, kgrid.x + kgrid.dx/2, '*linear');
    
    % set values outside of the interpolation range to original values 
    rho0_sgx(isnan(rho0_sgx)) = rho0(isnan(rho0_sgx)); 
    
else
    % rho0 is homogeneous or staggered grids are not used
    rho0_sgx = rho0;
end

% get the PML operators based on the reference sound speed and PML settings
pml_x = getPML(kgrid.Nx, kgrid.dx, kgrid.dt, c_ref, PML_x_size, PML_x_alpha, false, 1);
pml_x_sgx = getPML(kgrid.Nx, kgrid.dx, kgrid.dt, c_ref, PML_x_size, PML_x_alpha, true && use_sg, 1);

% define the k-space derivative operator
ddx_k = ifftshift(1i*kgrid.kx_vec);

% define the staggered grid shift operators (the option use_sg exists for
% debugging)
if use_sg
    shift_pos = ifftshift( exp(1i*kgrid.kx_vec*kgrid.dx/2) );
    shift_neg = ifftshift( exp(-1i*kgrid.kx_vec*kgrid.dx/2) );
else
    shift_pos = 1;
    shift_neg = 1;
end

% create k-space operator (the option use_kspace exists for debugging)
if use_kspace
    kappa = ifftshift( sinc(c_ref*dt*kgrid.k/2) );
else
    kappa = 1;
end

% clean up unused variables
clear ax* x0_min x0_max PML_x_alpha;

% =========================================================================
% GET SENSOR MASK INDEX
% ========================================================================

% run subscript to get the indices of the active sensor elements
if use_sensor
    kspaceFirstOrder_getSensorMaskIndex
end

% =========================================================================
% SCALE THE SOURCE TERMS
% =========================================================================
 
% run subscript to scale the source terms to the correct units
if scale_source_terms
    kspaceFirstOrder_scaleSourceTerms;
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
p = castZeros([kgrid.Nx, 1]);
rhox = castZeros([kgrid.Nx, 1]);
ux_sgx = castZeros([kgrid.Nx, 1]);
p_k = castZeros([kgrid.Nx, 1]);

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
    [x_sc, scale, prefix] = scaleSI(max(kgrid.x));  %#ok<ASGLU>
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
        p_k = fft(p);

        % compute rhox using an adiabatic equation of state
        rhox_mod = p./(c.^2);
        rhox(sensor_mask_index) = rhox_mod(sensor_mask_index);
        
    end   
    
    % calculate ux at the next time step using dp/dx at the current time
    % step
    if ~nonuniform_grid
        % calculate gradient using a regular grid
        ux_sgx = pml_x_sgx .* (  pml_x_sgx.*ux_sgx - dt./rho0_sgx .* real(ifft(ddx_k .* shift_pos .* kappa .* p_k))  );
    else
        % calculate gradients using a non-uniform grid via the mapped
        % pseudospectral method         
        ux_sgx = pml_x_sgx .* (  pml_x_sgx.*ux_sgx - dt./rho0_sgx .* kgrid.dxudxn_sgx .* real(ifft(ddx_k .* shift_pos .* kappa .* p_k))  );
    end    
    
    % add in the velocity source term
    if ux_source >= t_index
        if strcmp(source.u_mode, 'dirichlet')
            % enforce the source values as a dirichlet boundary condition
            ux_sgx(u_source_index) = source.ux(:, t_index);
        else
            % add the source values to the existing field values        
            ux_sgx(u_source_index) = ux_sgx(u_source_index) + source.ux(:, t_index);
        end
    end    

    % calculate du/dx at the next time step
    if ~nonuniform_grid
        % calculate gradients using a regular grid    
        duxdx = real(ifft(ddx_k .* shift_neg .* kappa .* fft(ux_sgx)));   
    else
        % calculate gradients using a non-uniform grid via the mapped
        % pseudospectral method 
        duxdx = kgrid.dxudxn .* real(ifft(ddx_k .* shift_neg .* kappa .* fft(ux_sgx)));
    end    
    
    % calculate rhox at the next time step
    if ~nonlinear
        % use linearised mass conservation equation
        rhox = pml_x .* ( pml_x .* rhox - dt.*rho0 .* duxdx );
    else
        % use nonlinear mass conservation equation (implicit calculation)
        rhox = pml_x .* ( ( pml_x .* rhox - dt.*rho0 .* duxdx) ./ (1 + 2*dt.*duxdx) );
    end      
    
    % add in the pre-scaled pressure source term as a mass source  
    if p_source >= t_index
        if strcmp(source.p_mode, 'dirichlet')
            % enforce the source values as a dirichlet boundary condition
            rhox(p_source_index) = source.p(:, t_index);
        else
            % add the source values to the existing field values
            rhox(p_source_index) = rhox(p_source_index) + source.p(:, t_index);
        end
    end    
    
    % equation of state
    if ~nonlinear
        switch equation_of_state
            case 'lossless'
                % compute p using an adiabatic equation of state
                p = c.^2.*rhox;
            case 'absorbing'
                % compute p using an absorbing equation of state 
                p = c.^2.*(rhox ...
                     + absorb_tau.*real(ifft( absorb_nabla1.*fft(rho0.*duxdx) )) ...
                     - absorb_eta.*real(ifft( absorb_nabla2.*fft(rhox) )) ...
                     );
        end
    else
        switch equation_of_state
            case 'lossless'
                % compute p using a nonlinear adiabatic equation of state
                p = c.^2.*( rhox + medium.BonA.*rhox.^2./(2*rho0) );
            case 'absorbing'
                % compute p using a nonlinear absorbing equation of state 
                p = c.^2.*( rhox ...
                    + absorb_tau.*real(ifft( absorb_nabla1.*fft(rho0.*duxdx) )) ...
                    - absorb_eta.*real(ifft( absorb_nabla2.*fft(rhox) )) ...                   
                    + medium.BonA.*rhox.^2./(2*rho0) ...
                    );
        end         
    end     
    
    % enforce initial conditions if source.p0 is defined instead of time
    % varying sources
    if t_index == 1 && isfield(source, 'p0')
    
        % add the initial pressure to rho as a mass source
        p = source.p0;
        rhox = source.p0./c.^2;
        
        % compute u(t = t1 - dt/2) based on u(dt/2) = -u(-dt/2) which
        % forces u(t = t1) = 0 
        ux_sgx = dt./rho0_sgx .* real(ifft(ddx_k .* shift_pos .* kappa .* fft(p))) / 2;

    end    
    
    % precompute fft of p here so p can be modified for visualisation
    p_k = fft(p);    
    
    
    % extract required sensor data from the pressure and particle velocity
    % fields if the number of time steps elapsed is greater than
    % sensor.record_start_index (defaults to 1) 
    if use_sensor && ~time_rev && (t_index >= sensor.record_start_index)
    
        % update index for data storage
        file_index = t_index - sensor.record_start_index + 1;
        
        % extract data from specified grid points
        if binary_sensor_mask
            
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
            end

            % store the maximum particle velocity
            if record.u_max
                if file_index == 1
                    sensor_data.ux_max = ux_sgx(sensor_mask_index);
                else
                    sensor_data.ux_max = max(sensor_data.ux_max, ux_sgx(sensor_mask_index));            
                end
            end        

            % store the minimum particle velocity
            if record.u_min
                if file_index == 1
                    sensor_data.ux_min = ux_sgx(sensor_mask_index);
                else
                    sensor_data.ux_min = min(sensor_data.ux_min, ux_sgx(sensor_mask_index));            
                end
            end              
            
            % store the rms particle velocity
            if record.u_rms
                sensor_data.ux_rms = sqrt((sensor_data.ux_rms.^2*(file_index - 1) + ux_sgx(sensor_mask_index).^2)./file_index);
            end        
            
            % store the time history of the acoustic intensity
            if record.I
                if file_index > 1
                    sensor_data.Ix(:, file_index-1) = 0.25*(sensor_data.ux_prev_t + ux_sgx(sensor_mask_index) + ux_sgx(sensor_mask_index_sgx)).*sensor_data.p_prev_t;
                end
            end
            
            % store the average acoustic intensity
            if record.I_avg
                if file_index > 1
                    if record.I
                        sensor_data.Ix_avg = sensor_data.Ix_avg + sensor_data.Ix(:, file_index-1);
                    else
                        sensor_data.Ix_avg = sensor_data.Ix_avg + 0.25*(sensor_data.ux_prev_t + ux_sgx(sensor_mask_index) + ux_sgx(sensor_mask_index_sgx)).*sensor_data.p_prev_t;
                    end
                end
            end
            
            % update the stored value of the particle velocity at the
            % previous time step if storing the intensity but not the time
            % varying particle velocity
            if (record.I || record.I_avg)
                sensor_data.ux_prev_t = ux_sgx(sensor_mask_index) + ux_sgx(sensor_mask_index_sgx);
                sensor_data.p_prev_t = p(sensor_mask_index);
            end
            
        % extract data from specified Cartesian coordinates using
        % linear interpolation   
        else

            % store the time history of the acoustic pressure
            if record.p
                sensor_data.p(:, file_index) = interp1(grid_x, p, sensor_x);
            end

            % store the maximum acoustic pressure
            if record.p_max
                if file_index == 1
                    sensor_data.p_max = interp1(grid_x, p, sensor_x);
                else
                    sensor_data.p_max = max(sensor_data.p_max, interp1(grid_x, p, sensor_x));
                end
            end     
            
            % store the minimum acoustic pressure
            if record.p_min
                if file_index == 1
                    sensor_data.p_min = interp1(grid_x, p, sensor_x);
                else
                    sensor_data.p_min = min(sensor_data.p_min, interp1(grid_x, p, sensor_x));
                end
            end               

            % store the rms acoustic pressure
            if record.p_rms
                sensor_data.p_rms = sqrt((sensor_data.p_rms.^2*(file_index - 1) + (interp1(grid_x, p, sensor_x)).^2)./file_index);
            end

            % store the time history of the particle velocity
            if record.u
                sensor_data.ux(:, file_index) = interp1(grid_x, ux_sgx, sensor_x);
            end

            % store the maximum particle velocity
            if record.u_max
                if file_index == 1
                    sensor_data.ux_max = interp1(grid_x, ux_sgx, sensor_x);
                else
                    sensor_data.ux_max = max(sensor_data.ux_max, interp1(grid_x, ux_sgx, sensor_x));
                end
            end  
            
            % store the minimum particle velocity
            if record.u_min
                if file_index == 1
                    sensor_data.ux_min = interp1(grid_x, ux_sgx, sensor_x);
                else
                    sensor_data.ux_min = min(sensor_data.ux_min, interp1(grid_x, ux_sgx, sensor_x));
                end
            end              

            % store the rms particle velocity
            if record.u_rms
                sensor_data.ux_rms = sqrt((sensor_data.ux_rms.^2*(file_index - 1) + (interp1(grid_x, ux_sgx, sensor_x)).^2)./file_index);
            end
            
            % store the time history of the acoustic intensity
            if record.I
                if file_index > 1
                    sensor_data.Ix(:, file_index-1) = 0.5*(sensor_data.ux_prev_t + interp1(grid_x, ux_sgx, sensor_x - kgrid.dx/2) ).* sensor_data.p_prev_t;
                end
            end
            
            % store the average acoustic intensity
            if record.I_avg
                if file_index > 1
                    if record.I
                        sensor_data.Ix_avg = sensor_data.Ix_avg + sensor_data.Ix(:, file_index-1);
                    else
                        sensor_data.Ix_avg = sensor_data.Ix_avg + 0.5*(sensor_data.ux_prev_t + interp1(grid_x, ux_sgx, sensor_x - kgrid.dx/2) ).* sensor_data.p_prev_t;
                    end
                end
            end
            
            % update the stored value of the particle velocity at the
            % previous time step if storing the intensity but not the time
            % varying particle velocity
            if (record.I || record.I_avg)
                sensor_data.ux_prev_t = interp1(grid_x, ux_sgx, sensor_x - kgrid.dx/2);
                sensor_data.p_prev_t = interp1(grid_x, p, sensor_x);
            end            
            
        end
        
% store the maximum acoustic pressure over all the grid
        % elements
        if record.p_max_all
            if file_index == 1
                sensor_data.p_max_all = p;
            else
                sensor_data.p_max_all = max(sensor_data.p_max_all, p);
            end
        end

        % store the minimum acoustic pressure over all the grid
        % elements
        if record.p_min_all
            if file_index == 1
                sensor_data.p_min_all = p;
            else
                sensor_data.p_min_all = min(sensor_data.p_min_all, p);
            end                
        end              
        
        % store the maximum particle velocity over all the grid
        % elements
        if record.u_max_all
            if file_index == 1
                sensor_data.ux_max_all = ux_sgx;
            else
                sensor_data.ux_max_all = max(sensor_data.ux_max_all, ux_sgx);                     
            end
        end        

        % store the minimum particle velocity over all the grid
        % elements
        if record.u_min_all
            if file_index == 1
                sensor_data.ux_min_all = ux_sgx;
            else
                sensor_data.ux_min_all = min(sensor_data.ux_min_all, ux_sgx);                    
            end
        end              
        
    end

    % estimate the time to run the simulation
    if t_index == ESTIMATE_SIM_TIME_STEPS
        disp(['  estimated simulation time ' scaleTime(etime(clock, loop_start_time)*index_end/t_index) '...']);
    end      
    
    % plot data if required
    if plot_sim && (rem(t_index, plot_freq) == 0 || t_index == 1 || t_index == index_end)  

        % update progress bar
        waitbar(t_index/length(t_array), pbar);
        drawnow;

        % ensure p is cast as a CPU variable and remove the PML from the
        % plot if required
        if strcmp(data_cast, 'gpuArray')
            p_plot = double(gather(p(x1:x2)));
        else
            p_plot = double(p(x1:x2));  
        end
                       
        % update plot
        if plot_scale_auto || plot_scale_log || t_index == 1
            
            % update plot scale if set to automatic or log
            if plot_scale_auto || plot_scale_log
                kspaceFirstOrder_adjustPlotScale;
            end
                        
            % replace entire plot
            if ~nonuniform_grid
                img_data = plot(kgrid.x(x1:x2)*scale, p_plot);
            else
                img_data = plot(kgrid.x_size*kgrid.xn(x1:x2)*scale, p_plot);
            end            
            
            % add display mask onto plot
            if ~(strcmp(display_mask, 'default') || strcmp(display_mask, 'off'))
                hold on;
                if ~nonuniform_grid
                    stairs(kgrid.x(x1:x2)*scale, display_mask(x1:x2).*(plot_scale(2) - plot_scale(1)) + plot_scale(1), 'k-');
                else
                    stairs(kgrid.xn(x1:x2)*scale, display_mask(x1:x2).*(plot_scale(2) - plot_scale(1)) + plot_scale(1), 'k-');
                end
                hold off
            end
            
            % set plot options
            xlabel(['x-position [' prefix 'm]']);
            if ~nonuniform_grid
                set(gca, 'YLim', plot_scale, 'XLim', kgrid.x([x1, x2])*scale);      
            else
                set(gca, 'YLim', plot_scale, 'XLim', kgrid.x_size*kgrid.xn([x1, x2])*scale);      
            end            
            
        else
            % just replace the y-data
            set(img_data, 'YData', p_plot);
        end
        
        % force plot update
        drawnow;
        
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


% assign the final intensity values
if use_sensor && ~time_rev && (record.I || record.I_avg)
    
    % update the values of the particle velocity at the next staggered time
    % point
    ux_sgx = pml_x_sgx .* (  pml_x_sgx.*ux_sgx - dt./rho0_sgx .* real(ifft(ddx_k .* shift_pos .* kappa .* p_k))  );
    
    % store the time history of the acoustic intensity
    if record.I
        if binary_sensor_mask
            sensor_data.Ix(:, file_index) = 0.25*(sensor_data.ux_prev_t + ux_sgx(sensor_mask_index) + ux_sgx(sensor_mask_index_sgx)).*sensor_data.p_prev_t;
        else
            sensor_data.Ix(:, file_index) = 0.5*(sensor_data.ux_prev_t + interp1(grid_x, ux_sgx, sensor_x - kgrid.dx/2) ).*sensor_data.p_prev_t;
        end
    end

    % store the average acoustic intensity
    if record.I_avg
        if record.I
            sensor_data.Ix_avg = sensor_data.Ix_avg + sensor_data.Ix(:, file_index);
        else
            if binary_sensor_mask
                sensor_data.Ix_avg = sensor_data.Ix_avg + 0.25*(sensor_data.ux_prev_t + ux_sgx(sensor_mask_index) + ux_sgx(sensor_mask_index_sgx)).*sensor_data.p_prev_t;
            else
                sensor_data.Ix_avg = sensor_data.Ix_avg + 0.5*(sensor_data.ux_prev_t + interp1(grid_x, ux_sgx, sensor_x - kgrid.dx/2) ).*sensor_data.p_prev_t;
            end
        end
    end
    
    % remove temporary variables
    sensor_data = rmfield(sensor_data, {'ux_prev_t', 'p_prev_t'});
    
end

% update command line status
disp(['  simulation completed in ' scaleTime(toc)]);

% =========================================================================
% CLEAN UP
% =========================================================================

% save the movie frames to disk
if record_movie
    kspaceFirstOrder_saveMovieFile;   
end

% clean up used figures
if plot_sim
    close(img);
    close(pbar);
    drawnow;
end

% reset the indexing variables to allow original grid size to be maintained
% (this is used to remove the PML from the user data if 'PMLInside' is set
% to false)
if (~plot_PML && PML_inside)
    x1 = x1 - PML_x_size;
    x2 = x2 + PML_x_size;
elseif (plot_PML && ~PML_inside)
    x1 = x1 + PML_x_size;
    x2 = x2 - PML_x_size;   
end

% save the final pressure field if in time reversal mode
if time_rev
    record.p_final = true;
end

% save the final acoustic pressure if required
if record.p_final
    sensor_data.p_final = p(x1:x2);
end

% save the final particle velocity if required
if record.u_final
    sensor_data.ux_final = ux_sgx(x1:x2);
end

% process the intensity data if required
if record.I_avg
    
    % divide the average values by the number of values recorded (the
    % variable num_recorded_time_points is set in createStorageVariables)
    sensor_data.Ix_avg = sensor_data.Ix_avg / num_recorded_time_points;
    
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

if ~use_sensor
    % if sensor is not used, return empty sensor data
    sensor_data = [];
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