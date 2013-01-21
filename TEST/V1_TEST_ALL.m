function V1_TEST_ALL
% Function to test the functionality of kspaceFirstOrder1D,
% kspaceFirstOrder2D, and kspaceFirstOrder3D. 
%
% ASCII generated using http://patorjk.com/software/taag/
%
% author: Bradley Treeby
% date: 4th September 2012
% last update: 4th October 2012

% -------------------------------------------------------------------------
% TEST MATRIX:
% -------------------------------------------------------------------------
% 1.1 Windows 7    2.1 1D   3.1 PML Inside   4.1 Even Grid Size  5.1 double
% 1.2 Windows XP   2.2 2D   3.2 PML Outside  4.2 Odd Grid Size   5.2 single    
% 1.3 Linux        2.3 3D                                        5.3 gsingle 
% 1.4 Mac                                                        5.4 GPUsingle
%                                                                5.5 gpuArray-single 
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% set the list of tests to run
% -------------------------------------------------------------------------

% 1: kspaceFirstOrder1D
% 2: kspaceFirstOrder2D
% 3: kspaceFirstOrder3D
TEST_DIM = 3;

% 1: Even grid size, PML inside
% 2: Even grid size, PML outside
% 3: Odd grid size, PML inside
% 4: Odd grid size, PML outside
TEST_TYPE = 2;

% overwrite test case if required (comment out to run all tests)
% -----------------------------------------------------------------------------------------------------
TEST_CASE = [...
% -----------------------------------------------------------------------------------------------------
%   LIN ABS HET SMH ALP SRC MNY DIR BIN P-U STR RCT FRQ RSI LOG PML LSC DIS PFQ PSC PLT LAY MOV TRV CPP   
%   01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25 
% ----------------------------------------------------------------------------------------------------- 
    0   0   0   0   0   6   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1; ...  
% -----------------------------------------------------------------------------------------------------
];

% run test
if exist('TEST_CASE', 'var')
    test_function(TEST_DIM, TEST_TYPE, TEST_CASE);
else
    test_function(TEST_DIM, TEST_TYPE);
end

end

% subfunction to run tests
function test_function(TEST_DIM, TEST_TYPE, TEST_CASE)
%#ok<*UNRCH>
%#ok<*NASGU>

% check for k-Wave Toolbox
if exist('getkWavePath', 'file') == 0
    addpath('/home/jaros/RSISE/K-Wave++/KSpaceFirstOrder3D_2.13/k-Wave Toolbox');
end

% get start time
start_time = clock;

% =========================================================================
% TEST OPTIONS
% =========================================================================

% tests to run
RUN_SOURCE_TESTS                        = true;
RUN_BIN_SENSOR_TESTS                    = true;
RUN_CART_SENSOR_LIN_INTERP_TESTS        = true;
RUN_CART_SENSOR_NEAREST_INTERP_TESTS    = true;
RUN_DISPLAY_TESTS                       = true;
RUN_TIME_REVERSAL_TESTS                 = true;

% force MATLAB k-Wave plotting to be off
FORCE_PLOT_OFF                          = false;

% what to do with the output data
SAVE_TEST_LOG                           = true;
SPLIT_LOG_AFTER_N_TESTS                 = 200;
CONTINUE_PAST_ERRORS                    = false;

% where to put all the output data
if ismac
    OUTPUT_FOLDER                       = '~/Documents/MATLAB/kWaveTests/';
elseif isunix
    OUTPUT_FOLDER                       = './';
else
    OUTPUT_FOLDER                       = 'C:\Users\Bradley\Documents\kWave Tests\';
end

% set to 'single' or 'gsingle' to speed up the matlab computations
DATA_CAST                               = 'off';  

% option to use non-uniform grid (NOT CURRENTLY USED)
USE_NONUNIFORM_GRID                     = false;

% set start index to skip tests on restart (set to 1 to run all tests)
TEST_CASE_START_INDEX                   = 1;

switch TEST_TYPE
    case 1
        GRID_SIZE_EVEN     = true;  % odd or even grid size
        PML_INSIDE         = true;  % location of the PML
    case 2
        GRID_SIZE_EVEN     = true;  % odd or even grid size
        PML_INSIDE         = false; % location of the PML
    case 3
        GRID_SIZE_EVEN     = false; % odd or even grid size
        PML_INSIDE         = true;  % location of the PML
    case 4
        GRID_SIZE_EVEN     = false; % odd or even grid size
        PML_INSIDE         = false; % location of the PML
end

% =========================================================================
% C++ TEST OPTIONS (3D simulations only)
% =========================================================================

% comparison with C++ code
RUN_CPP_COMPARISON_TESTS                = true;
PLOT_CPP_ERRORS                         = false;
SAVE_CPP_COMPARISON_PLOTS_TO_DISK       = false;
IMAGE_FOLDERNAME                        = 'CPP_COMPARISON_IMAGES/';

% switch to test MPI version instead of OpenMP
TEST_MPI                                = false;

% =========================================================================
% SIMULATION LITERALS
% =========================================================================

% set total number of grid points and the size of the perfectly matched
% layer (PML) 
if GRID_SIZE_EVEN
    if PML_INSIDE
        NX = 128;           % [grid points]
        NY = 64;            % [grid points]
        NZ = 32;            % [grid points]

        PML_X_SIZE = 20;    % [grid points]
        PML_Y_SIZE = 10;    % [grid points]
        PML_Z_SIZE = 10;    % [grid points]
    else
        NX = 128;           % [grid points]
        NY = 64;            % [grid points]
        NZ = 32;            % [grid points]

        PML_X_SIZE = 20;    % [grid points]
        PML_Y_SIZE = 10;    % [grid points]
        PML_Z_SIZE = 11;    % [grid points]               
    end
else
    if PML_INSIDE
        NX = 125;           % [grid points]
        NY = 63;            % [grid points]
        NZ = 35;            % [grid points]
        
        PML_X_SIZE = 20;    % [grid points]
        PML_Y_SIZE = 10;    % [grid points]
        PML_Z_SIZE = 10;    % [grid points]
    else
        NX = 135;           % [grid points]
        NY = 61;            % [grid points]
        NZ = 29;            % [grid points]
        
        PML_X_SIZE = 20;    % [grid points]
        PML_Y_SIZE = 10;    % [grid points]
        PML_Z_SIZE = 10;    % [grid points]
    end
end

% define the properties of the propagation medium
C0 = 1540;      % [m/s]
RHO0 = 1000;    % [kg/m^3]
ALPHA0 = 0.75;  % [dB/(MHz^Y cm)]
Y = 1.5;        % power law exponent
BONA = 6;       % parameter of nonlinearity
SC = 1.3;       % scale factor for heterogeneity

% =========================================================================
% CREATE TEST LOG
% =========================================================================

% get PC information
v = ver('matlab');
oper_sys = computer;
if isunix
    [~, comp] = unix('hostname', '-echo');
    comp = comp(1:end-1);
    comp_spacer = '_';
else
    comp = '';
    comp_spacer = '';
end

if SAVE_TEST_LOG 
    date_string = getDateString;
    diary([OUTPUT_FOLDER 'TESTLOG_' comp comp_spacer oper_sys '_MATLAB_' v.Version '_' v.Release '_' num2str(TEST_DIM) 'D_' DATA_CAST '_' date_string '.txt']);
end

% =========================================================================
% LIST OF SIMULATION OPTIONS
% =========================================================================

% -----------------------------------------------------------------------
% PARAMETER 1
%   0: Linear
%   1: Nonlinear
% PARAMETER 2
%   0: Lossless
%   1: Absorbing
% PARAMETER 3
%   0: Homogeneous
%   1: Heterogeneous c0 and rho0 only
%   2: Heterogeneous
% PARAMETER 4
%   0: 'Smooth' not set (defaults to [true, false, false])
%   1: 'Smooth', true
% PARAMETER 5
%   0: medium.alpha_mode not set
%   1: medium.alpha_mode = 'no_absorption'
%   2: medium.alpha_mode = 'no_dispersion'
% -----------------------------------------------------------------------
% PARAMETER 6
%   0: initial pressure
%   1: pressure source
%   2: velocity-x source
%   3: velocity-y source (2D and 3D only)
%   4: velocity-z source (3D only)
%   5: velocity-x/y/z source
%   6: transducer source (3D only)
% PARAMETER 7
%   0: single source
%   1: source many
% PARAMETER 8
%   0: additive source condition
%   1: dirichlet source condition
% -----------------------------------------------------------------------
% PARAMETER 9
%   0: binary sensor mask
%   1: Cartesian sensor mask, linear interpolation
%   2: Cartesian sensor mask, nearest neighbour interpolation
%   3: binary sensor mask with sensor directivity (2D only)
% PARAMETER 10
%   0: record only pressure (no sensor.record input)
%   1: record pressure and velocity
%   2: record max and rms pressure and velocity
%   3: record everything
% PARAMETER 11 (3D only)
%   0: 'StreamToDisk', false
%   1: 'StreamToDisk', true
%   2: 'StreamToDisk', 50
% PARAMETER 12
%   0: 'DataRecast', false
%   1: 'DataRecast', true
% PARAMETER 13
%   0: sensor.frequency_response not set
%   1: sensor.frequency_response = [0.5e6, 50]
% PARAMETER 14
%   0: sensor.record_start_index not set
%   1: sensor.record_start_index = 200;
% -----------------------------------------------------------------------
% PARAMETER 15
%   0: 'CreateLog', false (default)
%   1: 'CreateLog', true
% PARAMETER 16
%   0: 'PlotPML', true (default)
%   1: 'PlotPML', false (default)
% PARAMETER 17
%   0: 'LogScale', false (default)
%   1: 'LogScale', true
% PARAMETER 18
%   0: 'DisplayMask' not set (defaults to sensor.mask)
%   1: 'DisplayMask', 'off'
%   2: 'DisplayMask', makeBall
% PARAMETER 19
%   0: 'PlotFreq' not set (defaults to 10) 
%   1: 'PlotFreq', 30
% PARAMETER 20
%   0: 'PlotScale' set to -[p0, p0]
%   1: 'PlotScale' set to 'auto'
% PARAMETER 21
%   0: 'PlotSim' not set (defaults to true) 
%   1: 'PlotSim', false
%   2: 'MeshPlot', true (2D only)
% PARAMETER 22
%   0: 'PlotLayout' not set (defaults to false)
%   1: 'PlotLayout', true
% PARAMETER 23
%   0: 'RecordMovie' not set (defaults to false)
%   1: 'RecordMovie', true
%   2: 'RecordMovie', true / 'MovieType', 'image' (2D only)
% -----------------------------------------------------------------------
% PARAMETER 24
%   0: run only forward simulation
%   1: run time reversal simulation
%   2: run time reversal including attenuation compensation
% -----------------------------------------------------------------------
% PARAMETER 25 (3D only)
%   0: run MATLAB version only
%   1: run C++ version and compare outputs
% -----------------------------------------------------------------------

% set the number of variations for each parameters
variations = [2, 2, 3, 2, 3,...
    7, 2, 2,...
    4, 4, 3, 2, 2,...
    2, 2, 2, 3, 2, 2, 2, 2, 3,...
    3,...
    2];

% =========================================================================
% SET TEST CASES
% =========================================================================

% set pseudonym for running C++ comparisons - for tests that make sense to
% compare with the C++ code, option 25 is set to C, otherwise option 25 is
% set to 0
RUN_CPP_COMPARISON_TESTS = RUN_CPP_COMPARISON_TESTS && (TEST_DIM == 3);
C = double(RUN_CPP_COMPARISON_TESTS);

test_cases = [];
if RUN_SOURCE_TESTS
    test_cases = [test_cases;...
% -----------------------------------------------------------------------------------------------------
%   A: TEST SOURCE CONDITIONS IN HOMOGENEOUS AND HETEROGENEOUS MEDIA
% -----------------------------------------------------------------------------------------------------
%   LIN ABS HET SMH ALP SRC MNY DIR BIN P-U STR RCT FRQ RSI LOG PML LSC DIS PFQ PSC PLT LAY MOV TRV CPP   
%   01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25 
% -----------------------------------------------------------------------------------------------------
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   3   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   4   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   5   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   6   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   3   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   4   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   5   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   6   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   2   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   3   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   4   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   5   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   2   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   3   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   4   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   5   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   6   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...      
    0   0   0   0   0   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   2   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   3   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   4   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   5   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   2   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   3   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   4   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   5   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   2   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   3   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   4   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   5   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   2   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   3   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   4   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   5   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ... 
    1   0   1   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...     
    0   1   1   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...     
    1   1   1   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...     
    1   0   2   0   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...     
    0   1   2   0   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...     
    1   1   2   0   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ... 
    1   0   0   0   0   3   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...     
    0   1   0   0   0   3   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...     
    1   1   0   0   0   3   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...    
% -----------------------------------------------------------------------------------------------------
    ];
end

if RUN_BIN_SENSOR_TESTS
    test_cases = [test_cases;...
% -----------------------------------------------------------------------------------------------------
%   B: TEST BINARY SENSOR CONDITIONS
% -----------------------------------------------------------------------------------------------------
%   LIN ABS HET SMH ALP SRC MNY DIR BIN P-U STR RCT FRQ RSI LOG PML LSC DIS PFQ PSC PLT LAY MOV TRV CPP   
%   01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25 
% ----------------------------------------------------------------------------------------------------- 
    0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   0   0   0   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ... 
    0   0   0   0   0   0   0   0   0   3   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ... 
    0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   C; ...      
    0   0   0   0   0   0   0   0   0   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   0   0   0   0   3   0   1   0   0   0   0   0   0   0   0   0   0   0   0   C; ... 
    0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0; ...     
    0   0   0   0   0   0   0   0   0   0   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0; ...       
    0   0   0   0   0   0   0   0   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0; ...           
    0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   0   0   0   0   1   0   0   0   1   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   0   0   0   0   2   0   0   0   1   0   0   0   0   0   0   0   0   0   0   C; ... 
    0   0   0   0   0   0   0   0   0   3   0   0   0   1   0   0   0   0   0   0   0   0   0   0   C; ... 
    0   0   0   0   0   0   0   0   0   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   0   0   0   3   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
% -----------------------------------------------------------------------------------------------------
    ];
end        
    
if RUN_CART_SENSOR_LIN_INTERP_TESTS
    test_cases = [test_cases;...
% -----------------------------------------------------------------------------------------------------
%   C: TEST CARTESIAN SENSOR CONDITIONS WITH LINEAR INTERPOLATION
% -----------------------------------------------------------------------------------------------------
%   LIN ABS HET SMH ALP SRC MNY DIR BIN P-U STR RCT FRQ RSI LOG PML LSC DIS PFQ PSC PLT LAY MOV TRV CPP   
%   01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25 
% -----------------------------------------------------------------------------------------------------
    0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ...      
    0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   1   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   1   3   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   1   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0; ...      
    0   0   0   0   0   0   0   0   1   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   1   3   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   1   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   1   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0; ...     
    0   0   0   0   0   0   0   0   1   0   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0; ...       
    0   0   0   0   0   0   0   0   1   0   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0; ...           
    0   0   0   0   0   0   0   0   1   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   1   1   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   1   2   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   1   3   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   1   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ...   
% -----------------------------------------------------------------------------------------------------
    ];
end 

if RUN_CART_SENSOR_NEAREST_INTERP_TESTS
    test_cases = [test_cases;...
% -----------------------------------------------------------------------------------------------------
%   D: TEST CARTESIAN SENSOR CONDITIONS WITH NEAREST NEIGHBOUR INTERPOLATION
% -----------------------------------------------------------------------------------------------------
%   LIN ABS HET SMH ALP SRC MNY DIR BIN P-U STR RCT FRQ RSI LOG PML LSC DIS PFQ PSC PLT LAY MOV TRV CPP   
%   01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25 
% -----------------------------------------------------------------------------------------------------
    0   0   0   0   0   0   0   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ...      
    0   0   0   0   0   0   0   0   2   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   2   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   2   3   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   2   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   2   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0; ...      
    0   0   0   0   0   0   0   0   2   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   2   3   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   2   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   2   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0; ...     
    0   0   0   0   0   0   0   0   2   0   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0; ...       
    0   0   0   0   0   0   0   0   2   0   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0; ...           
    0   0   0   0   0   0   0   0   2   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   2   1   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   2   2   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   2   3   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   2   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ...  
% -----------------------------------------------------------------------------------------------------
    ];
end 

if RUN_DISPLAY_TESTS
    test_cases = [test_cases;...
% -----------------------------------------------------------------------------------------------------
%   E: TEST DISPLAY CONDITIONS
% -----------------------------------------------------------------------------------------------------
%   LIN ABS HET SMH ALP SRC MNY DIR BIN P-U STR RCT FRQ RSI LOG PML LSC DIS PFQ PSC PLT LAY MOV TRV CPP   
%   01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25 
% -----------------------------------------------------------------------------------------------------
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0; ...
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   1   0   0   0   0   0; ...
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0; ...
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   1   0   0   0   0   0; ...
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   2   0   0   0   0   0   0   0; ...
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   2   0   1   0   0   0   0   0; ...
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   1   0   0   0; ...
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   1   0   0; ...
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   1   0   0   1   0   0; ...
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   2   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   1   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   2   0   1   0   0   0   0   0; ...      
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0; ...    
    0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0; ...  
    0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   1   0   0   0; ...  
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0; ...  
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0; ...      
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   1   0   0; ...    
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   2   0   0   0   0   1   0   0; ...     
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   2   0   1   0   0   1   0   0; ... 
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0; ...     
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   2   0   0   0   0; ...         
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   2   0   0   0   0; ...         
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   2   0   0   0   0; ...             
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   1   2   0   0   0   0; ...             
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   2   0   1   0   0; ...         
% -----------------------------------------------------------------------------------------------------
    ]; 
end 

if RUN_TIME_REVERSAL_TESTS
    test_cases = [test_cases;...
% -----------------------------------------------------------------------------------------------------
%   E: TEST DISPLAY CONDITIONS
% -----------------------------------------------------------------------------------------------------
%   LIN ABS HET SMH ALP SRC MNY DIR BIN P-U STR RCT FRQ RSI LOG PML LSC DIS PFQ PSC PLT LAY MOV TRV CPP   
%   01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25 
% -----------------------------------------------------------------------------------------------------
    0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   1   C; ...     
    0   0   0   0   0   0   0   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0   1   0; ... 
    0   0   0   0   0   0   0   0   2   0   0   1   0   0   0   0   0   0   0   0   0   0   0   1   0; ...   
    0   0   1   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   1   C; ...     
    0   0   1   0   0   0   0   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0   1   0; ... 
    0   0   1   0   0   0   0   0   2   0   0   1   0   0   0   0   0   0   0   0   0   0   0   1   0; ...       
    0   1   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   2   0; ...     
    0   1   0   0   0   0   0   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0   2   0; ... 
    0   1   0   0   0   0   0   0   2   0   0   1   0   0   0   0   0   0   0   0   0   0   0   2   0; ...   
    0   1   1   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   2   0; ...     
    0   1   1   0   0   0   0   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0   2   0; ... 
    0   1   1   0   0   0   0   0   2   0   0   1   0   0   0   0   0   0   0   0   0   0   0   2   0; ...       
% -----------------------------------------------------------------------------------------------------
    ]; 
end 

% =========================================================================
% REMOVE TESTS THAT DON'T MAKE SENSE FOR PARTICULAR DIMENSIONS
% =========================================================================

% remove transducer source tests
if TEST_DIM ~= 3
    test_cases(find(test_cases(:, 6) == 6), :) = []; %#ok<*FNDSB>
end

% remove velocity-x/yz source tests
if TEST_DIM == 1
    test_cases(find(test_cases(:, 6) == 5), :) = [];
end

% remove velocity-z source tests
if TEST_DIM ~= 3
    test_cases(find(test_cases(:, 6) == 4), :) = [];
end

% remove velocity-y source tests
if TEST_DIM == 1
    test_cases(find(test_cases(:, 6) == 3), :) = [];
end

% remove directivity tests
if TEST_DIM ~= 2 || strcmp(DATA_CAST, 'gpuArray-single') || strcmp(DATA_CAST, 'GPUsingle')
    test_cases(find(test_cases(:, 9) == 3), :) = [];
end

% remove Cart interp nearest tests without recast
if ~(strcmp(DATA_CAST, 'off') || strcmp(DATA_CAST, 'single'))
    test_cases(find( (test_cases(:, 9) == 2) & (test_cases(:, 12) == 0)), :) = [];
end

% remove Cart sensor mask test (interp1 not supported by PCT or GPUmat)
if TEST_DIM == 1 && (strcmp(DATA_CAST, 'GPUsingle') || strcmp(DATA_CAST, 'gpuArray-single'))
    test_cases(find(test_cases(:, 9) == 1), :) = [];
end

% remove tests recording p_max and u_max (max not supported by GPUmat)
if strcmp(DATA_CAST, 'GPUsingle')
    test_cases(find(test_cases(:, 10) == 2), :) = [];
    test_cases(find(test_cases(:, 10) == 3), :) = [];
end

% remove stream to disk tests
if TEST_DIM ~= 3
    test_cases(find(test_cases(:, 11) ~= 0), :) = [];
end

% remove frequency response tests for GPU tests without recasting
if ~(strcmp(DATA_CAST, 'off') || strcmp(DATA_CAST, 'single'))
    test_cases(find( (test_cases(:, 13) == 1) & (test_cases(:, 12) == 0)), :) = [];
end

% remove mesh plot tests
if TEST_DIM ~= 2
    test_cases(find(test_cases(:, 21) == 2), :) = [];
end

% remove plot layout tests
if strcmp(DATA_CAST, 'GPUsingle')
    test_cases(find(test_cases(:, 22) == 1), :) = [];
end

% remove movie image tests
if TEST_DIM ~= 2
    test_cases(find(test_cases(:, 23) == 2), :) = [];
end

% =========================================================================
% REPLACE WITH A SINGLE TEST CASE IF DESIRED
% =========================================================================

if nargin == 3 && ~isempty(TEST_CASE)
    test_cases = TEST_CASE; 
end

% =========================================================================
% COMMAND LINE UPDATE
% =========================================================================
clc;
disp('   ');
disp('-----------------------------------------------------------------------------------------------');
disp('                 _      __        __                _____         _            ');
disp('                | | __  \ \      / /_ ___   _____  |_   _|__  ___| |_ ___ _ __ ');
disp('                | |/ /___\ \ /\ / / _` \ \ / / _ \   | |/ _ \/ __| __/ _ \ ''__|');
disp('                |   <_____\ V  V / (_| |\ V /  __/   | |  __/\__ \ ||  __/ |   ');
disp('                |_|\_\     \_/\_/ \__,_| \_/ \___|   |_|\___||___/\__\___|_|   ');
disp('  ');                                                                
disp('-----------------------------------------------------------------------------------------------');
disp('  ');
disp(['START DATE:                                      ' datestr(start_time)]);
if isunix
disp(['COMPUTER:                                        ' comp]);
end
disp(['O/S:                                             ' oper_sys]);
disp(['MATLAB VERSION:                                  ' v.Version ' ' v.Release]);
disp(['DATA CAST:                                       ' DATA_CAST]);
disp(['TEST DIM:                                        ' num2str(TEST_DIM)]);
switch TEST_DIM
    case 1
disp(['GRID SIZE:                                       ' num2str(NX)]);
disp(['PML SIZE:                                        ' num2str(PML_X_SIZE)]);
    case 2
disp(['GRID SIZE:                                       ' num2str(NX) ' by ' num2str(NY)]);
disp(['PML SIZE:                                        ' num2str(PML_X_SIZE) ' by ' num2str(PML_Y_SIZE)]);
    case 3
disp(['GRID SIZE:                                       ' num2str(NX) ' by ' num2str(NY) ' by ' num2str(NZ)]);
disp(['PML SIZE:                                        ' num2str(PML_X_SIZE) ' by ' num2str(PML_Y_SIZE) ' by ' num2str(PML_Z_SIZE)]);
end
disp(['PML INSIDE:                                      ' num2str(PML_INSIDE)]);
disp(['NONUNIFORM GRID:                                 ' num2str(USE_NONUNIFORM_GRID)]);
disp('  ');
disp(['RUN SOURCE TESTS:                                ' num2str(RUN_SOURCE_TESTS)]);
disp(['RUN BININARY SENSOR TESTS:                       ' num2str(RUN_BIN_SENSOR_TESTS)]);
disp(['RUN CARTESIAN SENSOR TESTS (LIN INTERP):         ' num2str(RUN_CART_SENSOR_LIN_INTERP_TESTS)]);
disp(['RUN CARTESIAN SENSOR TESTS (NN INTERP):          ' num2str(RUN_CART_SENSOR_NEAREST_INTERP_TESTS)]);
disp(['RUN DISPLAY TESTS:                               ' num2str(RUN_DISPLAY_TESTS)]);
disp(['RUN TIME REVERSAL TESTS:                         ' num2str(RUN_TIME_REVERSAL_TESTS)]);
disp(['RUN C++ COMPARISON TESTS:                        ' num2str(RUN_CPP_COMPARISON_TESTS)]);
disp('  ');
disp('  ');
disp(['The total number of possibe test combinations is: ' num2str(prod(variations))]);
disp(['This would take an estimated ' scaleTime(prod(variations)*60) ' to test at 1 min per simulation']);
disp(['The number of tested combinations is: ' num2str(size(test_cases, 1))]);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% reset the total number of errors
number_errors = 0;
number_cpp_errors = 0;
location_errors = [];
location_cpp_errors = [];

% extract parameters from test cast settings
for comp_index = TEST_CASE_START_INDEX:size(test_cases, 1)

    disp('  ');
    disp('---------------------------------------------------------------------------------------------------');
    disp(['INDEX OF CURRENT TEST = ' num2str(comp_index)])
    disp('---------------------------------------------------------------------------------------------------');
    disp('LIN ABS HET SMH ALP SRC MNY DIR BIN P-U STR RCT FRQ RSI LOG PML LSC DIS PFQ PSC PLT LAY MOV TRV CPP');
    disp('01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25');
    disp('---------------------------------------------------------------------------------------------------');
    disp(num2str(test_cases(comp_index, :), '%1d   '));
    disp('---------------------------------------------------------------------------------------------------');
    disp(' ');
    
    % -------------------------------------------------------------------
    
    % PARAMETER 1
    NONLINEAR = (test_cases(comp_index, 1) == 1);
    
    % PARAMETER 2
    ABSORBING = (test_cases(comp_index, 2) == 1);

    % PARAMETER 3
    if test_cases(comp_index, 3) == 0
        HETEROGENEOUS_RHO0 = false;
        HETEROGENEOUS_C0 = false;
        HETEROGENEOUS_BONA = false;      % not used if NONLINEAR = false
        HETEROGENEOUS_ALPHA0 = false;    % not used if ABSORBING = false
    elseif test_cases(comp_index, 3) == 1
        HETEROGENEOUS_RHO0 = true;
        HETEROGENEOUS_C0 = true;
        HETEROGENEOUS_BONA = false;      % not used if NONLINEAR = false
        HETEROGENEOUS_ALPHA0 = false;    % not used if ABSORBING = false
    elseif test_cases(comp_index, 3) == 2
        HETEROGENEOUS_RHO0 = true;
        HETEROGENEOUS_C0 = true;
        HETEROGENEOUS_BONA = true;       % not used if NONLINEAR = false
        HETEROGENEOUS_ALPHA0 = true;     % not used if ABSORBING = false
    else
        error('unknown setting');
    end

    % PARAMETER 4
    SMOOTH = (test_cases(comp_index, 4) == 1);
    
    % PARAMETER 5
    switch test_cases(comp_index, 5)
        case 0
            ALPHA_MODE = [];
        case 1
            ALPHA_MODE = 'no_absorption';
        case 2
            ALPHA_MODE = 'no_dispersion';
    end    
    
	% -------------------------------------------------------------------
    
    % PARAMETER 6
    SOURCE_TYPE = test_cases(comp_index, 6);

    % PARAMETER 7
    SOURCE_MANY = test_cases(comp_index, 7);
    
    % PARAMETER 8
    SOURCE_DIRICHLET = test_cases(comp_index, 8);
    
    % -------------------------------------------------------------------
        
    % PARAMETER 9
    switch test_cases(comp_index, 9)
        case 0
            BINARY_SENSOR_MASK = true;
            SENSOR_DIRECTIVITY = false;
        case 1
            BINARY_SENSOR_MASK = false;
            CART_INTERP = 'linear';
        case 2
            BINARY_SENSOR_MASK = false;
            CART_INTERP = 'nearest';
        case 3
            BINARY_SENSOR_MASK = true;
            SENSOR_DIRECTIVITY = true;
    end
    
    % PARAMETER 10
    switch test_cases(comp_index, 10)
        case 0
            SENSOR_RECORD_FIELDS = [];
        case 1
            SENSOR_RECORD_FIELDS = {'p', 'u'};
        case 2
            SENSOR_RECORD_FIELDS = {'p_max', 'p_rms', 'u_max', 'u_rms'};
        case 3
            if TEST_MPI
                SENSOR_RECORD_FIELDS = {'p', 'u', 'p_max', 'p_rms', 'u_max', 'u_rms', 'p_final', 'u_final'};
            else
                SENSOR_RECORD_FIELDS = {'p', 'u', 'p_max', 'p_rms', 'u_max', 'u_rms', 'p_final', 'u_final', 'I', 'I_avg'};
            end
    end
    
    % PARAMETER 11
    switch test_cases(comp_index, 11)
        case 0
            STREAM_TO_DISK = [];
        case 1
            STREAM_TO_DISK = true;
        case 2
            STREAM_TO_DISK = 20;
    end
    
    % PARAMETER 12
    DATA_RECAST = (test_cases(comp_index, 12) == 1);
        
    % PARAMETER 13
    switch test_cases(comp_index, 13)
        case 0
            FREQUENCY_RESPONSE = [];
        case 1
            FREQUENCY_RESPONSE = [0.5e6, 50];
    end
    
    % PARAMETER 14
    switch test_cases(comp_index, 14)
        case 0
            RECORD_START_INDEX = 1;
        case 1
            RECORD_START_INDEX = 200;
    end
    
    % -------------------------------------------------------------------
    
    % PARAMETER 15
    switch test_cases(comp_index, 15)
        case 0
            CREATE_LOG = [];
        case 1
            CREATE_LOG = true;
    end
    
    % PARAMETER 16
    switch test_cases(comp_index, 16)
        case 0
            PLOT_PML = [];
        case 1
            PLOT_PML = false;
    end
    
    % PARAMETER 17
    switch test_cases(comp_index, 17)
        case 0
            PLOT_LOG_SCALE = [];
        case 1
            PLOT_LOG_SCALE = false;
    end
    
    % PARAMETER 18
    switch test_cases(comp_index, 18)
        case 0
            DISPLAY_MASK = [];
        case 1
            DISPLAY_MASK = 'off';
        case 2
            DISPLAY_MASK = 'ball';
    end
    
    % PARAMETER 19
    switch test_cases(comp_index, 19)
        case 0
            PLOT_FREQ = [];
        case 1
            PLOT_FREQ = 30;
    end    
    
    % PARAMETER 20
    switch test_cases(comp_index, 20)
        case 0
            PLOT_SCALE = [];
        case 1
            PLOT_SCALE = 'auto';
    end 
    
    % PARAMETER 21
    switch test_cases(comp_index, 21)
        case 0
            PLOT_SIM = [];
            MESH_PLOT = [];
        case 1
            PLOT_SIM = false;
            MESH_PLOT = [];
        case 2
            PLOT_SIM = [];
            MESH_PLOT = true;
            
    end 
    
    % PARAMETER 22
    switch test_cases(comp_index, 22)
        case 0
            PLOT_LAYOUT = [];
        case 1
            PLOT_LAYOUT = true;
    end     
        
    % PARAMETER 23
    switch test_cases(comp_index, 23)
        case 0
            SAVE_MOVIE = [];
            MOVIE_TYPE_IMAGE = [];
        case 1
            SAVE_MOVIE = true;
            MOVIE_TYPE_IMAGE = [];
        case 2
            SAVE_MOVIE = true;
            MOVIE_TYPE_IMAGE = true;
    end  
    
    % -------------------------------------------------------------------
    
    % PARAMETER 24
    RUN_TIME_REVERSAL = test_cases(comp_index, 24);
    
    % -------------------------------------------------------------------
    
    % PARAMETER 25
    RUN_CPP_CODE = test_cases(comp_index, 25);
       
    % -------------------------------------------------------------------
    
    % run simulation
    try
        run_simulation;
    catch ME
        if CONTINUE_PAST_ERRORS
            disp(' ');
            disp('  SIMULATION ERROR:');
            disp('  -----------------');
            disp(['  message identifier: ' ME.identifier]);
            disp(['  occured in: ' ME.stack(1).name ' at line ' num2str(ME.stack(1).line)]);
            disp(['  message: ' ME.message]);
            disp(' ');
            number_errors = number_errors + 1;
            location_errors = [location_errors, comp_index];
        else
            diary off;
            rethrow(ME);
        end
    end
    
    % split test log in parts (some problems with really long log files)
    if SAVE_TEST_LOG && ~rem(comp_index, SPLIT_LOG_AFTER_N_TESTS)
        diary off;
        diary([OUTPUT_FOLDER 'TESTLOG_' comp comp_spacer oper_sys '_MATLAB_' v.Version '_' v.Release '_' DATA_CAST '_' date_string '_PART' num2str(floor(comp_index / SPLIT_LOG_AFTER_N_TESTS) + 1) '.txt']);        
    end

end

disp(' ');
disp('---------------------------------------------------------------------------------------------------');
disp(['TOTAL RUNTIME ERRORS:        ' num2str(number_errors)])
disp(['INDEX OF RUNTIME ERRORS:     ' num2str(location_errors)]);
if RUN_CPP_COMPARISON_TESTS
disp(['C++ SIMS WITH L_INF > 1e-5:  ' num2str(number_cpp_errors)]);
disp(['INDEX OF C++ ERRORS:         ' num2str(location_cpp_errors)]);
end
disp(['ELAPSED TIME:                ' scaleTime(etime(clock, start_time))]);
disp('---------------------------------------------------------------------------------------------------');

% switch off log
if SAVE_TEST_LOG 
    diary off;
end

% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
% ========================================================================
% ------------------------------------------------------------------------
% ========================================================================
% \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

function run_simulation()
% Nested subfunction to set options and run simulations. For 3D
% simulations, the outputs from the MATLAB and C++ version of the code are
% compared. (Note, nested subfunctions can see all variables created
% above.) 

% create empty input arguments;
input_args = {};

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% create heterogeneous region
switch TEST_DIM
    case 1
        heterog = ones(NX, 1);
        heterog(round(NX/2) - 9:round(NX/2) + 10) = 1;
    case 2
        heterog = makeDisc(NX, NY, round(NX/2), round(NY/2), 10);
    case 3
        heterog = makeBall(NX, NY, NZ, round(NX/2), round(NY/2), round(NZ/2), 10);
end

% set density
if HETEROGENEOUS_RHO0
    medium.density = RHO0*ones(size(heterog));
    medium.density(heterog == 1) = RHO0*SC;
else
    medium.density = RHO0; 
end

% set sound speed
if HETEROGENEOUS_C0
    medium.sound_speed = C0*ones(size(heterog));
    medium.sound_speed(heterog == 1) = C0*SC;
else
    medium.sound_speed = C0;
end
    
% set BONA
if NONLINEAR
    if HETEROGENEOUS_BONA
        medium.BonA = BONA*ones(size(heterog));
        medium.BonA(heterog == 1) = BONA*SC;
    else
        medium.BonA = BONA;
    end
end

% set absorption terms
if ABSORBING
    if HETEROGENEOUS_ALPHA0
        medium.alpha_coeff  = ALPHA0*ones(size(heterog));
        medium.alpha_coeff(heterog == 1) = ALPHA0*SC;
        medium.alpha_power = Y;
    else
        medium.alpha_coeff = ALPHA0;      
        medium.alpha_power = Y;
    end
    
    % set absorption mode
    if ~isempty(ALPHA_MODE)
        medium.alpha_mode = ALPHA_MODE;
    end
end

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set desired grid size in the x-direction not including the PML
x = 40e-3;                  % [m]

% calculate the spacing between the grid points
dx = x/NX;                  % [m]
dy = dx;                    % [m]
dz = dx;                    % [m]

% create the k-space grid
switch TEST_DIM
    case 1
        kgrid = makeGrid(NX, dx);
    case 2
        kgrid = makeGrid(NX, dx, NY, dy);
    case 3
        kgrid = makeGrid(NX, dx, NY, dy, NZ, dz);
end

% =========================================================================
% DEFINE NONUNIFORM GRID SETTINGS
% =========================================================================

% if USE_NONUNIFORM_GRID
%     
%    ...
%     
% end

% =========================================================================
% DEFINE THE TIME ARRAY
% =========================================================================

% create the time array
t_end = 15e-6;                  % [s]
kgrid.t_array = makeTime(kgrid, medium.sound_speed, [], t_end);

% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = 0.5e6;    	% [Hz]
tone_burst_cycles = 5;

% create the input signal using toneBurst 
input_signal = source_strength*toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% =========================================================================
% DEFINE THE SOURCE CONDITIONS
% =========================================================================

% 0: initial pressure
% 1: pressure source
% 2: velocity-x source
% 3: velocity-y source
% 4: velocity-z source
% 5: velocity-x/y/z source
% 6: transducer source

switch SOURCE_TYPE
    case 0
        
        % -----------------------------
        % OPTION 0: Initial Pressure
        % ----------------------------- 
        
        % create ball shaped initial pressure distribution
        source_radius = 5;
        switch TEST_DIM
            case 1
                source.p0 = ones(NX, 1);
                source.p0(round(NX/2) - source_radius + 1:round(NX/2) + source_radius) = source_strength;
            case 2
                source.p0 = source_strength*makeDisc(NX, NY, round(NX/2), round(NY/2), source_radius);
            case 3
                source.p0 = source_strength*makeBall(NX, NY, NZ, round(NX/2), round(NY/2), round(NZ/2), source_radius);
        end
        
    case 1
        
        % -----------------------------
        % OPTION 1: Pressure Source
        % -----------------------------
        
        % set the source mask to be a small rectangle
        source_radius = 15;
        switch TEST_DIM
            case 1
                source.p_mask = zeros(NX, 1);
                source.p_mask(PML_X_SIZE + 1) = 1;
                source.p_mask(PML_X_SIZE + source_radius) = 1;
            case 2
                source.p_mask = zeros(NX, NY);
                source.p_mask(PML_X_SIZE + 1, round(NY/2) - source_radius + 1:round(NY/2) + source_radius) = 1;
            case 3
                source.p_mask = zeros(NX, NY, NZ);
                source.p_mask(PML_X_SIZE + 1, round(NY/2) - source_radius + 1:round(NY/2) + source_radius, round(NZ/2) - round(source_radius/2) + 1:round(NZ/2) + round(source_radius/2)) = 1;
        end

        % assign the source term
        if SOURCE_MANY
            switch TEST_DIM
                case 1
                    focus_position = 0;
                    source.p = focus(kgrid, input_signal, source.p_mask, focus_position, C0);
                case 2
                    focus_position = [-5*dx, 0];
                    source.p = focus(kgrid, input_signal, source.p_mask, focus_position, C0);
                case 3
                    focus_position = [0, 0, 0];
                    source.p = focus(kgrid, input_signal, source.p_mask, focus_position, C0);
            end
        else
            source.p = input_signal;
        end
        
    case 2
        
        % -----------------------------
        % OPTION 2: Velocity X Source
        % -----------------------------
        
        % set the source mask to be a small rectangle
        source_radius = 10;
        switch TEST_DIM
            case 1
                source.u_mask = zeros(NX, 1);
                source.u_mask(PML_X_SIZE + 1) = 1;
                source.u_mask(PML_X_SIZE + source_radius) = 1;
            case 2
                source.u_mask = zeros(NX, NY);
                source.u_mask(PML_X_SIZE + 1, round(NY/2) - source_radius + 1:round(NY/2) + source_radius) = 1;
            case 3
                source.u_mask = zeros(NX, NY, NZ);
                source.u_mask(PML_X_SIZE + 1, round(NY/2) - source_radius + 1:round(NY/2) + source_radius, round(NZ/2) - source_radius/2 + 1:round(NZ/2) + source_radius/2) = 1;
        end

        % assign the source term scaled by the impedance
        if SOURCE_MANY
            switch TEST_DIM
                case 1
                    focus_position = 0;
                    source.ux = focus(kgrid, input_signal./(C0*RHO0), source.u_mask, focus_position, C0);
                case 2
                    focus_position = [-5*dx, 0];
                    source.ux = focus(kgrid, input_signal./(C0*RHO0), source.u_mask, focus_position, C0);
                case 3
                    focus_position = [0, 0, 0];
                    source.ux = focus(kgrid, input_signal./(C0*RHO0), source.u_mask, focus_position, C0);
            end
        else
            source.ux = input_signal./(C0*RHO0);
        end
        
    case 3
        
        % -----------------------------
        % OPTION 3: Velocity Y Source
        % -----------------------------
        
        % set the source mask to be a small rectangle
        source_radius = 10;
        switch TEST_DIM
            case 2
                source.u_mask = zeros(NX, NY);
                source.u_mask(round(NX/2) - source_radius + 1:round(NX/2) + source_radius, PML_Y_SIZE + 1) = 1;
            case 3
                source.u_mask = zeros(NX, NY, NZ);
                source.u_mask(round(NX/2) - source_radius + 1:round(NX/2) + source_radius, PML_Y_SIZE + 1, round(NZ/2) - source_radius + 1:round(NZ/2) + source_radius) = 1;
            otherwise
                error('source.uy not supported in 1D');
        end

        % assign the source term scaled by the impedance
        if SOURCE_MANY
            switch TEST_DIM
                case 1
                    focus_position = 0;
                    source.uy = focus(kgrid, input_signal./(C0*RHO0), source.u_mask, focus_position, C0);
                case 2
                    focus_position = [-5*dx, 0];
                    source.uy = focus(kgrid, input_signal./(C0*RHO0), source.u_mask, focus_position, C0);
                case 3
                    focus_position = [0, 0, 0];
                    source.uy = focus(kgrid, input_signal./(C0*RHO0), source.u_mask, focus_position, C0);
            end
        else
            source.uy = input_signal./(C0*RHO0);
        end
        
    case 4
        
        % -----------------------------
        % OPTION 4: Velocity Z Source
        % -----------------------------        
        
        % set the source mask to be a small rectangle
        source_radius = 10;
        switch TEST_DIM
            case 3
                source.u_mask = zeros(NX, NY, NZ);
                source.u_mask(round(NX/2) - source_radius + 1:round(NX/2) + source_radius, round(NY/2) - source_radius + 1:round(NY/2) + source_radius, PML_Z_SIZE + 1) = 1;
            otherwise
                error('source.uz only supported in 3D');
        end

        % assign the source term scaled by the impedance
        if SOURCE_MANY
            switch TEST_DIM
                case 1
                    focus_position = 0;
                    source.uz = focus(kgrid, input_signal./(C0*RHO0), source.u_mask, focus_position, C0);
                case 2
                    focus_position = [-5*dx, 0];
                    source.uz = focus(kgrid, input_signal./(C0*RHO0), source.u_mask, focus_position, C0);
                case 3
                    focus_position = [0, 0, 0];
                    source.uz = focus(kgrid, input_signal./(C0*RHO0), source.u_mask, focus_position, C0);
            end
        else
            source.uz = input_signal./(C0*RHO0);
        end
        
    case 5
        
        % -----------------------------
        % OPTION 5: Velocity XYZ Source
        % ----------------------------- 
        switch TEST_DIM
            case 2
                source.u_mask = makeDisc(NX, NY, round(NX/2), round(NY/2), 3);
                source.ux = input_signal./(C0*RHO0);
                source.uy = input_signal./(C0*RHO0);
            case 3
                source.u_mask = makeBall(NX, NY, NZ, round(NX/2), round(NY/2), round(NZ/2), 3);
                source.ux = input_signal./(C0*RHO0);
                source.uy = input_signal./(C0*RHO0);
                source.uz = input_signal./(C0*RHO0);
            otherwise
                error('x/y/z velocity source not supported in 1D');
        end
        
    case 6
        
        % -----------------------------
        % OPTION 6: Transducer Source
        % -----------------------------
        
        % scale the source magnitude by the source_strength divided by the
        % impedance (the source is assigned to the particle velocity)
        input_signal = input_signal./(C0*RHO0);

        % physical properties of the transducer
        transducer.number_elements = 32;    % total number of transducer elements
        transducer.element_width = 1;       % width of each element [grid points]
        transducer.element_length = 12;     % length of each element [grid points]
        transducer.element_spacing = 0;     % spacing (kerf  width) between the elements [grid points]
        transducer.radius = inf;            % radius of curvature of the transducer [m]

        % calculate the width of the transducer in grid points
        transducer_width = transducer.number_elements*transducer.element_width ...
            + (transducer.number_elements - 1)*transducer.element_spacing;

        % use this to position the transducer in the middle of the computational grid
        transducer.position = round([PML_X_SIZE + 1, round(NY/2) - transducer_width/2, round(NZ/2) - transducer.element_length/2]);

        % properties used to derive the beamforming delays
        transducer.sound_speed = 1540;              % sound speed [m/s]
        transducer.focus_distance = 20e-3;          % focus distance [m]
        transducer.elevation_focus_distance = 19e-3;% focus distance in the elevation plane [m]
        transducer.steering_angle = 0;              % steering angle [degrees]

        % apodization
        transducer.transmit_apodization = 'Rectangular';    
        transducer.receive_apodization = 'Rectangular';

        % define the transducer elements that are currently active
        transducer.active_elements = ones(transducer.number_elements, 1);

        % append input signal used to drive the transducer
        transducer.input_signal = input_signal;

        % create the transducer using the defined settings
        transducer = makeTransducer(kgrid, transducer);
        
        % pass the transducer to the source
        source = transducer;
        
end

% dirichlet or additive source condition
if (SOURCE_TYPE == 1) && SOURCE_DIRICHLET
    source.p_mode = 'dirichlet';
end

if (SOURCE_TYPE > 1) && (SOURCE_TYPE < 6) && SOURCE_DIRICHLET
    source.u_mode = 'dirichlet';
end

% =========================================================================
% DEFINE SENSOR MASK
% =========================================================================

if BINARY_SENSOR_MASK

    % define a sensor mask through the central plane
    switch TEST_DIM
        case 1
            sensor.mask = zeros(NX, 1);
            sensor.mask(round(NX/2):round(NX/2)+1) = 1;
        case 2
            sensor.mask = zeros(NX, NY);
            sensor.mask(:, round(NZ/2)) = 1;
            
            if SENSOR_DIRECTIVITY
                sensor.directivity_angle = zeros(NX, NY);
                sensor.directivity_angle(sensor.mask == 1) = 0;
                sensor.directivity_size = 3*kgrid.dx;
            end
        case 3
            sensor.mask = zeros(NX, NY, NZ);
            sensor.mask(:, :, round(NZ/2)) = 1;
    end

else
    
    % define a series of Cartesian points to collect the data
    switch TEST_DIM
        case 1
            x = (-22:2:22)*dx;          % [m]
            sensor.mask = x;
        case 2
            x = (-22:2:22)*dx;          % [m]
            y = (-11:1:11)*dy;          % [m]
            sensor.mask = [x; y];
        case 3
            x = (-22:2:22)*dx;          % [m]
            y = (-11:1:11)*dy;          % [m]
            z = 10*dz*ones(size(x));    % [m]
            sensor.mask = [x; y; z];
    end
    
    % add interpolation option to input arguments
    input_args = [input_args, {'CartInterp', CART_INTERP}];
    
end

% sensor record option
if ~isempty(SENSOR_RECORD_FIELDS)
    sensor.record = SENSOR_RECORD_FIELDS;
end

% sensor frequency response
if ~isempty(FREQUENCY_RESPONSE)
    sensor.frequency_response = FREQUENCY_RESPONSE;
end

% sensor record start time
if RECORD_START_INDEX > 1
    sensor.record_start_index = RECORD_START_INDEX;
end

% =========================================================================
% SET OPTIONAL INPUT ARGUMENTS
% =========================================================================

% stream to disk option
if ~isempty(STREAM_TO_DISK)
    input_args = [input_args, {'StreamToDisk', STREAM_TO_DISK}];
end

% save logfile
if ~isempty(CREATE_LOG)
    input_args = [input_args, {'CreateLog', CREATE_LOG}];
end

% plot pml
if ~isempty(PLOT_PML)
    input_args = [input_args, {'PlotPML', PLOT_PML}];
end

% logscale plot
if ~isempty(PLOT_LOG_SCALE)
    input_args = [input_args, {'LogScale', PLOT_LOG_SCALE}];
end

% display mask
if ~isempty(DISPLAY_MASK)
    if strcmp(DISPLAY_MASK, 'ball')
        switch TEST_DIM
            case 1
                DISPLAY_MASK = zeros(NX, 1);
                DISPLAY_MASK(round(NX/2):round(NX/2)+1) = 1;
            case 2
                DISPLAY_MASK = makeDisc(NX, NY, round(NX/2), round(NY/2), 10);
            case 3
                DISPLAY_MASK = makeBall(NX, NY, NZ, round(NX/2), round(NY/2), round(NZ/2), 10);
        end
    end
    input_args = [input_args, {'DisplayMask', DISPLAY_MASK}];
end

% plot frequency
if ~isempty(PLOT_FREQ)
    input_args = [input_args, {'PlotFreq', PLOT_FREQ}];
end

% plot scale
if isempty(PLOT_SCALE)
    if TEST_DIM == 1
        PLOT_SCALE = [-source_strength*1.2, source_strength*1.2];
    else
        PLOT_SCALE = [-source_strength/2, source_strength/2];
    end
end
    
% plot on or off
if FORCE_PLOT_OFF
    input_args = [input_args, {'PlotSim', false}];
elseif ~isempty(PLOT_SIM)
    input_args = [input_args, {'PlotSim', PLOT_SIM}];
end
    
% plot layout
if ~isempty(PLOT_LAYOUT)
    input_args = [input_args, {'PlotLayout', PLOT_LAYOUT}];
end

% save movie
if ~isempty(SAVE_MOVIE)
    movie_name = [OUTPUT_FOLDER getDateString '_kWave_Test_' num2str(comp_index)];
    input_args = [input_args, {'RecordMovie', SAVE_MOVIE, 'MovieName', movie_name}];
end

% smoothing
if SMOOTH
    input_args = [input_args, {'Smooth', SMOOTH}];
end

% mesh plot
if ~isempty(MESH_PLOT)
    input_args = [input_args, {'MeshPlot', true}];
end

if ~isempty(MOVIE_TYPE_IMAGE)
    input_args = [input_args, {'MovieType', 'image'}];
end

% PML size
switch TEST_DIM
    case 1
        input_args = [input_args, {'PMLSize', PML_X_SIZE}];
    case 2
        input_args = [input_args, {'PMLSize', [PML_X_SIZE, PML_Y_SIZE]}];
    case 3
        input_args = [input_args, {'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE]}];
end

% set additional input settings
input_args = [input_args, {'PMLInside', PML_INSIDE,...
    'DataCast', DATA_CAST, 'DataRecast', DATA_RECAST, 'PlotScale', PLOT_SCALE}];

% =========================================================================
% RUN THE MATLAB SIMULATION
% =========================================================================

% run the simulation using k-Wave
switch TEST_DIM
    case 1
        sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, input_args{:});
    case 2
        sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    case 3
        sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
end

% run the time reversal simulation 
if RUN_TIME_REVERSAL
   
    % clear the source distribution
    if isfield(source, 'p0');
        source.p0 = 0;
    end
    
    % assign the time reversal data
    if isstruct(sensor_data)
        sensor.time_reversal_boundary_data = sensor_data.p;
    else
        sensor.time_reversal_boundary_data = sensor_data;
    end

    % attenuation compensation options
    if RUN_TIME_REVERSAL == 2
        
        % define the cutoff frequency for the filter
        f_cutoff = 1e6;

        % create the filter to regularise the absorption parameters
        medium.alpha_filter = getAlphaFilter(kgrid, medium, f_cutoff);

        % reverse the sign of the absorption proportionality coefficient
        medium.alpha_sign = [-1, 1];        % [absorption, dispersion];
        
    end
    
    % run the reconstruction
    switch TEST_DIM
        case 1
            p0_recon = kspaceFirstOrder1D(kgrid, medium, source, sensor, input_args{:});
        case 2
            p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
        case 3
            p0_recon = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
    end
    
end

% =========================================================================
% RUN THE C++ SIMULATION
% =========================================================================

% compare the output of the MATLAB code with the C++ if required
if RUN_CPP_CODE 

    % run the simulation again using Jiri's C++ version
    if TEST_MPI
        sensor_data_cpp = kspaceFirstOrder3DMPI(kgrid, medium, source, sensor, input_args{:});
    else
        sensor_data_cpp = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});
    end
    
    if RUN_TIME_REVERSAL

        % clear the source distribution
        if isfield(source, 'p0');
            source.p0 = 0;
        end

        % assign the time reversal data
        if isstruct(sensor_data)
            sensor.time_reversal_boundary_data = sensor_data.p;
        else
            sensor.time_reversal_boundary_data = sensor_data;
        end

        % run the reconstruction
        if TEST_MPI
            p0_recon_cpp = kspaceFirstOrder3DMPI(kgrid, medium, source, sensor, input_args{:});
        else
            p0_recon_cpp = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});
        end

    end

    % set plot axis
    x_axis = (kgrid.x_vec - min(kgrid.x_vec(:)))*1e3;
    y_axis = kgrid.y_vec*1e3;
    
    % number of time steps
    if isfield(sensor, 'record_start_index')
        NT = kgrid.Nt - sensor.record_start_index + 1;
    else
        NT = kgrid.Nt;
    end
    
    % check to see if data should be recast to the CPU
    if ismember(DATA_CAST, {'gsingle', 'GPUsingle', 'gpuArray-single'})
        recast = @(x) single(x);
    else
        recast = @(x) x;
    end
    
    % set command line output
    disp(' ');
    disp(' ');
    disp('  C++ ACCURACY COMPARED TO MATLAB:');
    disp('  --------------------------------');
    
    % compare outputs for p0_recon
    if RUN_TIME_REVERSAL
        
        % reshape the output data
        mat = reshape(recast(p0_recon), [NX, NY, NZ]);               
        cpp = reshapeCPP(recast(p0_recon_cpp), [NX, NY, NZ], false, PML_INSIDE, PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE);
        
        % display error norms
        [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'p0_recon', comp_index, number_cpp_errors, location_cpp_errors);
        
        % plot
        if PLOT_CPP_ERRORS
            mat = squeeze(mat(round(end/2), :, :));
            cpp = squeeze(cpp(round(end/2), :, :));
            plot_title = 'p0 recon';
            plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);
        end
        
    else

        % compare the outputs for p
        if (isfield(sensor, 'record') && ismember('p', sensor.record)) || ~isfield(sensor, 'record')

            % reshape the output data
            if isfield(sensor, 'record')
                mat = reshape(recast(sensor_data.p), [NX, NY, NT]);
                cpp = reshape(recast(sensor_data_cpp.p), [NX, NY, NT]);
            else
                mat = reshape(recast(sensor_data), [NX, NY, NT]);
                cpp = reshape(recast(sensor_data_cpp), [NX, NY, NT]);
            end

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data.p', comp_index, number_cpp_errors, location_cpp_errors);

            % calculate the distribution of fundamental and harmonic
            if NONLINEAR && PLOT_CPP_ERRORS

                % compute the frequency axis
                freq = (0:ceil((NT + 1)/2) - 1) ./ (kgrid.dt*NT);

                % compute the index at which the source frequency and its harmonics occur
                [~, f0_index] = findClosest(freq, tone_burst_freq);
                [~, f1_index] = findClosest(freq, tone_burst_freq*2);

                % preallocate the beam pattern variables
                beam_pattern_f0 = zeros(NX, NY);
                beam_pattern_f1 = zeros(NX, NY);
                beam_pattern_total = zeros(NX, NY);
                beam_pattern_f0_cpp = zeros(NX, NY);
                beam_pattern_f1_cpp = zeros(NX, NY);
                beam_pattern_total_cpp = zeros(NX, NY);

                % compute the amplitude spectrum of the time series recorded at each sensor
                % point, and then extract the corresponding amplitudes at the fundamental
                % frequency and second harmonic.
                for x_index = 1:NX
                    for y_index = 1:NY

                        % compute the amplitude spectrum
                        amp_spect = spect(squeeze(mat(x_index, y_index, :)), 1/kgrid.dt); %, 'Window', 'Hanning');
                        amp_spect_cpp = spect(squeeze(cpp(x_index, y_index, :)), 1/kgrid.dt);

                        % extract the amplitude at the source frequency and store
                        beam_pattern_f0(x_index, y_index) = amp_spect(f0_index);
                        beam_pattern_f0_cpp(x_index, y_index) = amp_spect_cpp(f0_index);

                        % extract the amplitude at the source frequency and store
                        beam_pattern_f1(x_index, y_index) = amp_spect(f1_index); 
                        beam_pattern_f1_cpp(x_index, y_index) = amp_spect_cpp(f1_index);

                        % extract the integral of the total amplitude spectrum
                        beam_pattern_total(x_index, y_index) = sum(amp_spect(:));
                        beam_pattern_total_cpp(x_index, y_index) = sum(amp_spect_cpp(:));

                    end
                end

                % produce plots
                mat = beam_pattern_total;
                cpp = beam_pattern_total_cpp;
                plot_title = 'Total Pressure';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);

                mat = beam_pattern_f0;
                cpp = beam_pattern_f0_cpp;          
                plot_title = 'Fundamental';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);

                mat = beam_pattern_f1;
                cpp = beam_pattern_f1_cpp; 
                plot_title = '2nd Harmonic';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);

            elseif PLOT_CPP_ERRORS

                % compute the maximum value to plot
                mat = max(mat, [], 3);
                cpp = max(cpp, [], 3);
                plot_title = 'Total Pressure';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);

            end
        end

        % compare the outputs for p_max
        if isfield(sensor, 'record') && ismember('p_max', sensor.record)

            % reshape the output data
            mat = reshape(recast(sensor_data.p_max), [NX, NY]);
            cpp = reshape(recast(sensor_data_cpp.p_max), [NX, NY]);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data.p_max', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if PLOT_CPP_ERRORS
                plot_title = 'p max';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);
            end

        end

        % compare the outputs for p_rms
        if isfield(sensor, 'record') && ismember('p_rms', sensor.record)

            % reshape the output data
            mat = reshape(recast(sensor_data.p_rms), [NX, NY]);
            cpp = reshape(recast(sensor_data_cpp.p_rms), [NX, NY]);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data.p_rms', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if PLOT_CPP_ERRORS
                plot_title = 'p rms';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);
            end

        end    

        % compare the outputs for p_final
        if isfield(sensor, 'record') && ismember('p_final', sensor.record)

            % reshape the output data
            mat = reshape(recast(sensor_data.p_final), [NX, NY, NZ]);
            cpp = reshapeCPP(recast(sensor_data_cpp.p_final), [NX, NY, NZ], false, PML_INSIDE, PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data.p_final', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if PLOT_CPP_ERRORS
                mat = squeeze(mat(round(end/2), :, :));
                cpp = squeeze(cpp(round(end/2), :, :));
                plot_title = 'p final';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);
            end

        end      

        % compare outputs for u
        if isfield(sensor, 'record') && ismember('u', sensor.record)

            % reshape the output data
            mat = reshape(recast(sensor_data.ux), [NX, NY, NT]);
            cpp = reshape(recast(sensor_data_cpp.ux), [NX, NY, NT]);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data.ux', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if PLOT_CPP_ERRORS
                mat = max(mat, [], 3);
                cpp = max(cpp, [], 3);
                plot_title = 'ux';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);
            end

            % reshape the output data
            mat = reshape(recast(sensor_data.uy), [NX, NY, NT]);
            cpp = reshape(recast(sensor_data_cpp.uy), [NX, NY, NT]);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data.uy', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if PLOT_CPP_ERRORS
                mat = max(mat, [], 3);
                cpp = max(cpp, [], 3);
                plot_title = 'uy';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);
            end

            % reshape the output data
            mat = reshape(recast(sensor_data.uz), [NX, NY, NT]);
            cpp = reshape(recast(sensor_data_cpp.uz), [NX, NY, NT]);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data.uz', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if PLOT_CPP_ERRORS
                mat = max(mat, [], 3);
                cpp = max(cpp, [], 3);
                plot_title = 'uz';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);
            end

        end

        % compare the outputs for u_max
        if isfield(sensor, 'record') && ismember('u_max', sensor.record)

            % reshape the output data
            mat = reshape(recast(sensor_data.ux_max), [NX, NY]);
            cpp = reshape(recast(sensor_data_cpp.ux_max), [NX, NY]);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data.ux_max', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if PLOT_CPP_ERRORS
                plot_title = 'ux max';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);
            end

            % reshape the output data
            mat = reshape(recast(sensor_data.uy_max), [NX, NY]);
            cpp = reshape(recast(sensor_data_cpp.uy_max), [NX, NY]);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data.uy_max', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if PLOT_CPP_ERRORS
                plot_title = 'uy max';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);
            end

            % reshape the output data
            mat = reshape(recast(sensor_data.uz_max), [NX, NY]);
            cpp = reshape(recast(sensor_data_cpp.uz_max), [NX, NY]);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data.uz_max', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if PLOT_CPP_ERRORS
                plot_title = 'uz max';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);
            end        

        end

        % compare the outputs for u_rms
        if isfield(sensor, 'record') && ismember('u_rms', sensor.record)

            % reshape the output data
            mat = reshape(recast(sensor_data.ux_rms), [NX, NY]);
            cpp = reshape(recast(sensor_data_cpp.ux_rms), [NX, NY]);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data.ux_rms', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if PLOT_CPP_ERRORS
                plot_title = 'ux rms';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);
            end

            % reshape the output data
            mat = reshape(recast(sensor_data.uy_rms), [NX, NY]);
            cpp = reshape(recast(sensor_data_cpp.uy_rms), [NX, NY]);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data.uy_rms', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if PLOT_CPP_ERRORS
                plot_title = 'uy rms';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);
            end

            % reshape the output data
            mat = reshape(recast(sensor_data.uz_rms), [NX, NY]);
            cpp = reshape(recast(sensor_data_cpp.uz_rms), [NX, NY]);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data.uz_rms', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if PLOT_CPP_ERRORS
                plot_title = 'uz rms';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);
            end        

        end    

        % compare the outputs for u_final
        if isfield(sensor, 'record') && ismember('u_final', sensor.record)

            % reshape the output data
            mat = reshape(recast(sensor_data.ux_final), [NX, NY, NZ]);
            cpp = reshapeCPP(recast(sensor_data_cpp.ux_final), [NX, NY, NZ], false, PML_INSIDE, PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data.ux_final', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if PLOT_CPP_ERRORS
                mat = squeeze(mat(round(end/2), :, :));
                cpp = squeeze(cpp(round(end/2), :, :));
                plot_title = 'ux final';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);
            end

            mat = reshape(recast(sensor_data.uy_final), [NX, NY, NZ]);
            cpp = reshapeCPP(recast(sensor_data_cpp.uy_final), [NX, NY, NZ], false, PML_INSIDE, PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data.uy_final', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if PLOT_CPP_ERRORS
                mat = squeeze(mat(round(end/2), :, :));
                cpp = squeeze(cpp(round(end/2), :, :));
                plot_title = 'uy final';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);
            end

            mat = reshape(recast(sensor_data.uz_final), [NX, NY, NZ]);
            cpp = reshapeCPP(recast(sensor_data_cpp.uz_final), [NX, NY, NZ], false, PML_INSIDE, PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data.uz_final', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if PLOT_CPP_ERRORS
                mat = squeeze(mat(round(end/2), :, :));
                cpp = squeeze(cpp(round(end/2), :, :));
                plot_title = 'uz final';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);
            end        

        end          

        % compare outputs for I_avg
        if isfield(sensor, 'record') && ismember('I_avg', sensor.record)

            % reshape the output data
            mat = reshape(recast(sensor_data.Ix_avg), [NX, NY]);
            cpp = reshape(recast(sensor_data_cpp.Ix_avg), [NX, NY]);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data.Ix_avg', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if PLOT_CPP_ERRORS
                plot_title = 'Ix avg';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);
            end

            % reshape the output data
            mat = reshape(recast(sensor_data.Iy_avg), [NX, NY]);
            cpp = reshape(recast(sensor_data_cpp.Iy_avg), [NX, NY]);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data.Iy_avg', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if PLOT_CPP_ERRORS
                plot_title = 'Iy avg';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);
            end

            % reshape the output data
            mat = reshape(recast(sensor_data.Iz_avg), [NX, NY]);
            cpp = reshape(recast(sensor_data_cpp.Iz_avg), [NX, NY]);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data.Iz_avg', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if PLOT_CPP_ERRORS
                plot_title = 'Iz avg';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, SAVE_CPP_COMPARISON_PLOTS_TO_DISK, [OUTPUT_FOLDER IMAGE_FOLDERNAME]);
            end

        end   
    end 
end

% end for subfunction
end

% end for parent function
end


% subfunction to reshape CPP output data depending on position of PML
function cpp = reshapeCPP(cpp, data_sz, z_is_nt, pml_inside, pml_x_size, pml_y_size, pml_z_size)

if ~pml_inside
    if z_is_nt
        cpp = reshape(cpp, [data_sz(1) + 2*pml_x_size, data_sz(2) + 2*pml_y_size, data_sz(3)]);
        cpp = cpp(1 + pml_x_size:end - pml_x_size, 1 + pml_y_size:end - pml_y_size, :);
    else
        cpp = reshape(cpp, [data_sz(1) + 2*pml_x_size, data_sz(2) + 2*pml_y_size, data_sz(3) + 2*pml_z_size]);
        cpp = cpp(1 + pml_x_size:end - pml_x_size, 1 + pml_y_size:end - pml_y_size, 1 + pml_z_size:end - pml_z_size);
    end
else
    cpp = reshape(cpp, data_sz);
end

end

% subfunction to calculate and display error norms
function [num_errors, loc_errors] = errorNorms(mat, cpp, error_title, test_index, num_errors, loc_errors)

% calculate and display the error norms
L2 = sqrt( sum(abs(mat(:).^2 - cpp(:).^2)) / sum(mat(:).^2) );
LINF = max(abs(mat(:) - cpp(:)));
disp(' ');
disp(['  Error in ' error_title]);
disp(['    L2 = ' num2str(L2)]);
disp(['    LINF = ' num2str(LINF) ' (' num2str(LINF/max(abs(mat(:)))) ')' ]);

if (LINF/max(abs(mat(:)))) > 1e-5
    num_errors = num_errors + 1;
    if isempty(loc_errors)
        loc_errors = test_index;
    elseif loc_errors(end) ~= test_index
        loc_errors = [loc_errors, test_index];
    end
end

end

% subfunction to produce plots
function plotErrors(x_axis, y_axis, matlab_data, cpp_data, test_index, plot_title, save_plots, image_folder)
    
    % produce plot of errors
    h = figure;
    subplot(1, 4, 1);
    imagesc(y_axis, x_axis, matlab_data/1e6);
    set(gca, 'FontSize', 10);
    xlabel('y-position [mm]');
    ylabel('x-position [mm]');
    title(['k-Wave (' plot_title ')']);
    colormap(jet(256));
    c = colorbar;
    ylabel(c, 'Pressure [MPa]');
    axis image;

    subplot(1, 4, 2);
    imagesc(y_axis, x_axis, cpp_data/1e6);
    set(gca, 'FontSize', 10);
    xlabel('y-position [mm]');
    ylabel('x-position [mm]');
    title(['C++ (' plot_title ')']);
    colormap(jet(256));
    c = colorbar;
    ylabel(c, 'Pressure [MPa]');
    axis image;

    subplot(1, 4, 3);
    imagesc(y_axis, x_axis, 100*abs(matlab_data - cpp_data)./matlab_data);
    set(gca, 'FontSize', 10);
    xlabel('y-position [mm]');
    ylabel('x-position [mm]');
    title('Local Error');
    colormap(jet(256));
    c = colorbar;
    ylabel(c, '[%]');
    axis image;

    subplot(1, 4, 4);
    imagesc(y_axis, x_axis, 100*abs(matlab_data - cpp_data)./max(abs(matlab_data(:))));
    set(gca, 'FontSize', 10);
    xlabel('y-position [mm]');
    ylabel('x-position [mm]');
    title('Global Error');
    colormap(jet(256));
    c = colorbar;
    ylabel(c, '[%]');
    axis image;
    
    % scale figure
    scaleFig(2, 1);
    
    % save plots to disk and close
    if save_plots
        print(h, '-dpng','-r300', [image_folder 'test' num2str(test_index) '-' plot_title  '.png']);
        close all hidden; 
    end
    
end