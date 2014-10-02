% Script to call kWavetester to test the functionality of kspaceFirstOrder1D,
% kspaceFirstOrder2D, and kspaceFirstOrder3D. 
%
% author: Bradley Treeby
% date: 26th August 2014
% last update: 26th August 2014

clear all;

% =========================================================================
% GENERAL OPTIONS
% =========================================================================

% 1: kspaceFirstOrder1D
% 2: kspaceFirstOrder2D
% 3: kspaceFirstOrder3D
options.test_dim = 3;

% 1: Even grid size, PML inside
% 2: Even grid size, PML outside
% 3: Odd grid size, PML inside
% 4: Odd grid size, PML outside
% 5: Custom grid and PML size
options.test_type = 5;

% custom grid size (used only if options.test_type = 5)
options.Nx                                      = 64;
options.Ny                                      = 64;
options.Nz                                      = 64;
options.pml_x_size                              = 6;
options.pml_y_size                              = 6;
options.pml_z_size                              = 6;
options.pml_inside                              = true;

% =========================================================================
% MATLAB OPTIONS
% =========================================================================

% data format for MATLAB simulations 'off', 'single', or 'gpuArray-single'
options.data_cast                               = 'single';  

% option to use non-uniform grid (NOT CURRENTLY USED)
options.use_nonuniform_grid                     = false;

% force MATLAB k-Wave plotting to be off
options.force_plot_off                          = false;

% =========================================================================
% C++ OPTIONS
% =========================================================================

% comparison with C++ code
options.run_cpp_comparison_tests                = true;
options.plot_cpp_errors                         = false;
options.save_cpp_comparison_plots_to_disk       = false;
options.image_foldername                        = 'CPP_COMPARISON_IMAGES/';

% switch to test CUDA version instead of OpenMP
options.use_gpu_code                            = false;

% path to cpp binary
options.cpp_binary_path                         = getkWavePath('binaries');

% name of cpp binary to test
if isunix
    options.cpp_binary_name                     = 'kspaceFirstOrder3D-OMP';
else
    options.cpp_binary_name                     = 'kspaceFirstOrder3D-OMP.exe';
end

% option to just save the HDF5 input file and not run the test 
% (only used if options.custom_test and options.run_cpp_comparison_tests
% are true)
options.cpp_save_to_disk_only                   = false;
options.cpp_save_to_disk_filename               = 'test_input.h5';

% =========================================================================
% SET LOCATION OF OUTPUT FOLDER
% =========================================================================

% where to put all the output data
options.output_folder                       = '/';

% =========================================================================
% TEST OPTIONS
% =========================================================================

% tests to run -> these automatically run different sets of input and
% output combinations, alternatively, set options.custom_test to true
options.run_source_tests                        = true;     % outputs are compared with C++ code
options.run_bin_sensor_tests                    = true;     % outputs are compared with C++ code
options.run_cuboid_corner_tests                 = true;     % outputs are compared with C++ code
options.run_cart_sensor_lin_interp_tests        = true;     % outputs are NOT compared with C++ code
options.run_cart_sensor_nearest_interp_tests    = true;     % outputs are NOT compared with C++ code
options.run_display_tests                       = true;     % outputs are NOT compared with C++ code
options.run_time_reversal_tests                 = true;     % outputs are compared with C++ code

% set this to true to run a single test as defined below
options.custom_test                             = true;

% what to do with the output data
options.save_test_log                           = false;
options.split_log_after_n_tests                 = 500;
options.continue_past_errors                    = false;

% set start index to skip tests on restart (set to 1 to run all tests)
options.test_case_start_index                   = 1;

% overwrite test case if required (only used if options.custom_test = true)
% -----------------------------------------------------------------------------------------------------
options.custom_test_case = [...
% -----------------------------------------------------------------------------------------------------
%   LIN ABS HET SMH ALP SRC MNY DIR BIN P-U STR RCT FRQ RSI LOG PML LSC DIS PFQ PSC PLT LAY MOV TRV CPP   
%   01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25 
% ----------------------------------------------------------------------------------------------------- 
    0   0   0   0   0   0   0   0   4   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1; ...  
% -----------------------------------------------------------------------------------------------------
];

% -----------------------------------------------------------------------
% PARAMETER 01 (LIN)
%   0: Linear
%   1: Nonlinear
% PARAMETER 02 (ABS)
%   0: Lossless
%   1: Absorbing
% PARAMETER 03 (HET)
%   0: Homogeneous
%   1: Heterogeneous c0 and rho0 only
%   2: Heterogeneous
% PARAMETER 04 (SMH)
%   0: 'Smooth' not set (defaults to [true, false, false])
%   1: 'Smooth', true
% PARAMETER 05 (ALP) 
%   0: medium.alpha_mode not set
%   1: medium.alpha_mode = 'no_absorption'
%   2: medium.alpha_mode = 'no_dispersion'
% -----------------------------------------------------------------------
% PARAMETER 06 (SRC)
%   0: initial pressure
%   1: pressure source
%   2: velocity-x source
%   3: velocity-y source (2D and 3D only)
%   4: velocity-z source (3D only)
%   5: velocity-x/y/z source
%   6: transducer source (3D only)
% PARAMETER 07 (MNY)
%   0: single source
%   1: source many
% PARAMETER 08 (DIR)
%   0: additive source condition
%   1: dirichlet source condition
% -----------------------------------------------------------------------
% PARAMETER 09 (BIN)
%   0: binary sensor mask
%   1: Cartesian sensor mask, linear interpolation
%   2: Cartesian sensor mask, nearest neighbour interpolation
%   3: binary sensor mask with sensor directivity (2D only)
%   4: cuboid corners sensor mask
% PARAMETER 10 (P-U)
%   0: record only pressure (no sensor.record input)
%   1: record pressure and velocity
%   2: record max, min and rms pressure and velocity
%   3: record everything BUT u_non_staggered or intensity
%   4: record everything
% PARAMETER 11 (STR) (3D only)
%   0: 'StreamToDisk', false
%   1: 'StreamToDisk', true
%   2: 'StreamToDisk', 50
% PARAMETER 12 (RCT)
%   0: 'DataRecast', false
%   1: 'DataRecast', true
% PARAMETER 13 (FRQ)
%   0: sensor.frequency_response not set
%   1: sensor.frequency_response = [0.5e6, 50]
% PARAMETER 14 (RSI)
%   0: sensor.record_start_index not set
%   1: sensor.record_start_index = 200;
% -----------------------------------------------------------------------
% PARAMETER 15 (LOG)
%   0: 'CreateLog', false (default)
%   1: 'CreateLog', true
% PARAMETER 16 (PML)
%   0: 'PlotPML', true (default)
%   1: 'PlotPML', false (default)
% PARAMETER 17 (LSC)
%   0: 'LogScale', false (default)
%   1: 'LogScale', true
% PARAMETER 18 (DIS)
%   0: 'DisplayMask' not set (defaults to sensor.mask)
%   1: 'DisplayMask', 'off'
%   2: 'DisplayMask', makeBall
% PARAMETER 19 (PFQ)
%   0: 'PlotFreq' not set (defaults to 10) 
%   1: 'PlotFreq', 30
% PARAMETER 20 (PSC)
%   0: 'PlotScale' set to -[p0, p0]
%   1: 'PlotScale' set to 'auto'
% PARAMETER 21 (PLT)
%   0: 'PlotSim' not set (defaults to true) 
%   1: 'PlotSim', false
%   2: 'MeshPlot', true (2D only)
% PARAMETER 22 (LAY)
%   0: 'PlotLayout' not set (defaults to false)
%   1: 'PlotLayout', true
% PARAMETER 23 (MOV)
%   0: 'RecordMovie' not set (defaults to false)
%   1: 'RecordMovie', true
%   2: 'RecordMovie', true / 'MovieType', 'image' (2D only)
% -----------------------------------------------------------------------
% PARAMETER 24 (TRV)
%   0: run only forward simulation
%   1: run time reversal simulation
%   2: run time reversal including attenuation compensation
% -----------------------------------------------------------------------
% PARAMETER 25 (CPP) (3D only)
%   0: run MATLAB version only
%   1: run C++ version and compare outputs
% -----------------------------------------------------------------------

% =========================================================================
% RUN TEST
% =========================================================================

kwaveTester(options);