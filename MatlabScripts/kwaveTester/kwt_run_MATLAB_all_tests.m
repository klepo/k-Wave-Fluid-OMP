% Script to call kWavetester to test the functionality of kspaceFirstOrder1D,
% kspaceFirstOrder2D, and kspaceFirstOrder3D. 
%
% author: Bradley Treeby
% date: 26th August 2014
% last update: 26th August 2014

clear all;

% =========================================================================
% MATLAB OPTIONS
% =========================================================================

% option to use non-uniform grid (NOT CURRENTLY USED)
options.use_nonuniform_grid                     = false;

% force MATLAB k-Wave plotting to be off
options.force_plot_off                          = false;

% run GPU test
RUN_GPU_TESTS                                   = false;

% =========================================================================
% C++ OPTIONS
% =========================================================================

% comparison with C++ code
options.run_cpp_comparison_tests                = false;
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

% =========================================================================
% SET LOCATION OF OUTPUT FOLDER
% =========================================================================

% computer name
[~, comp] = unix('hostname', '-echo');
comp = comp(1:end-1);

% where to put all the output data
if ismac
    options.output_folder                       = '~/Documents/MATLAB/kWaveTests/';
elseif isunix
    if strcmp(comp, 'rayleigh')
        options.output_folder                   = '~/Documents/kwave_tests/';
    elseif strcmp(comp, 'nyborg')
        options.output_folder                   = '/allusers/k-Wave-Tests/';
    elseif strcmp(comp, 'helmholtz')
        options.output_folder                   = '~/kWaveTests/';
    else
        options.output_folder                   = '~/kWaveTests/';
    end
else
    options.output_folder                       = 'C:\Users\Bradley Treeby\Documents\MATLAB\k-Wave Tests\Results\';
end

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
options.custom_test                             = false;

% what to do with the output data
options.save_test_log                           = true;
options.split_log_after_n_tests                 = 500;
options.continue_past_errors                    = false;

% set start index to skip tests on restart (set to 1 to run all tests)
options.test_case_start_index                   = 1;

% =========================================================================
% RUN 1D TESTS
% =========================================================================

% 1: kspaceFirstOrder1D
options.test_dim = 1;

% ----

% data format for MATLAB simulations
options.data_cast = 'off';  

% 1: Even grid size, PML inside
options.test_type = 1;
kwaveTester(options);
close all hidden;

% 2: Even grid size, PML outside
options.test_type = 2;
kwaveTester(options);
close all hidden;

% 3: Odd grid size, PML inside
options.test_type = 3;
kwaveTester(options);
close all hidden;

% 4: Odd grid size, PML outside
options.test_type = 4;
kwaveTester(options);
close all hidden;

% ----

% data format for MATLAB simulations
options.data_cast = 'single';  

% 1: Even grid size, PML inside
options.test_type = 1;
kwaveTester(options);
close all hidden;

% 2: Even grid size, PML outside
options.test_type = 2;
kwaveTester(options);
close all hidden;

% 3: Odd grid size, PML inside
options.test_type = 3;
kwaveTester(options);
close all hidden;

% 4: Odd grid size, PML outside
options.test_type = 4;
kwaveTester(options);
close all hidden;

% ----

if RUN_GPU_TESTS

    % data format for MATLAB simulations
    options.data_cast = 'gpuArray-single';  

    % 1: Even grid size, PML inside
    options.test_type = 1;
    kwaveTester(options);
    close all hidden;

    % 2: Even grid size, PML outside
    options.test_type = 2;
    kwaveTester(options);
    close all hidden;

    % 3: Odd grid size, PML inside
    options.test_type = 3;
    kwaveTester(options);
    close all hidden;

    % 4: Odd grid size, PML outside
    options.test_type = 4;
    kwaveTester(options);
    close all hidden;

end

% =========================================================================
% RUN 2D TESTS
% =========================================================================

% 1: kspaceFirstOrder1D
options.test_dim = 2;

% ----

% data format for MATLAB simulations
options.data_cast = 'off';  

% 1: Even grid size, PML inside
options.test_type = 1;
kwaveTester(options);
close all hidden;

% 2: Even grid size, PML outside
options.test_type = 2;
kwaveTester(options);
close all hidden;

% 3: Odd grid size, PML inside
options.test_type = 3;
kwaveTester(options);
close all hidden;

% 4: Odd grid size, PML outside
options.test_type = 4;
kwaveTester(options);
close all hidden;

% ----

% data format for MATLAB simulations
options.data_cast = 'single';  

% 1: Even grid size, PML inside
options.test_type = 1;
kwaveTester(options);
close all hidden;

% 2: Even grid size, PML outside
options.test_type = 2;
kwaveTester(options);
close all hidden;

% 3: Odd grid size, PML inside
options.test_type = 3;
kwaveTester(options);
close all hidden;

% 4: Odd grid size, PML outside
options.test_type = 4;
kwaveTester(options);
close all hidden;

% ----

if RUN_GPU_TESTS

    % data format for MATLAB simulations
    options.data_cast = 'gpuArray-single';  

    % 1: Even grid size, PML inside
    options.test_type = 1;
    kwaveTester(options);
    close all hidden;

    % 2: Even grid size, PML outside
    options.test_type = 2;
    kwaveTester(options);
    close all hidden;

    % 3: Odd grid size, PML inside
    options.test_type = 3;
    kwaveTester(options);
    close all hidden;

    % 4: Odd grid size, PML outside
    options.test_type = 4;
    kwaveTester(options);
    close all hidden;

end

% =========================================================================
% RUN 3D TESTS
% =========================================================================

% 1: kspaceFirstOrder1D
options.test_dim = 3;

% ----

% data format for MATLAB simulations
options.data_cast = 'off';  

% 1: Even grid size, PML inside
options.test_type = 1;
kwaveTester(options);
close all hidden;

% 2: Even grid size, PML outside
options.test_type = 2;
kwaveTester(options);
close all hidden;

% 3: Odd grid size, PML inside
options.test_type = 3;
kwaveTester(options);
close all hidden;

% 4: Odd grid size, PML outside
options.test_type = 4;
kwaveTester(options);
close all hidden;

% ----

% data format for MATLAB simulations
options.data_cast = 'single';  

% 1: Even grid size, PML inside
options.test_type = 1;
kwaveTester(options);
close all hidden;

% 2: Even grid size, PML outside
options.test_type = 2;
kwaveTester(options);
close all hidden;

% 3: Odd grid size, PML inside
options.test_type = 3;
kwaveTester(options);
close all hidden;

% 4: Odd grid size, PML outside
options.test_type = 4;
kwaveTester(options);
close all hidden;

% ----

if RUN_GPU_TESTS

    % data format for MATLAB simulations
    options.data_cast = 'gpuArray-single';  

    % 1: Even grid size, PML inside
    options.test_type = 1;
    kwaveTester(options);
    close all hidden;

    % 2: Even grid size, PML outside
    options.test_type = 2;
    kwaveTester(options);
    close all hidden;

    % 3: Odd grid size, PML inside
    options.test_type = 3;
    kwaveTester(options);
    close all hidden;

    % 4: Odd grid size, PML outside
    options.test_type = 4;
    kwaveTester(options);
    close all hidden;

end