function sensor_data = kspaceFirstOrder3DC(varargin)
%KSPACEFIRSTORDER3DC   3D time-domain simulation of wave propagation using C++ code
%
% DESCRIPTION:
%       kspaceFirstOrder3DC provides an blind interface to the C++ version
%       of kspaceFirstOrder3D (called kspaceFirstOrder3D-OMP). The function
%       works by appending the optional input 'SaveToDisk' to the user
%       inputs and then calling kspaceFirstOrder3D to save the input files
%       to disk. The contents of sensor.record (if set) are parsed as input
%       flags, and the C++ code is run using the system command. The output
%       files are then automatically loaded from disk and returned in the
%       same fashion as kspaceFirstOrder3D. The input and output files are
%       saved to the temporary directory native to the operating system,
%       and are deleted after the function runs.
%
%       This function requires the C++ binary/executable of
%       kspaceFirstOrder3D-OMP to be downloaded from
%       http://www.k-wave.org/download.php and placed in the "binaries"
%       directory of the k-Wave toolbox. 
% 
%       Note, not all input options are currently supported, and all
%       display options are ignored (only command line outputs are given).
%       See the k-Wave user manual for more information.
%
%       This function is not recommended for large simulations, as the
%       input variables will reside twice in main memory (once in MATLAB,
%       and once in C++). For large simulations, the C++ code should be
%       called outside of MATLAB. See the k-Wave manual for more
%       information.
%
% USAGE:
%       see kspaceFirstOrder3D
%
%
% ABOUT:
%       author      - Bradley Treeby and Jiri Jaros
%       date        - 3rd February 2012
%       last update - 1st October 2012
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2012 Bradley Treeby and Ben Cox
%
% See also kspaceFirstOrder3D

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

% check the binaries exist and are in the correct place before doing
% anything else
if (isunix && ~exist([getkWavePath 'binaries/kspaceFirstOrder3D-OMP'], 'file')) || ...
        (~isunix && ~exist([getkWavePath 'binaries\kspaceFirstOrder3D-OMP.exe'], 'file')) 
    error('To use C++ code, the C++ binaries for your operating system must be downloaded from www.k-wave.org/download.php and placed in the binaries folder.');
end

% set the filename inputs to store data in the default temp directory
date_string = getDateString;
input_filename = [tempdir 'kwave_input_data' date_string '.h5'];
output_filename = [tempdir 'kwave_output_data' date_string '.h5'];   

% set empty options string
options_string = '';

% assign pseudonyms for input structures
source = varargin{3};
sensor = varargin{4};

% check if performing time reversal, and replace inputs to explicitly use a
% source with a dirichlet boundary condition
if isfield(sensor, 'time_reversal_boundary_data')
        
    % define a new source structure
    clear source;
    source.p_mask = sensor.mask;
    source.p = flipdim(sensor.time_reversal_boundary_data, 2);
    source.p_mode = 'dirichlet';
    
    % define a new sensor structure
    clear sensor;
    sensor.mask = ones(varargin{1}.Nx, varargin{1}.Ny, varargin{1}.Nz);
    sensor.record = {'p_final'};
    
    % set time reversal flag
    time_rev = true;

else
    
    % set time reversal flag
    time_rev = false;
    
end

% check if sensor.record is given
if isfield(sensor, 'record')
    
    % set the options string to record the required output fields
    if ismember('p', sensor.record)
        options_string = [options_string ' --p_raw'];
    end
    if ismember('p_max', sensor.record)
        options_string = [options_string ' --p_max'];
    end
    if ismember('p_rms', sensor.record)
        options_string = [options_string ' --p_rms'];
    end    
    if ismember('p_final', sensor.record)
        options_string = [options_string ' --p_final'];
    end    
    if ismember('u', sensor.record)
        options_string = [options_string ' --u_raw'];
    end
    if ismember('u_max', sensor.record)
        options_string = [options_string ' --u_max'];
    end
    if ismember('u_rms', sensor.record)
        options_string = [options_string ' --u_rms'];
    end
    if ismember('u_final', sensor.record)
        options_string = [options_string ' --u_final'];
    end      
    if ismember('I_avg', sensor.record)
        options_string = [options_string ' --I_avg'];
    end
    
    if ismember('I', sensor.record)
        disp('WARNING: output parameter ''I'' is not supported by the C++ code and will not be returned.');
    end
    
end
    
% check if sensor.record_start_imdex is given
if isfield(sensor, 'record_start_index')
    options_string = [options_string ' -s ' num2str(sensor.record_start_index)];
end

% extract the optional input arguments 
if nargin > 4
    input_args = varargin(5:end);
else
    input_args = {};
end
 
% append the save to disk parameter
input_args = [input_args {'SaveToDisk', input_filename}];

% run the MATLAB code first to generate the input file and save to disk
kspaceFirstOrder3D(varargin{1:2}, source, sensor, input_args{:});
 
% get the location of this m-file to locate the C++ executable (the
% executables are assumed to be in the binaries folder)
path = mfilename('fullpath');
path = path(1:end-length(mfilename));

% run the simulation in C++ and print outputs to the MATLAB command line
if isunix
    
    % prepend spaces in linux pathnames with \ to allow cd to work
    path = strrep(path, ' ', '\ ');
    
    % clear the library path to prevent a conflict with the FFTW libraries
    % loaded automatically by MATLAB, and run linux binary
    run_string = ['system(''export LD_LIBRARY_PATH=; cd ' path 'binaries; ./kspaceFirstOrder3D-OMP -i ' input_filename ' -o ' output_filename options_string ''' ,''-echo'');'];
    eval(run_string);
    
else
    % run Windows binary, placing the input and output filenames in double
    % quotations to avoid problems with spaces
    run_string = ['system(''cd ' path 'binaries & kspaceFirstOrder3D-OMP.exe -i "' input_filename '" -o "' output_filename '" ' options_string ''' ,''-echo'');'];
    eval(run_string);
end

% load the C++ data back from disk using h5read
if time_rev
    sensor_data = h5read(output_filename, '/p_final');
elseif isfield(sensor, 'record')    
    if ismember('p', sensor.record)
        sensor_data.p = h5read(output_filename, '/p');
    end
    if ismember('p_max', sensor.record)
        sensor_data.p_max = h5read(output_filename, '/p_max');
    end
    if ismember('p_rms', sensor.record)
        sensor_data.p_rms = h5read(output_filename, '/p_rms');
    end    
    if ismember('p_final', sensor.record)
        sensor_data.p_final = h5read(output_filename, '/p_final');
    end    
    if ismember('u', sensor.record)
        sensor_data.ux = h5read(output_filename, '/ux');
        sensor_data.uy = h5read(output_filename, '/uy');
        sensor_data.uz = h5read(output_filename, '/uz');
    end
    if ismember('u_max', sensor.record)
        sensor_data.ux_max = h5read(output_filename, '/ux_max');
        sensor_data.uy_max = h5read(output_filename, '/uy_max');
        sensor_data.uz_max = h5read(output_filename, '/uz_max');
    end
    if ismember('u_rms', sensor.record)
        sensor_data.ux_rms = h5read(output_filename, '/ux_rms');
        sensor_data.uy_rms = h5read(output_filename, '/uy_rms');
        sensor_data.uz_rms = h5read(output_filename, '/uz_rms');
    end 
    if ismember('u_final', sensor.record)
        sensor_data.ux_final = h5read(output_filename, '/ux_final');
        sensor_data.uy_final = h5read(output_filename, '/uy_final');
        sensor_data.uz_final = h5read(output_filename, '/uz_final');
    end    
    if ismember('I_avg', sensor.record)
        sensor_data.Ix_avg = h5read(output_filename, '/Ix_avg');
        sensor_data.Iy_avg = h5read(output_filename, '/Iy_avg');
        sensor_data.Iz_avg = h5read(output_filename, '/Iz_avg');
    end
else
    sensor_data = h5read(output_filename, '/p');
end

% delete the input and output files
delete(input_filename);
delete(output_filename);