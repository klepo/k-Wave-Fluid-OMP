% DESCRIPTION:
%       subscript to cast loop variables
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 26th November 2010
%       last update - 2nd October 2012
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2012 Bradley Treeby and Ben Cox

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

% update command line status
disp(['  casting variables to ' data_cast ' type...']);
    
% create list of variable names to cast
cast_variables = {'kappa', 'c', 'dt', 'rho0'};

% create a separate list for indexing variables
cast_index_variables = {};

% add variables specific to simulations in certain dimensions
switch kgrid.dim
    case 1
        cast_variables = [cast_variables, {'ddx_k', 'shift_pos', 'shift_neg', ...
            'pml_x', 'pml_x_sgx',...
            'rho0_sgx'}];              
    case 2
        cast_variables = [cast_variables, {'ddx_k_shift_pos', 'ddx_k_shift_neg', 'ddy_k_shift_pos', 'ddy_k_shift_neg',...
            'pml_x', 'pml_y', 'pml_x_sgx', 'pml_y_sgy',...
            'rho0_sgx', 'rho0_sgy'}];
    case 3
        cast_variables = [cast_variables, {'ddx_k_shift_pos', 'ddy_k_shift_pos', 'ddz_k_shift_pos',... 
            'ddx_k_shift_neg', 'ddy_k_shift_neg', 'ddz_k_shift_neg',...
            'pml_x', 'pml_y', 'pml_z', 'pml_x_sgx', 'pml_y_sgy', 'pml_z_sgz',...        
            'rho0_sgx', 'rho0_sgy', 'rho0_sgz'}];          
end

% add sensor mask variables
if use_sensor
    cast_index_variables = [cast_index_variables, {'sensor_mask_index'}];
    if binary_sensor_mask && record.I
        switch kgrid.dim
            case 1
                cast_index_variables = [cast_index_variables, {'sensor_mask_index_sgx'}];
            case 2
                cast_index_variables = [cast_index_variables, {'sensor_mask_index_sgx', 'sensor_mask_index_sgy'}];
            case 3
                cast_index_variables = [cast_index_variables, {'sensor_mask_index_sgx', 'sensor_mask_index_sgy', 'sensor_mask_index_sgz'}];
        end
    end
end

% additional variables only used if the medium is absorbing
if strcmp(equation_of_state, 'absorbing')
    cast_variables = [cast_variables, {'absorb_nabla1', 'absorb_nabla2', 'absorb_eta', 'absorb_tau'}];
end

% additional variables only used if the propagation is nonlinear
if nonlinear
    cast_variables = [cast_variables, {'medium.BonA'}];
end

% additional variables only used if there is an initial pressure source
if isfield(source, 'p0')
    cast_variables = [cast_variables, {'source.p0'}];
end

% additional variables only used if there is a time varying pressure source term
if p_source
    cast_variables = [cast_variables, {'source.p'}];
    cast_index_variables = [cast_index_variables, {'p_source_index'}];
end

% additional variables only used if there is a time varying velocity source term
if ux_source || uy_source || uz_source
    cast_index_variables = [cast_index_variables, {'u_source_index'}];
end
if ux_source
    cast_variables = [cast_variables, {'source.ux'}];
end
if uy_source
    cast_variables = [cast_variables, {'source.uy'}];
end
if uz_source
    cast_variables = [cast_variables, {'source.uz'}];
end        

% addition variables only used if there is a transducer source
if transducer_source
    cast_variables = [cast_variables, {'transducer_input_signal'}];
    cast_index_variables = [cast_index_variables, {'u_source_index', 'delay_mask', 'transducer_source', 'transducer_transmit_apodization'}];
end

% addition variables only used if there is a transducer sensor with an
% elevation focus
if transducer_sensor && transducer_receive_elevation_focus
    cast_index_variables = [cast_index_variables, {'sensor_data_buffer', 'transducer_receive_mask'}];
end

% additional variables only used with nonuniform grids
if nonuniform_grid
    switch kgrid.dim
        case 1
            cast_index_variables = [cast_index_variables, {'kgrid.dxudxn'}];
        case 2
            cast_index_variables = [cast_index_variables, {'kgrid.dxudxn', 'kgrid.dyudyn'}];    
        case 3
            cast_index_variables = [cast_index_variables, {'kgrid.dxudxn', 'kgrid.dyudyn', 'kgrid.dzudzn'}];
    end
end

% additional variables only used for Cartesian sensor masks with linear
% interpolation
if use_sensor && ~binary_sensor_mask && ~time_rev
    if  kgrid.dim == 1
        cast_variables = [cast_variables, {'grid_x', 'sensor_x'}];
    else
        cast_variables = [cast_variables, {'tri', 'bc'}];
        if record.I || record.I_avg
            switch kgrid.dim
                case 2
                    cast_variables = [cast_variables, {'tri_sgx', 'bc_sgx', 'tri_sgy', 'bc_sgy'}];
                case 3
                    cast_variables = [cast_variables, {'tri_sgx', 'bc_sgx', 'tri_sgy', 'bc_sgy', 'tri_sgz', 'bc_sgz'}];
            end
        end
    end
end

% cast variables
for cast_index = 1:length(cast_variables)
    eval([cast_variables{cast_index} ' = ' data_cast '(' data_cast_prepend '(' cast_variables{cast_index} '));']);
end

% cast index variables only if casting to the GPU
if strncmp(data_cast, 'g', 1) || strncmp(data_cast, 'kWaveGPU', 8)
    for cast_index = 1:length(cast_index_variables)
        eval([cast_index_variables{cast_index} ' = ' data_cast '(' data_cast_prepend '(' cast_index_variables{cast_index} '));']);
    end
end