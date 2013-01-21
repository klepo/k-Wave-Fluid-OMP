% DESCRIPTION:
%       subscript to create storage variables
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 1st August 2011
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

% =========================================================================
% PREPARE DATA MASKS AND STORAGE VARIABLES
% =========================================================================

if use_sensor
    
    % check sensor mask based on the Cartesian interpolation setting
    if ~binary_sensor_mask && strcmp(cartesian_interp, 'nearest')
        
        % extract the data using the binary sensor mask created in
        % inputChecking, but switch on Cartesian reorder flag so that the
        % final data is returned in the correct order (not in time
        % reversal mode).
        binary_sensor_mask = true;
        if ~time_rev
            reorder_data = true;
        end

        % check if any duplicate points have been discarded in the
        % conversion from a Cartesian to binary mask
        num_discarded_points = length(sensor_x) - sum(sensor.mask(:));
        if num_discarded_points ~= 0
            disp(['  WARNING: ' num2str(num_discarded_points) ' duplicated sensor points discarded (nearest neighbour interpolation)']);
        end        

    end
          
    % preallocate output variables
    if ~time_rev

        % get the number of sensor points
        if binary_sensor_mask
            num_sensor_points = sum(sensor.mask(:));
        else
            num_sensor_points = length(sensor_x);
        end

        % calculate the number of time points that are stored.
        % - if streaming data to disk, reduce to the size of the
        %   sensor_data matrix based on the value of stream_to_disk
        % - if a user input for sensor.record_start_index is given, reduce
        %   the size of the sensor_data matrix based on the value given
        if kgrid.dim == 3 && stream_to_disk
            num_recorded_time_points = stream_to_disk;

            % initialise the file index variable
            stream_data_index = 1;
        else
            num_recorded_time_points = kgrid.Nt - sensor.record_start_index + 1;
        end

        % create storage and scaling variables - all variables are saved as
        % fields of a structure called sensor_data. If only p is being stored
        % (i.e., if no user input is given for sensor.record), then
        % sensor_data.p is copied to sensor_data at the end of the simulation

        % time history of the acoustic pressure
        if record.p
            sensor_data.p = castZeros([num_sensor_points, num_recorded_time_points]);
        end

        % maximum pressure
        if record.p_max
            sensor_data.p_max = castZeros([num_sensor_points, 1]);
        end

        % rms pressure
        if record.p_rms
            sensor_data.p_rms = castZeros([num_sensor_points, 1]);
        end

        % time history of the acoustic particle velocity
        if record.u
            % pre-allocate the velocity fields based on the number of
            % dimensions in the simulation
            switch kgrid.dim
                case 1
                    sensor_data.ux = castZeros([num_sensor_points, num_recorded_time_points]);
                case 2
                    sensor_data.ux = castZeros([num_sensor_points, num_recorded_time_points]);
                    sensor_data.uy = castZeros([num_sensor_points, num_recorded_time_points]);
                case 3
                    sensor_data.ux = castZeros([num_sensor_points, num_recorded_time_points]);
                    sensor_data.uy = castZeros([num_sensor_points, num_recorded_time_points]);
                    sensor_data.uz = castZeros([num_sensor_points, num_recorded_time_points]); 
            end                
        end

        % maximum particle velocity
        if record.u_max
            % pre-allocate the velocity fields based on the number of
            % dimensions in the simulation
            switch kgrid.dim
                case 1            
                    sensor_data.ux_max = castZeros([num_sensor_points, 1]);
                case 2
                    sensor_data.ux_max = castZeros([num_sensor_points, 1]);
                    sensor_data.uy_max = castZeros([num_sensor_points, 1]);
                case 3
                    sensor_data.ux_max = castZeros([num_sensor_points, 1]);
                    sensor_data.uy_max = castZeros([num_sensor_points, 1]);
                    sensor_data.uz_max = castZeros([num_sensor_points, 1]);
            end
        end

        % rms particle velocity
        if record.u_rms
            % pre-allocate the velocity fields based on the number of
            % dimensions in the simulation
            switch kgrid.dim
                case 1            
                    sensor_data.ux_rms = castZeros([num_sensor_points, 1]);
                case 2
                    sensor_data.ux_rms = castZeros([num_sensor_points, 1]);
                    sensor_data.uy_rms = castZeros([num_sensor_points, 1]);
                case 3
                    sensor_data.ux_rms = castZeros([num_sensor_points, 1]);
                    sensor_data.uy_rms = castZeros([num_sensor_points, 1]);
                    sensor_data.uz_rms = castZeros([num_sensor_points, 1]);
            end
        end

        % instantaneus acoustic intensity
        if record.I
            % pre-allocate the intensity fields based on the number of
            % dimensions in the simulation
            switch kgrid.dim
                case 1
                    sensor_data.Ix = castZeros([num_sensor_points, num_recorded_time_points]);
                case 2
                    sensor_data.Ix = castZeros([num_sensor_points, num_recorded_time_points]);
                    sensor_data.Iy = castZeros([num_sensor_points, num_recorded_time_points]);
                case 3
                    sensor_data.Ix = castZeros([num_sensor_points, num_recorded_time_points]);
                    sensor_data.Iy = castZeros([num_sensor_points, num_recorded_time_points]);
                    sensor_data.Iz = castZeros([num_sensor_points, num_recorded_time_points]); 
            end               
        end

        % average acoustic intensity
        if record.I_avg
            % pre-allocate the intensity fields based on the number of
            % dimensions in the simulation
            switch kgrid.dim
                case 1            
                    sensor_data.Ix_avg = castZeros([num_sensor_points, 1]);
                case 2
                    sensor_data.Ix_avg = castZeros([num_sensor_points, 1]);
                    sensor_data.Iy_avg = castZeros([num_sensor_points, 1]);
                case 3
                    sensor_data.Ix_avg = castZeros([num_sensor_points, 1]);
                    sensor_data.Iy_avg = castZeros([num_sensor_points, 1]);
                    sensor_data.Iz_avg = castZeros([num_sensor_points, 1]);
            end
        end

        % temporary buffer to store p, ux, uy, uz at the previous time step
        % if storing the intensity
        if (record.I || record.I_avg)
            % pre-allocate the velocity fields based on the number of
            % dimensions in the simulation
            switch kgrid.dim
                case 1            
                    sensor_data.ux_prev_t = castZeros([num_sensor_points, 1]);
                case 2
                    sensor_data.ux_prev_t = castZeros([num_sensor_points, 1]);
                    sensor_data.uy_prev_t = castZeros([num_sensor_points, 1]);
                case 3
                    sensor_data.ux_prev_t = castZeros([num_sensor_points, 1]);
                    sensor_data.uy_prev_t = castZeros([num_sensor_points, 1]);
                    sensor_data.uz_prev_t = castZeros([num_sensor_points, 1]);
            end
            sensor_data.p_prev_t = castZeros([num_sensor_points, 1]);
        end
        
        % object of the kWaveTransducer class is being used as a sensor            
        if transducer_sensor
            if transducer_receive_elevation_focus
                % if there is elevation focusing, a buffer is
                % needed to store a short time history at each
                % sensor point before averaging
                sensor_data_buffer_size = max(sensor.elevation_beamforming_delays) + 1;
                if sensor_data_buffer_size > 1
                    sensor_data_buffer = castZeros([num_sensor_points, sensor_data_buffer_size]);
                else
                    clear sensor_data_buffer sensor_data_buffer_size;
                    transducer_receive_elevation_focus = false;
                end
            end

            % the grid points can be summed on the fly and so the
            % sensor is the size of the number of active elements 
            sensor_data.transducer = castZeros([sensor.number_active_elements, num_recorded_time_points]);
        end

        % precomputate the triangulation points if a Cartesian sensor mask
        % is used with linear interpolation (tri and bc are the Delaunay
        % triangulation and Barycentric coordinates)
        if ~binary_sensor_mask
            if kgrid.dim == 1
                
                % assign pseudonym for Cartesain grid points in 1D (this
                % is later used for data casting)
                grid_x = kgrid.x_vec;
                
            else

                % update command line status
                disp('  calculating Delaunay triangulation...');

                % compute triangulation
                if kgrid.dim == 2
                    [tri, bc] = gridDataFast2D(kgrid.x, kgrid.y, sensor_x, sensor_y);
                elseif kgrid.dim == 3
                    [tri, bc] = gridDataFast3D(kgrid.x, kgrid.y, kgrid.z, sensor_x, sensor_y, sensor_z);
                end
                
            end
        end

        % if recording the intensity, the particle velocity must be intepolated
        % onto the regular grid before multiplying by the pressure. For a
        % binary sensor mask, this is done by finding the linear index for the
        % neighbouring x, y, and z points for each element of the sensor mask.
        % The velocity on the regular grid points can then be calculated from
        % the average of two adjacent staggered grid points.
        if record.I_avg || record.I
            if binary_sensor_mask
                switch kgrid.dim   
                    case 1
                        
                        % get the indices of grid points offset in each direction
                        sensor_mask_index_sgx = sensor_mask_index - 1;
                        sensor_mask_index_sgx(sensor_mask_index_sgx < 1) = 1;
                        
                    case 2
                        % get x, y indices of the sensor mask
                        sz = [kgrid.Nx, kgrid.Ny];
                        [indx, indy] = ind2sub(sz, sensor_mask_index);

                        % get the indices of grid points offset in each direction
                        indx_off = indx - 1;
                        indx_off(indx_off < 1) = 1;
                        sensor_mask_index_sgx = sub2ind(sz, indx_off, indy);
                        indy_off = indy - 1;
                        indy_off(indy_off < 1) = 1;                    
                        sensor_mask_index_sgy = sub2ind(sz, indx, indy_off);

                        % cleanup unused variables
                        clear sz indx indy indx_off indy_off;                
                    case 3
                        % get x, y, z indices of the sensor mask
                        sz = [kgrid.Nx, kgrid.Ny, kgrid.Nz];
                        [indx, indy, indz] = ind2sub(sz, sensor_mask_index);

                        % get the indices of grid points offset in each direction
                        indx_off = indx - 1;
                        indx_off(indx_off < 1) = 1;
                        sensor_mask_index_sgx = sub2ind(sz, indx_off, indy, indz);
                        indy_off = indy - 1;
                        indy_off(indy_off < 1) = 1;                    
                        sensor_mask_index_sgy = sub2ind(sz, indx, indy_off, indz);
                        indz_off = indz - 1;
                        indz_off(indz_off < 1) = 1;                    
                        sensor_mask_index_sgz = sub2ind(sz, indx, indy, indz_off);

                        % cleanup unused variables
                        clear sz indx indy indz indx_off indy_off indz_off;
                end       
            elseif kgrid.dim ~= 1;

                % update command line status
                disp('  calculating Delaunay triangulation for Intensity calculation...');

                % precompute triangulation points for interpolating the particle
                % velocity at the sensor positions onto the regular grid
                switch kgrid.dim
                    case 2
                        [tri_sgx, bc_sgx] = gridDataFast2D(kgrid.x, kgrid.y, sensor_x - kgrid.dx/2, sensor_y);
                        [tri_sgy, bc_sgy] = gridDataFast2D(kgrid.x, kgrid.y, sensor_x, sensor_y - kgrid.dy/2);
                    case 3
                        [tri_sgx, bc_sgx] = gridDataFast3D(kgrid.x, kgrid.y, kgrid.z, sensor_x - kgrid.dx/2, sensor_y, sensor_z);
                        [tri_sgy, bc_sgy] = gridDataFast3D(kgrid.x, kgrid.y, kgrid.z, sensor_x, sensor_y - kgrid.dy/2, sensor_z);
                        [tri_sgz, bc_sgz] = gridDataFast3D(kgrid.x, kgrid.y, kgrid.z, sensor_x, sensor_y, sensor_z - kgrid.dz/2);
                end
            end
        end  
    end
end