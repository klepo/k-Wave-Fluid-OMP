%%
% author: Bradley Treeby
% date: 11th Feb 2014

clear all;

% add toolbox
addpath('../../../k-wave-matlab');

SAVE_TO_DISK = false;
USE_CUBOID = true;

sensor.record = {'p', 'p_max','p_min','p_rms'};

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 64;            % number of grid points in the x direction
Ny = 64;            % number of grid points in the y direction
Nz = 64;            % number of grid points in the z direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
dz = 0.1e-3;        % grid point spacing in the z direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

% define the properties of the propagation medium
medium.sound_speed = 1500*ones(Nx, Ny, Nz);	% [m/s]
medium.sound_speed(1:Nx/2, :, :) = 1800;    % [m/s]
medium.density = 1000*ones(Nx, Ny, Nz);     % [kg/m^3]
medium.density(:, Ny/4:end, :) = 1200;      % [kg/m^3]

% create initial pressure distribution using makeBall
ball_magnitude = 10;    % [Pa]
ball_x_pos = 38;        % [grid points]
ball_y_pos = 32;        % [grid points]
ball_z_pos = 32;        % [grid points]
ball_radius = 5;        % [grid points]
ball_1 = ball_magnitude*makeBall(Nx, Ny, Nz, ball_x_pos, ball_y_pos, ball_z_pos, ball_radius);

ball_magnitude = 10;    % [Pa]
ball_x_pos = 20;        % [grid points]
ball_y_pos = 20;        % [grid points]
ball_z_pos = 20;        % [grid points]
ball_radius = 3;        % [grid points]
ball_2 = ball_magnitude*makeBall(Nx, Ny, Nz, ball_x_pos, ball_y_pos, ball_z_pos, ball_radius);

source.p0 = ball_1 + ball_2;

if USE_CUBOID
    % define list of cuboid corners
    sensor.mask = [40, 30, 30, 50, 40, 40; 10, 30, 30, 20, 40, 40; 1, 1, 12, 64, 64, 13; 7,13,15,25,27,39].';    
    %sensor.mask = [10, 30, 30, 20, 40, 40].';
else
    % define binary mask
    sensor.mask = zeros(Nx, Ny, Nz);
    sensor.mask(:, Ny/2, :) = 1;
end

sensor_data_matlab = kspaceFirstOrder3D(kgrid, medium, source, sensor, 'DataCast', 'single');
sensor_data_cxx = kspaceFirstOrder3DC(kgrid, medium, source, sensor, 'DataCast', 'single', 'BinaryPath', '/home/jarosjir/Devices/kazi/Sources/K-Wave++/k-wave-fluid-omp');
%kspaceFirstOrder3D(kgrid, medium, source, sensor, 'SaveToDisk', 'cuboid_corners.h5');
% run the simulation
%if ~SAVE_TO_DISK
%    sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, 'DataCast', 'single');
%else
%    %kspaceFirstOrder3D(kgrid, medium, source, sensor, 'SaveToDisk', 'cuboid_corners.h5');
%    sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, 'DataCast', 'single', 'BinaryPath', '/home/jarosjir/Devices/kazi/Sources/K-Wave++/k-wave-fluid-omp');
%    
%   %return
%end


%%
% =========================================================================
% VISUALISATION
% =========================================================================

data_matlab = sensor_data_matlab.p;
data_cxx    = sensor_data_cxx.p;

%plot the simulated sensor data
figure;
imagesc(data_matlab, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;

figure;
imagesc(data_cxx, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;
    
if size(data_matlab,1) *size(data_matlab,2) ~= size(data_cxx,1) *size(data_cxx,2)    
    disp('Number of ellements is not the same')
end

y1 = reshape(data_matlab, size(data_matlab,1) *size(data_matlab,2) , []);
y2 = reshape(data_cxx,    size(data_cxx,1)    *size(data_cxx,2), []);
y2 = y2./max(y1(:));
y1 = y1./max(y1(:));
figure; 
subplot(2, 1, 1), plot(y1.', 'k-');
hold on;
plot(y2.', 'r-');
subplot(2, 1, 2), plot((y1 - y2).');