% Script to compare MATLAB and C++ versions of 3D code.
%
% author: Bradley Treeby
% date: 16th February 2012
%
% update: 20th June 2012 for version 2.10
% update: 15th August 2012 for version 2.12
% update: 19th August 2012 to include linear / lossless
% update: 28th August 2012 to automatically include permutations

clear all;
addpath('k-Wave');

% =========================================================================
% SIMULATION LITERALS
% =========================================================================

% set total number of grid points not including the PML
NX = 128;    % [grid points]
NY = 64;     % [grid points]
NZ = 32;     % [grid points]

% option to use non-uniform grid
USE_NONUNIFORM_GRID = false;

% medium settings
HETEROGENEOUS_RHO0   = true;
HETEROGENEOUS_C0     = true;
HETEROGENEOUS_BONA   = true;      % not used if NONLINEAR = false
HETEROGENEOUS_ALPHA0 = true;    % not used if ABSORBING = false

% sensor settings
MASK_PLANE = 'xy';     % set to 'xy' or 'xz' to generate the beam pattern in different planes
RECORD_VELOCITY = false;
RECORD_P_AVG = false;
RECORD_INTENSITY = false;

% source type
% 1: transducer source
% 2: pressure source
% 3: velocity-x source
% 4: velocity-y source
% 5: velocity-z source
% 6: velocity-x/y/z source
% 7: initial pressure
SOURCE_TYPE = 7;    

% option to generate a single or multiple time series
SOURCE_MANY = false;

% calculation settings
DATA_CAST = 'single';  % set to 'single' or 'gsingle' to speed up the matlab computations
NONLINEAR = true;
ABSORBING = true;

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% define the properties of the propagation medium
c0 = 1540;      % [m/s]
rho0 = 1000;    % [kg/m^3]
alpha0 = 0.75;  % [dB/(MHz^y cm)]
y = 1.5;
BonA = 6;

% create heterogeneous region
sc = 1.3;
heterog = makeBall(NX, NY, NZ, NX/2, NY/2, NZ/2, 30);

% set density
if HETEROGENEOUS_RHO0
    medium.density = rho0*ones(NX, NY, NZ);
    medium.density(heterog == 1) = rho0*sc;
else
    medium.density = rho0; 
end

% set sound speed
if HETEROGENEOUS_C0
    medium.sound_speed = c0*ones(NX, NY, NZ);
    medium.sound_speed(heterog == 1) = c0*sc;
else
    medium.sound_speed = c0;
end
    
% set BonA
if NONLINEAR
    if HETEROGENEOUS_BONA
        medium.BonA = BonA*ones(NX, NY, NZ);
        medium.BonA(heterog == 1) = BonA*sc;
    else
        medium.BonA = BonA;
    end
end

% set absorption terms
if ABSORBING
    if HETEROGENEOUS_ALPHA0
        medium.alpha_coeff  = alpha0*ones(NX, NY, NZ);
        medium.alpha_coeff(heterog == 1) = alpha0*sc;
        medium.alpha_power = y;
    else
        medium.alpha_coeff = alpha0;      
        medium.alpha_power = y;
    end
end

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
PML_X_SIZE = 20;            % [grid points]
PML_Y_SIZE = 10;            % [grid points]
PML_Z_SIZE = 10;            % [grid points]

% set desired grid size in the x-direction not including the PML
x = 40e-3;                  % [m]

% calculate the spacing between the grid points
dx = x/NX;                  % [m]
dy = dx;                    % [m]
dz = dx;                    % [m]

% create the k-space grid
kgrid = makeGrid(NX, dx, NY, dy, NZ, dz);

% =========================================================================
% DEFINE NONUNIFORM GRID SETTINGS
% =========================================================================

if USE_NONUNIFORM_GRID
    
    % non-uniform grid settings
    gamma = 0;      % position of staggering
    beta = 0.15;     % clustering parameter where abs(beta) < 1
    x = kgrid.x_vec - kgrid.x_vec(1);
    
    % create nonuniform grid gradient
    dxu = 2*pi/(NX - 1);
    xu = 0:dxu:2*pi;
    dxudxn = (1 + beta^2 - 2*beta*cos(xu - gamma))./(1 - beta^2);
    
    dyu = 2*pi/(NY - 1);
    yu = 0:dyu:2*pi;
    dyudyn = (1 + beta^2 - 2*beta*cos(yu - gamma))./(1 - beta^2);
    
    dzu = 2*pi/(NZ - 1);
    zu = 0:dzu:2*pi;
    dzudzn = (1 + beta^2 - 2*beta*cos(zu - gamma))./(1 - beta^2);    
    
    % grid point mapping
    xn = (exp(1i*(xu - gamma)) - beta) ./ (1 - beta*exp(1i*(xu - gamma)));
    xn = unwrap(real(log(xn)/1i));
    xn = xn - min(xn(:));
    xn = xn .* (kgrid.x_size / xn(end));
    
    % plot to check
    figure;
    subplot(2, 1, 1), plot(xn);
    subplot(2, 1, 2), plot(dxudxn);
    
    % assign non-uniform grid properties
    kgrid.dxudxn = dxudxn;
    kgrid.dyudyn = dyudyn;
    kgrid.dzudzn = dzudzn;
    
    % manually scale source based on dx (average in each direction)
    % Scaling is: 1/c^2 then 1/N where N = 3 then 2*dt*c/dx   
    av_dx = (xn(PML_X_SIZE + 2) - xn(PML_X_SIZE))/2;
    av_spacing = (av_dx + dy)/2;
    input_signal = input_signal .* (2*kgrid.dt./(3*medium.sound_speed*av_spacing));
    
end

% =========================================================================
% DEFINE THE TIME ARRAY
% =========================================================================

% create the time array
t_end = 45e-6;                  % [s]
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

% 1: transducer source
% 2: pressure source
% 3: velocity-x source
% 4: velocity-y source
% 5: velocity-z source
% 6: velocity-x/y/z source
% 7: initial pressure

switch SOURCE_TYPE
    case 1
        
        % -----------------------------
        % OPTION 1: Transducer Source
        % -----------------------------
        
        % scale the source magnitude by the source_strength divided by the
        % impedance (the source is assigned to the particle velocity)
        input_signal = input_signal./(c0*rho0);

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
        transducer.position = round([PML_X_SIZE + 1, NY/2 - transducer_width/2, NZ/2 - transducer.element_length/2]);

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
        
    case 2
        
        % -----------------------------
        % OPTION 2: Pressure Source
        % -----------------------------
        
        % set the source mask to be a small rectangle
        source_radius = 10;
        source.p_mask = zeros(NX, NY, NZ);
        source.p_mask(PML_X_SIZE + 1, NY/2 - source_radius + 1:NY/2 + source_radius,  NZ/2 - source_radius/2 + 1:NZ/2 + source_radius/2) = 1;

        % assign the source term
        if SOURCE_MANY
            source.p = focus(input_signal, squeeze(source.p_mask(PML_X_SIZE + 1, :, :)), 10*dx, c0, dx, kgrid.dt);
        else
            source.p = input_signal;
        end
        
    case 3
        
        % -----------------------------
        % OPTION 3: Velocity X Source
        % -----------------------------
        
        % set the source mask to be a small rectangle
        source_radius = 10;
        source.u_mask = zeros(NX, NY, NZ);
        source.u_mask(PML_X_SIZE + 1, NY/2 - source_radius + 1:NY/2 + source_radius,  NZ/2 - source_radius/2 + 1:NZ/2 + source_radius/2) = 1;

        % assign the source term scaled by the impedance
        if SOURCE_MANY
            source.ux = focus(input_signal./(c0*rho0), squeeze(source.u_mask(PML_X_SIZE + 1, :, :)), 10*dx, c0, dx, kgrid.dt);
        else
            source.ux = input_signal./(c0*rho0);
        end
        
    case 4
        
        % -----------------------------
        % OPTION 4: Velocity Y Source
        % -----------------------------
        
        % set the source mask to be a small rectangle
        source_radius = 10;
        source.u_mask = zeros(NX, NY, NZ);
        source.u_mask(NX/2 - source_radius + 1:NX/2 + source_radius, PML_Y_SIZE + 1, NZ/2 - source_radius + 1:NZ/2 + source_radius) = 1;

        % assign the source term scaled by the impedance
        if SOURCE_MANY
            source.uy = focus(input_signal./(c0*rho0), squeeze(source.u_mask(:, PML_Y_SIZE + 1, :)), 10*dx, c0, dx, kgrid.dt);
        else
            source.uy = input_signal./(c0*rho0);
        end
        
    case 5
        
        % -----------------------------
        % OPTION 5: Velocity Z Source
        % -----------------------------        
        
        % set the source mask to be a small rectangle
        source_radius = 10;
        source.u_mask = zeros(NX, NY, NZ);
        source.u_mask(NX/2 - source_radius + 1:NX/2 + source_radius, NY/2 - source_radius + 1:NY/2 + source_radius, PML_Z_SIZE + 1) = 1;

        % assign the source term scaled by the impedance
        if SOURCE_MANY
            source.uz = focus(input_signal./(c0*rho0), squeeze(source.u_mask(:, :, PML_Z_SIZE + 1)), 10*dx, c0, dx, kgrid.dt);
        else
            source.uz = input_signal./(c0*rho0);
        end
        
    case 6
        
        % -----------------------------
        % OPTION 6: Velocity XYZ Source
        % ----------------------------- 
        
        
    case 7
        
        % -----------------------------
        % OPTION 7: Initial Pressure
        % ----------------------------- 
        
        % create ball shaped initial pressure distribution
        source.p0 = source_strength*makeBall(NX, NY, NZ, NX/2, NY/2, NZ/2, 5);
        
end

% =========================================================================
% DEFINE SENSOR MASK
% =========================================================================

% define a sensor mask through the central plane
sensor.mask = zeros(NX, NY, NZ);
switch MASK_PLANE
    case 'xy'
        % define mask
        sensor.mask(:, :, NZ/2) = 1;
        
        % store y axis properties        
        Nj = NY;
        j_vec = kgrid.y_vec;
        j_label = 'y';
        
    case 'xz'
        % define mask
        sensor.mask(:, NY/2, :) = 1;
        
        % store z axis properties
        Nj = NZ;
        j_vec = kgrid.z_vec;
        j_label = 'z';
        
    case 'statistics'
        % define mask 
        sensor.mask = ones(NX, NY, NZ);
        sensor.record_mode = 'statistics';
end 

% =========================================================================
% RUN THE MATLAB SIMULATION
% =========================================================================

% set the input settings
input_args = {'DisplayMask', 'off', 'PlotPML', false, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE], ...
    'DataCast', DATA_CAST, 'PlotScale', [-source_strength/2, source_strength/2]};

% return velocity
if RECORD_VELOCITY
    input_args = [input_args {'ReturnVelocity', true}];
end

% run the simulation using k-Wave
[sensor_data] = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% RUN THE C++ SIMULATION
% =========================================================================

% run the simulation again using Jiri's C++ version
[sensor_data_cpp] = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% COMPARE THE ACCURACY OF THE BEAM PATTERNS
% =========================================================================

if strcmp(MASK_PLANE, 'statistics');
    
    % trim the PML from the c++ data
    sensor_data_cpp.p_max = sensor_data_cpp.p_max(1 + PML_X_SIZE:end - PML_X_SIZE, 1 + PML_Y_SIZE:end - PML_Y_SIZE, 1 + PML_Z_SIZE:end - PML_Z_SIZE);
    
    % reshape the MATLAB data
    sensor_data.p_max = reshape(sensor_data.p_max, [NX, NY, NZ]);
    
    % plot the beam patterns
    figure;
    subplot(1, 3, 1);
    beamPlot(sensor_data.p_max);
    colorbar;
    title('k-Wave');
    
    % plot the beam patterns
    subplot(1, 3, 2);
    beamPlot(sensor_data_cpp.p_max);
    colorbar;
    title('C++');
    
    % plot the error
    subplot(1, 3, 3);
    beamPlot(100*abs(sensor_data.p_max - sensor_data_cpp.p_max)./max(abs(sensor_data.p_max(:))));
    colorbar;
    title('Global Error [%]');
    
else
    
    % reshape the sensor data to its original position so that it can be
    % indexed as sensor_data(x, j, t)
    if RECORD_VELOCITY
        sensor_data.p = reshape(sensor_data.p, [NX, Nj, length(kgrid.t_array)]);
        sensor_data.ux = reshape(sensor_data.ux, [NX, Nj, length(kgrid.t_array)]);
        sensor_data.uy = reshape(sensor_data.uy, [NX, Nj, length(kgrid.t_array)]);
        sensor_data.uz = reshape(sensor_data.uz, [NX, Nj, length(kgrid.t_array)]);
        sensor_data_cpp.p = reshape(sensor_data_cpp.p, [NX, Nj, length(kgrid.t_array)]);
        sensor_data_cpp.ux = reshape(sensor_data_cpp.ux, [NX, Nj, length(kgrid.t_array)]);
        sensor_data_cpp.uy = reshape(sensor_data_cpp.uy, [NX, Nj, length(kgrid.t_array)]);
        sensor_data_cpp.uz = reshape(sensor_data_cpp.uz, [NX, Nj, length(kgrid.t_array)]);
    else
        sensor_data = reshape(sensor_data, [NX, Nj, length(kgrid.t_array)]);
        sensor_data_cpp = reshape(sensor_data_cpp, [NX, Nj, length(kgrid.t_array)]);
    end

    % compute the frequency axis
    freq = (0:ceil((kgrid.Nt + 1)/2) - 1) ./ (kgrid.dt*length(kgrid.t_array));

    % compute the index at which the source frequency and its harmonics occur
    [f0_value, f0_index] = findClosest(freq, tone_burst_freq);
    [f1_value, f1_index] = findClosest(freq, tone_burst_freq*2);

    % preallocate the beam pattern variables
    beam_pattern_f0 = zeros(NX, Nj);
    beam_pattern_f1 = zeros(NX, Nj);
    beam_pattern_total = zeros(NX, Nj);
    beam_pattern_f0_cpp = zeros(NX, Nj);
    beam_pattern_f1_cpp = zeros(NX, Nj);
    beam_pattern_total_cpp = zeros(NX, Nj);

    % compute the amplitude spectrum of the time series recorded at each sensor
    % point, and then extract the corresponding amplitudes at the fundamental
    % frequency and second harmonic.
    for x_index = 1:NX
        for j_index = 1:Nj

            % compute the amplitude spectrum
            if RECORD_VELOCITY
                amp_spect = spectrum(squeeze(sensor_data.p(x_index, j_index, :)), 1/kgrid.dt); %, 'Window', 'Hanning');
                amp_spect_cpp = spectrum(squeeze(sensor_data_cpp.p(x_index, j_index, :)), 1/kgrid.dt);
            else
                amp_spect = spectrum(squeeze(sensor_data(x_index, j_index, :)), 1/kgrid.dt); %, 'Window', 'Hanning');
                amp_spect_cpp = spectrum(squeeze(sensor_data_cpp(x_index, j_index, :)), 1/kgrid.dt);
            end

            % extract the amplitude at the source frequency and store
            beam_pattern_f0(x_index, j_index) = amp_spect(f0_index);
            beam_pattern_f0_cpp(x_index, j_index) = amp_spect_cpp(f0_index);

            % extract the amplitude at the source frequency and store
            beam_pattern_f1(x_index, j_index) = amp_spect(f1_index);      
            beam_pattern_f1_cpp(x_index, j_index) = amp_spect_cpp(f1_index);

            % extract the integral of the total amplitude spectrum
            beam_pattern_total(x_index, j_index) = sum(amp_spect(:));
            beam_pattern_total_cpp(x_index, j_index) = sum(amp_spect_cpp(:));

        end
    end

    for plot_index = 1:(3 + 3*double(RECORD_VELOCITY))

        switch plot_index
            case 1
                % plot of total beam pattern
                k_wave = beam_pattern_total;
                cpp = beam_pattern_total_cpp;
                plot_title = ' (Total)';
            case 2
                % plot of fundamental
                k_wave = beam_pattern_f0;
                cpp = beam_pattern_f0_cpp;          
                plot_title = ' (Fundamental)';
            case 3
                % plot of second harmonic
                k_wave = beam_pattern_f1;
                cpp = beam_pattern_f1_cpp; 
                plot_title = ' (2nd Harmonic)';
            case 4
                k_wave = max(sensor_data.ux, [], 3);
                cpp = max(sensor_data_cpp.ux, [], 3);
                plot_title = ' (ux)';
            case 5
                k_wave = max(sensor_data.uy, [], 3);
                cpp = max(sensor_data_cpp.uy, [], 3);
                plot_title = ' (uy)';            
            case 6
                k_wave = max(sensor_data.uz, [], 3);
                cpp = max(sensor_data_cpp.uz, [], 3);
                plot_title = ' (uz)';            
        end

        figure;
        subplot(1, 4, 1);
        imagesc(j_vec*1e3, (kgrid.x_vec - min(kgrid.x_vec(:)))*1e3, k_wave/1e6);
        xlabel([j_label '-position [mm]']);
        ylabel('x-position [mm]');
        title(['k-Wave' plot_title]);
        colormap(jet(256));
        c = colorbar;
        ylabel(c, 'Pressure [MPa]');
        axis image;

        subplot(1, 4, 2);
        imagesc(j_vec*1e3, (kgrid.x_vec - min(kgrid.x_vec(:)))*1e3, cpp/1e6);
        xlabel([j_label '-position [mm]']);
        ylabel('x-position [mm]');
        title(['C++' plot_title]);
        colormap(jet(256));
        c = colorbar;
        ylabel(c, 'Pressure [MPa]');
        axis image;

        subplot(1, 4, 3);
        imagesc(j_vec*1e3, (kgrid.x_vec - min(kgrid.x_vec(:)))*1e3, 100*abs(k_wave - cpp)./k_wave);
        xlabel([j_label '-position [mm]']);
        ylabel('x-position [mm]');
        title('Local Error');
        colormap(jet(256));
        c = colorbar;
        ylabel(c, '[%]');
        axis image;

        subplot(1, 4, 4);
        imagesc(j_vec*1e3, (kgrid.x_vec - min(kgrid.x_vec(:)))*1e3, 100*abs(k_wave - cpp)./max(abs(k_wave(:))));
        xlabel([j_label '-position [mm]']);
        ylabel('x-position [mm]');
        title('Global Error');
        colormap(jet(256));
        c = colorbar;
        ylabel(c, '[%]');
        axis image;

        scaleFig(2, 1);

    end
end