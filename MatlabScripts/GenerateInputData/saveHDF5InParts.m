function saveHDF5InParts(simulation_case, grid_size, filename)
% Function to save a HDF5 test file in parts
%
% Simulation case:
%
%       1:      nonlinear / absorbing
%               c0, rho0, BonA, alpha0 all heterogeneous
%       2:      nonlinear / absorbing
%               c0, rho0 heterogeneous / BonA, alpha0 scalar
%       3:      nonlinear / absorbing
%               c0, rho0, BonA, alpha0 all scalar
%       4:      linear / lossless
%               c0, rho0 heterogeneous
%       5:      linear / lossless
%               c0, rho0 scalar
%
% All simulations use a pressure source.
%
% All simulations use a binary sensor mask with default save option (p)
%
% author: Bradley Treeby
% date: 14th August 2012
% last update: 12th September 2012

start_time = clock;
disp('-------------------------------------------------------------------');
disp('WRITING HDF5 FILE IN PARTS');
disp('-------------------------------------------------------------------');

% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
% ========================================================================
% Filename
% ========================================================================
% \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

% remove file if it already exists
if exist(filename, 'file')
    delete(filename);
end

% load HDF5 constants
run([getkWavePath 'private/getH5Literals']);

% set simulation flags
switch simulation_case
    case 1
        flags.nonlinear         = true;
        flags.absorbing         = true;
        flags.c0_heterog        = true;
        flags.rho0_heterog      = true;
        flags.alpha0_heterog    = true;
        flags.BonA_heterog      = true;
    case 2
        flags.nonlinear         = true;
        flags.absorbing         = true;
        flags.c0_heterog        = true;
        flags.rho0_heterog      = true;
        flags.alpha0_heterog    = false;
        flags.BonA_heterog      = false;
    case 3
        flags.nonlinear         = true;
        flags.absorbing         = true;
        flags.c0_heterog        = false;
        flags.rho0_heterog      = false;
        flags.alpha0_heterog    = false;
        flags.BonA_heterog      = false;    
    case 4
        flags.nonlinear         = false;
        flags.absorbing         = false;
        flags.c0_heterog        = true;
        flags.rho0_heterog      = true;
        flags.alpha0_heterog    = false;
        flags.BonA_heterog      = false; 
    case 5
        flags.nonlinear         = false;
        flags.absorbing         = false;
        flags.c0_heterog        = false;
        flags.rho0_heterog      = false;
        flags.alpha0_heterog    = false;
        flags.BonA_heterog      = false; 
end

% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
% ========================================================================
% 1. Set and save grid properties
% ========================================================================
% \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

% set grid size variables (these are not cleared immediately as they are
% needed to calculate other things further down)
Nx_ref = grid_size(1);
Ny_ref = grid_size(2);
Nz_ref = grid_size(3);
Nt_ref = 1000;
dt_ref = 1e-8;
dx_ref = 1e-4;
dy_ref = 1e-4;
dz_ref = 1e-4;

% set vals
Nx = Nx_ref;
Ny = Ny_ref;
Nz = Nz_ref;
Nt = Nt_ref;
dt = dt_ref;
dx = dx_ref;
dy = dy_ref;
dz = dz_ref;

% :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::

% update command line status
tic; fprintf('Writing grid parameters... ');

% integer variables
variable_names = {'Nx', 'Ny', 'Nz', 'Nt'};

% change all the index variables to be in 64-bit unsigned integers (long in C++)
for index = 1:length(variable_names)

    % cast matrix to 64-bit unsigned integer
    eval([variable_names{index} ' = ' INTEGER_DATA_TYPE_MATLAB '(' variable_names{index} ');']);

    % write to HDF5 file
    writeMatrix(filename, eval(variable_names{index}), variable_names{index});

end

% clear holder
clear('variables_names');

% :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::

% float variables
variable_names = {'dt', 'dx', 'dy', 'dz'};

% change all the variables to be in single precision (float in C++), then
% add to HDF5 File
for index = 1:length(variable_names)

    % cast matrix to single precision
    eval([variable_names{index} ' = ' MATRIX_DATA_TYPE_MATLAB '(' variable_names{index} ');']);

    % write to HDF5 file
    writeMatrix(filename, eval(variable_names{index}), variable_names{index});

end

% clear holder
clear('variables_names');

% reset variables back in case they are needed later
Nx = Nx_ref;
Ny = Ny_ref;
Nz = Nz_ref;
Nt = Nt_ref;
dt = dt_ref;
dx = dx_ref;
dy = dy_ref;
dz = dz_ref;

% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
% ========================================================================
% 2. Set and save medium properties
% ========================================================================
% \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

% set material properties
rho0_scalar = 1000;
c0_scalar = 1500;
BonA_scalar = 6;
alpha_coeff_dB = 0.75;
alpha_power = 0.25;
heterog_sc = 1.2;
c_ref = c0_scalar*heterog_sc;

% :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::

% update command line status
toc; tic; fprintf('Writing c0... ');

% set sound speed
if flags.c0_heterog
    c0 = c0_scalar*ones(Nx, Ny, Nz, 'single');
    c0(Nx/2:end, :, :) = c0_scalar*heterog_sc;
else
    c0 = c0_scalar;
end

% cast matrix to single precision
eval(['c0 = ' MATRIX_DATA_TYPE_MATLAB '(c0);']);

% save
writeMatrix(filename, c0, 'c0');

% clear
clear('c0');

% :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::

% update command line status
toc; tic; fprintf('Writing rho0... ');

% set density
if flags.rho0_heterog
    rho0 = rho0_scalar*ones(Nx, Ny, Nz, 'single');
    rho0(Nx/2:end, :, :) = rho0_scalar*heterog_sc;
else
    rho0 = rho0_scalar;
end

% cast matrix to single precision
eval(['rho0 = ' MATRIX_DATA_TYPE_MATLAB '(rho0);']);

% save
writeMatrix(filename, rho0, 'rho0');

% update command line status
toc; tic; fprintf('Writing rho0_sgx... ');

% save
writeMatrix(filename, rho0, 'rho0_sgx');

% update command line status
toc; tic; fprintf('Writing rho0_sgy... ');

% save
writeMatrix(filename, rho0, 'rho0_sgy');

% update command line status
toc; tic; fprintf('Writing rho0_sgz... ');

% save
writeMatrix(filename, rho0, 'rho0_sgz');

% clear
clear('rho0');

% :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::

if flags.nonlinear

    % update command line status
    toc; tic; fprintf('Writing BonA... ');

    % set BonA
    if flags.BonA_heterog
        BonA = BonA_scalar*ones(Nx, Ny, Nz, 'single');
        BonA(Nx/2:end, :, :) = BonA_scalar*heterog_sc;
    else
        BonA = BonA_scalar;
    end

    % cast matrix to single precision
    eval(['BonA = ' MATRIX_DATA_TYPE_MATLAB '(BonA);']);

    % save
    writeMatrix(filename, BonA, 'BonA');

    % clear
    clear('BonA');
    
    % set flag
    nonlinear_flag = 1;

else
    
    % set flag
    nonlinear_flag = 0;
    
end

% :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::

if flags.absorbing

    % update command line status
    toc; tic; fprintf('Writing alpha_coeff... ');

    % set alpha_coeff
    if flags.alpha0_heterog
        alpha_coeff = alpha_coeff_dB*ones(Nx, Ny, Nz, 'single');
        alpha_coeff(Nx/2:end, :, :) = alpha_coeff_dB*heterog_sc;    
    else
        alpha_coeff = alpha_coeff_dB;
    end

    % cast matrix to single precision
    eval(['alpha_coeff = ' MATRIX_DATA_TYPE_MATLAB '(alpha_coeff);']);

    % save
    writeMatrix(filename, alpha_coeff, 'alpha_coeff');

    % clear
    clear('alpha_coeff');
    

% :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::

    % cast matrix to single precision
    eval(['alpha_power = ' MATRIX_DATA_TYPE_MATLAB '(alpha_power);']);
    
    % save
    writeMatrix(filename, alpha_power, 'alpha_power');
    
    % clear
    clear('alpha_power');
    
    % set flag
    absorbing_flag = 1;
    
else
    
    % set flag
    absorbing_flag = 0;

end

% :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::

% cast matrix to single precision
eval(['c_ref = ' MATRIX_DATA_TYPE_MATLAB '(c_ref);']);

% save
writeMatrix(filename, c_ref, 'c_ref');

% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
% ========================================================================
% 3. Set and save source flags
% ========================================================================
% \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

% maximum supported frequency
f_max = 1500/(2*dx_ref);

% create the input signal
input_signal = toneBurst(1/dt_ref, f_max/5, 3);

% set source flags (do not change these without actually changing the
% source matrices further down)
ux_source_flag = 0;
uy_source_flag = 0;
uz_source_flag = 0;
p0_source_flag = 0;
p_source_flag = length(input_signal);          % set to the length of the input file
transducer_source_flag = 0;
nonuniform_grid_flag = 0;

% :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::

% update command line status
toc; tic; fprintf('Writing source flags... ');

% integer variables
variable_names = {'ux_source_flag', 'uy_source_flag', 'uz_source_flag', ...
    'p_source_flag', 'p0_source_flag', 'transducer_source_flag',...
    'nonuniform_grid_flag', 'nonlinear_flag', 'absorbing_flag'};

% change all the index variables to be in 64-bit unsigned integers (long in C++)
for index = 1:length(variable_names)

    % cast matrix to 64-bit unsigned integer
    eval([variable_names{index} ' = ' INTEGER_DATA_TYPE_MATLAB '(' variable_names{index} ');']);

    % write to HDF5 file
    writeMatrix(filename, eval(variable_names{index}), variable_names{index});

end

% clear variables
clear(variable_names{:});

% clear holder
clear('variables_names');

% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
% ========================================================================
% 4. Set and save source variables
% ========================================================================
% \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

% update command line status
toc; tic; fprintf('Writing p source... ');

% set source mode to additive (1 is Additive, 0 is Dirichlet)
p_source_mode = 1;

% set source to use a single time series for all points
p_source_many = 0;

% set source mask to be a rectangle
source_offset = 20;
source_width = 40;
source_length = 80;

% set source indices to be a rectangle
p_source_index = sub2ind([Nx, Ny, Nz], (source_offset + 1)*ones(source_width*source_length, 1), ...
    repmat(Ny/2 - source_width/2 + 1:Ny/2 + source_width/2, 1, source_length).',...
    reshape(repmat(Nz/2 - source_length/2 + 1:Nz/2 + source_length/2, source_width, 1), [], 1));

% set source input
p_source_input = input_signal;

% clean up
clear source_offset source_width source_length input_signal f_max;

% :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::

% integer variables
variable_names = {'p_source_mode', 'p_source_many', 'p_source_index'};

% change all the index variables to be in 64-bit unsigned integers (long in C++)
for index = 1:length(variable_names)

    % cast matrix to 64-bit unsigned integer
    eval([variable_names{index} ' = ' INTEGER_DATA_TYPE_MATLAB '(' variable_names{index} ');']);

    % write to HDF5 file
    writeMatrix(filename, eval(variable_names{index}), variable_names{index});

end

% clear variables
clear(variable_names{:});

% clear holder
clear('variables_names');

% :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::

% cast matrix to single precision
eval(['p_source_input = ' MATRIX_DATA_TYPE_MATLAB '(p_source_input);']);

% save
writeMatrix(filename, p_source_input, 'p_source_input');

% clear
clear('p_source_input');
    
% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
% ========================================================================
% 5. Set and save sensor variables
% ========================================================================
% \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

% update command line status
toc; tic; fprintf('Writing sensor_mask_index... ');

% set the sensor mask indices to be central plane
start = (Nz/2 - 1)*Nx*Ny;
sensor_mask_index = (start + 1:start+Nx*Ny).';

% cast matrix to 64-bit unsigned integer
eval(['sensor_mask_index = ' INTEGER_DATA_TYPE_MATLAB '(sensor_mask_index);']);

% save
writeMatrix(filename, sensor_mask_index, 'sensor_mask_index');

% clear
clear('sensor_mask_index');

% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
% ========================================================================
% 7. Set and save k-space and shift variables
% ========================================================================
% \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

% update command line status
toc; tic; fprintf('Writing k-space and shift variables... ');

% create the wavenumber grids (assuming Nx, Ny and Nz are even)
nx = ((-Nx/2:Nx/2-1)/Nx).';
nx(floor(Nx/2) + 1) = 0;
kx_vec = (2*pi/dx).*nx; 

ny = ((-Ny/2:Ny/2-1)/Ny).';
ny(floor(Ny/2) + 1) = 0;
ky_vec = (2*pi/dy).*ny; 

nz = ((-Nz/2:Nz/2-1)/Nz).';
nz(floor(Nz/2) + 1) = 0;
kz_vec = (2*pi/dz).*nz; 

% force the vector operators be in the correct direction (Nx, 1, 1), (1, Ny, 1), (1, 1, Nz) 
ky_vec = ky_vec.'; 
kz_vec = permute(kz_vec, [2 3 1]);

% :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::

% create vector shift variables
ddx_k_shift_pos = ifftshift( 1i*kx_vec .* exp(1i*kx_vec*dx/2), 1);
ddx_k_shift_neg = ifftshift( 1i*kx_vec .* exp(-1i*kx_vec*dx/2), 1);
ddy_k_shift_pos = ifftshift( 1i*ky_vec .* exp(1i*ky_vec*dy/2), 2);
ddy_k_shift_neg = ifftshift( 1i*ky_vec .* exp(-1i*ky_vec*dy/2), 2);
ddz_k_shift_pos = ifftshift( 1i*kz_vec .* exp(1i*kz_vec*dz/2), 3);
ddz_k_shift_neg = ifftshift( 1i*kz_vec .* exp(-1i*kz_vec*dz/2), 3);
    
% create reduced variables for use with real-to-complex FFT (saves memory)
Nx_r = floor(Nx/2) + 1;
ddx_k_shift_pos_r = ddx_k_shift_pos(1:Nx_r);
ddx_k_shift_neg_r = ddx_k_shift_neg(1:Nx_r);

% float variables
variable_names = {...
    'ddx_k_shift_pos_r', 'ddx_k_shift_neg_r',...
    'ddy_k_shift_pos', 'ddy_k_shift_neg', ...
    'ddz_k_shift_pos', 'ddz_k_shift_neg'};

% change all the variables to be in single precision (float in C++), then
% add to HDF5 File
for index = 1:length(variable_names)

    % cast matrix to single precision
    eval([variable_names{index} ' = ' MATRIX_DATA_TYPE_MATLAB '(' variable_names{index} ');']);

    % write to HDF5 file
    writeMatrix(filename, eval(variable_names{index}), variable_names{index});

end

% clear saved files
clear(variable_names{:});
clear('ddx_k_shift_pos', 'ddx_k_shift_neg');

% clear additional variables
clear('kx_vec', 'ky_vec', 'kz_vec', 'nx', 'ny', 'nz');

% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
% ========================================================================
% 8. Set and save PML properties
% ========================================================================
% \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

% update command line status
toc; tic; fprintf('Writing pml variables... ');

% PML sizes
pml_x_size = 20;
pml_y_size = 20;
pml_z_size = 20;

% additional PML properties
pml_x_alpha = 2;
pml_y_alpha = 2;
pml_z_alpha = 2;

% vector PML variables
% execin(fullFunctionName, varargin)
pml_x = execin([getkWavePath '/private/getPML.m'], Nx, dx, dt, c_ref, pml_x_size, pml_x_alpha, false, 1);
pml_x_sgx = execin([getkWavePath '/private/getPML.m'], Nx, dx, dt, c_ref, pml_x_size, pml_x_alpha, true, 1);
pml_y = execin([getkWavePath '/private/getPML.m'], Ny, dy, dt, c_ref, pml_y_size, pml_y_alpha, false, 2);
pml_y_sgy = execin([getkWavePath '/private/getPML.m'], Ny, dy, dt, c_ref, pml_y_size, pml_y_alpha, true, 2);
pml_z = execin([getkWavePath '/private/getPML.m'], Nz, dz, dt, c_ref, pml_z_size, pml_z_alpha, false, 3);
pml_z_sgz = execin([getkWavePath '/private/getPML.m'], Nz, dz, dt, c_ref, pml_z_size, pml_z_alpha, true, 3);

% :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::

% integer variables
variable_names = {...
    'pml_x_size' , 'pml_y_size' , 'pml_z_size'};

% change all the index variables to be in 64-bit unsigned integers (long in C++)
for index = 1:length(variable_names)

    % cast matrix to 64-bit unsigned integer
    eval([variable_names{index} ' = ' INTEGER_DATA_TYPE_MATLAB '(' variable_names{index} ');']);

    % write to HDF5 file
    writeMatrix(filename, eval(variable_names{index}), variable_names{index});

end

% clear saved files
clear(variable_names{:});

% clear additional variables
clear c_ref

% :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::

% float variables
variable_names = {...
    'pml_x_sgx', 'pml_y_sgy', 'pml_z_sgz', ...
    'pml_x', 'pml_y', 'pml_z', ...
    'pml_x_alpha', 'pml_y_alpha', 'pml_z_alpha'};

% change all the variables to be in single precision (float in C++), then
% add to HDF5 File
for index = 1:length(variable_names)

    % cast matrix to single precision
    eval([variable_names{index} ' = ' MATRIX_DATA_TYPE_MATLAB '(' variable_names{index} ');']);

    % write to HDF5 file
    writeMatrix(filename, eval(variable_names{index}), variable_names{index});

end

% clear saved files
clear(variable_names{:});

% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
% ========================================================================
% File attributes
% ========================================================================
% \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

% update command line status
toc; tic; fprintf('Writing attributes... ');

% set file description
file_description = 'example simulation data using k-Wave';
kwave_ver = execin([getkWavePath '/private/getkWaveVersion.m']);

% set additional file attributes
h5writeatt(filename, '/', FILE_MAJOR_VER_ATT_NAME, HDF_FILE_MAJOR_VERSION);
h5writeatt(filename, '/', FILE_MINOR_VER_ATT_NAME, HDF_FILE_MINOR_VERSION);
h5writeatt(filename, '/', CREATED_BY_ATT_NAME, ['Save in parts using k-Wave ' kwave_ver]);
h5writeatt(filename, '/', FILE_DESCR_ATT_NAME, file_description);
h5writeatt(filename, '/', FILE_TYPE_ATT_NAME, HDF_INPUT_FILE);
h5writeatt(filename, '/', FILE_CREATION_DATE_ATT_NAME,  getDateString);

toc;

% update command line status
disp('-------------------------------------------------------------------');
disp(['Total computation time ' scaleTime(etime(clock, start_time))]);
disp('-------------------------------------------------------------------');

% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
% ========================================================================
% SPARE PARTS
% ========================================================================
% \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/


% % :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::
% 
% % convert the absorption coefficient to nepers.(rad/s)^-y.m^-1
% alpha_coeff = db2neper(alpha_coeff_dB, alpha_power);
% 
% % set absorb_tau
% absorb_tau = -2*alpha_coeff.*c0_scalar.^(alpha_power - 1)*ones(Nx, Ny, Nz);
% 
% % cast matrix to single precision
% eval(['absorb_tau = ' MATRIX_DATA_TYPE_MATLAB '(absorb_tau);']);
% 
% % save
% writeMatrix(filename, absorb_tau, 'absorb_tau');
% 
% % clear
% clear(absorb_tau);
% 
% % :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::
% 
% % set absorb_eta
% absorb_eta = 2*alpha_coeff.*c0_scalar.^(alpha_power)*tan(pi*alpha_power/2)*ones(Nx, Ny, Nz);
% 
% % cast matrix to single precision
% eval(['absorb_eta = ' MATRIX_DATA_TYPE_MATLAB '(absorb_eta);']);
% 
% % save
% writeMatrix(filename, absorb_eta, 'absorb_eta');
% 
% % clear
% clear(absorb_eta);
% 
% % :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::


% % :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::
% 
% % define plaid grids of the wavenumber components centered about 0
% k = zeros(Nx, Ny, Nz);
% k = bsxfun(@plus, (kx_vec.^2), k);
% k = bsxfun(@plus, (ky_vec.^2), k);
% k = bsxfun(@plus, (kz_vec.^2), k);
% k = sqrt(k);
% 
% % :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::
% 
% % create the k-space operator
% c_ref = c0_scalar;
% kappa_r = ifftshift( sinc(c_ref*dt*k/2) );
% kappa_r = kappa_r(1:Nx_r, :, :);
% 
% % cast matrix to single precision
% eval(['kappa_r = ' MATRIX_DATA_TYPE_MATLAB '(kappa_r);']);
% 
% % save
% writeMatrix(filename, kappa_r, 'kappa_r');
% 
% % clear
% clear(kappa_r);
% 
% % :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::
% 
% % compute the absorbing fractional Laplacian operator and coefficient
% absorb_nabla1_r = (k).^(alpha_power-2); 
% absorb_nabla1_r(isinf(absorb_nabla1_r)) = 0;
% absorb_nabla1_r = ifftshift(absorb_nabla1_r);
% absorb_nabla1_r = absorb_nabla1_r(1:Nx_r, :, :); 
% 
% % cast matrix to single precision
% eval(['absorb_nabla1_r = ' MATRIX_DATA_TYPE_MATLAB '(absorb_nabla1_r);']);
% 
% % save
% writeMatrix(filename, absorb_nabla1_r, 'absorb_nabla1_r');
% 
% % clear
% clear(absorb_nabla1_r);
% 
% % :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::
% 
% % compute the dispersive fractional Laplacian operator and coefficient
% absorb_nabla2_r = (kgrid.k).^(alpha_power-1); 
% absorb_nabla2_r(isinf(absorb_nabla2_r)) = 0;
% absorb_nabla2_r = ifftshift(absorb_nabla2_r);            
% absorb_nabla2_r = absorb_nabla2_r(1:Nx_r, :, :);
% 
% % cast matrix to single precision
% eval(['absorb_nabla2_r = ' MATRIX_DATA_TYPE_MATLAB '(absorb_nabla2_r);']);
% 
% % save
% writeMatrix(filename,absorb_nabla2_r, 'absorb_nabla2_r');
% 
% % clear
% clear(absorb_nabla2_r);
% 
% % :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::