% Script to generate HDF5 test data
%
% author: Bradley Treeby
% date: 12th September 2012
% last update: 12th September 2012

clear all;

% add toolbox
addpath('../../k-Wave Toolbox');

% folder to put output data
output_folder = './';

% generate required data (copy and paste as many times as required,
% changing the simulation case and grid size)

for simulation_case = 1:5 
    grid_size = [128, 128, 128];
    file_name = [output_folder 'input_data_' num2str(grid_size(1)) '_' num2str(grid_size(2)) '_' num2str(grid_size(3)) '_case_' num2str(simulation_case) '.h5'];
    saveHDF5InParts(simulation_case, grid_size, file_name)
end 

