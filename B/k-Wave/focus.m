function signal_mat = focus(kgrid, input_signal, source_mask, focus_position, sound_speed)
% 21st February 2012
% 30th September 2012
% Bradley Treeby

% allow focus position to also be a focus distance on axis?

% check that kgrid.t_array is defined
if strcmp(kgrid.t_array, 'auto')
    error('kgrid.t_array must be defined');
end

% calculate the distance from every point in the source mask to the focus
% position
switch kgrid.dim
    case 1
        dist = abs(kgrid.x(source_mask == 1) - focus_position(1));
    case 2
        dist = sqrt( (kgrid.x(source_mask == 1) - focus_position(1)).^2 + ...
            (kgrid.y(source_mask == 1) - focus_position(2)).^2 );
    case 3
        dist = sqrt( (kgrid.x(source_mask == 1) - focus_position(1)).^2 + ...
            (kgrid.y(source_mask == 1) - focus_position(2)).^2 + ...
            (kgrid.z(source_mask == 1) - focus_position(3)).^2 );
end

% convert the distance to units of time points
dist = round(dist./(kgrid.dt*sound_speed));

% convert time points to relative delays
dist = -(dist - max(dist(:)));
max_delay = max(dist(:));

% create an input matrix
signal_mat = zeros(length(dist), length(input_signal) + max_delay);

% assign the input signal
for source_index = 1:length(dist)
    delay = dist(source_index);
    signal_mat(source_index, :) = [zeros(1, delay), input_signal, zeros(1, max_delay - delay)];
end



% function to generate beam forming delays for an abitrary 2D surface based
% on a focus distance in the third dimension assuming the focus is on axis

% sz = size(source_mask);
% 
% % calculate the distance from every point in the disk to the focus point in
% % grid points
% dist_mat = makePixelMap(sz(1), sz(2), 'Shift', [0 0]);
% dist_mat = sqrt((dist_mat*dx).^2 + focus_dist.^2);
% dist_mat = dist_mat./c0;   % convert to [s]
% dist_mat = round(dist_mat./dt);    % convert to time points
% dist_mat = dist_mat.*source_mask;   % remove non-piston points
% 
% % convert to delays
% dist_mat(source_mask ~= 0) = abs(dist_mat(source_mask ~= 0) - max(dist_mat(:)));
% max_delay = max(dist_mat(:));
% num_grid_points = sum(source_mask(:));
% 
% % create an input matrix
% signal_mat = zeros(num_grid_points, length(input_signal) + max_delay);
% 
% % extract indices from mask
% mask_index = find(source_mask ~= 0);
% 
% % assign the driving signal
% for gp_index = 1:num_grid_points 
%     delay = dist_mat(mask_index(gp_index));
%     signal_mat(gp_index, :) = [zeros(1, delay), input_signal, zeros(1, max_delay - delay)];
% end