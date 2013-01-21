% DESCRIPTION:
%       subscript to get the indices of the active elements in sensor.mask
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 6th September 2012
%       last update - 6th September 2012
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
% along with k-Wave. If not, see <http://www.gnu.org/license

% create mask indices (this works for both normal sensor and transducer
% inputs)
sensor_mask_index = find(sensor.mask ~= 0);

% convert the data type depending on the number of indices (this saves
% memory)
eval(['sensor_mask_index = ' index_data_type '(sensor_mask_index);']); 