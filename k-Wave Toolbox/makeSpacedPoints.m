function data = makeSpacedPoints(X1, X2, N, spacing)
%MAKESPACEDPOINTS     Create vector of log or linear spaced points
%
% DESCRIPTION:
%       makeSpacedPoints generates a row vector of either logarithmically
%       or linearly spaced points between X1 and X2. When spacing is set to
%       'linear', the function is identical to the inbuilt  linspace
%       function. When spacing is set to 'log', the function is similar to
%       the inbuilt logspace function, except that X1 and X2 define the
%       start and end numbers, not decades. For logarithmically spaced
%       data, X1 must be > 0. If N < 2, X2 is returned. 
%
% USAGE:
%       data = makeSpacedData(X1, X2)
%       data = makeSpacedData(X1, X2, N)
%       data = makeSpacedData(X1, X2, N, spacing)
%
% INPUTS:
%       X1          - starting data value
%       X2          - ending data value (where X2 > X1)
%
% OPTIONAL INPUTS:
%       N           - number of points in the vector (default = 100)
%       spacing     - 'log' or 'linear' spaced values (default = 'linear')
%
% OUTPUTS:
%       data        - row vector of equally spaced points
%
% ABOUT:
%       author      - Bradley E. Treeby
%       date        - 14th July 2005
%       last update - 1st December 2012
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2013 Bradley Treeby and Ben Cox

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

% check for number of points input
if nargin < 3
    N = 100;
end

% check for spacing input
if nargin < 4
    spacing = 'linear';
end

% check if the end point is larger than the start point
if X2 <= X1
    error('X2 must be larger than X1');
end

% force N to be an integer
N = round(N);

if (N < 2)
    % return the end point if N < 2
    data = X2;
elseif (N == 2)
    % return the start and end points if N = 2
    data(1) = X1;
    data(2) = X2;
else

    % update X1 and X2 values for log spaced variables
    if (strcmp(spacing, 'log'))
        
        % check that X1 is greater than 0
        if X1 <= 0
            error('X1 must be > 0 for log spaced data');
        end
        
        X1 = log10(X1);
        X2 = log10(X2);
    end
    
    % create step variable and data range
    step = (X2 - X1)/(N - 1);
    data = X1:step:X2;
    
    % update log spaced data
    if (strcmp(spacing, 'log'))
        data = 10.^(data);
    end
end