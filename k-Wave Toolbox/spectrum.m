function varargout = spectrum(varargin)
%SPECTRUM   Compute the single sided amplitude and phase spectrums.
%
% DESCRIPTION:
%       The k-Wave spectrum function has been renamed to spect to avoid
%       naming conflicts with the signal processing toolbox.
%
%       This is a wrapper function to ensure backwards compatability. In
%       future releases, this wrapper will be removed. 
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 26th September 2012
%       last update - 28th September 2012
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

% give warning
disp('WARNING: The k-Wave spectrum function has been renamed to spect to avoid a naming conflict with the Signal Processing Toolbox. Please update usage.');

% call spect
[f, func_as, func_ps] = spect(varargin{:});

% assign the outputs
if nargout == 1
    varargout(1) = {func_as};
elseif nargout == 2;
    varargout(1) = {f};
    varargout(2) = {func_as};
elseif nargout == 3;
    varargout(1) = {f};
    varargout(2) = {func_as};
    varargout(3) = {func_ps};    
end