function varargout = execin(fullFunctionName, varargin)
%EXECIN  Execute a function or script in different directory.
%   EXECIN(FUNCNAME) will execute function FUNCNAME, which is a string that
%   includes the full path of the function.
%
%   [Y1, Y2, ...] = EXECIN(FUNCNAME, X1, X2, ...) allows input and output
%   arguments that are normally allowed by the function FUNCNAME.
%
%   Example:
%     [s, out] = execin('C:\mywork dir\testfunction.m', x1, x2);

%   VERSIONS:
%     v1.0 - First version
%     v1.1 - added ability to run scripts
%     v1.2 - minor syntax fix
%
% Jiro Doke
% Copyright 2005-2010 The MathWorks, Inc.

% Copyright (c) 2005-2010, MathWorks
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the  nor the names 
%       of its contributors may be used to endorse or promote products derived 
%       from this software without specific prior written permission.
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.


[p, n] = fileparts(fullFunctionName);

% use SETAPPDATA instead of assigning to a variable in case the script
% clears all variables.
setappdata(0, 'curDir', pwd);
setappdata(0, 'fullFunctionName', fullFunctionName);

try
  cd(p);
  if nargin > 1 % if there are additional arguments, assume it is a function
    [varargout{1:nargout}] = feval(n, varargin{:});
  else
    if nargout % if there are output arguments, assume it is a function
      [varargout{1:nargout}] = feval(n);
    else % use EVAL in case it is a script
      eval(n)
    end
  end    
  cd(getappdata(0, 'curDir')); rmappdata(0, 'curDir'); rmappdata(0, 'fullFunctionName');
catch ME
  cd(getappdata(0, 'curDir')); rmappdata(0, 'curDir');
  fprintf('%s\n', ME.message);
  fullFunctionName = getappdata(0, 'fullFunctionName'); rmappdata(0, 'fullFunctionName');
  error('execin:ErrorInExecution', 'Could not execute ''%s''', fullFunctionName);
end