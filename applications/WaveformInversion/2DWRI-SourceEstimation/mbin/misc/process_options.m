function [varargout] = process_options(args, varargin)
% PROCESS_OPTIONS - Sets default options for unset parameters.
%
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
% Usage:
%   [var_1, ..., var_n] = process_options(args, varname_1, varval_1, ..., varname_n,varval_n)
%
%   opts = process_options(args, varname_1, varval_1, ..., varname_n,varval_n);
%
% Input:
%   args                      - cell array of input arguments, alternating
%                               between variable names and values
%   varname_1, ..., varname_n - strings corresponding to variable names
%   varval_1, ..., varval_n   - default values for the variable
%                               names, if not specified
%
% Output:
%   var_1, ..., var_n  - values assigned to variables
%   opts               - struct containing variables + assigned values

n = length(varargin);
if (mod(n, 2) == 1)
  error('Each option must be a string/value pair.');
end

if (nargout > 1 && nargout < (n / 2))
  error('Insufficient number of output arguments given');
else
  numout = n / 2;
end
 

varargout = cell(1, numout);
nameout = cell(1,numout);
for i=2:2:n
  varargout{i/2} = varargin{i};
  nameout{i/2} = varargin{i-1};
end

for i=1:2:length(args)
  for j=1:2:n
    if strcmpi(args{i}, varargin{j})
      varargout{(j + 1)/2} = args{i + 1};      
      break;
    end
  end
  
end

if nargout == 1
    output = struct;
    for i=1:numout
        output.( nameout{i} ) = varargout{i};
    end
    varargout = cell(1,1);
    varargout{1} = output;
end

