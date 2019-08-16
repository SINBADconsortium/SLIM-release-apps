function checkFor_validArrayIndex(x, dimension_size)
% Checks whether x is a valid array index for sizes upto dimension_size
%
%    checkFor_validArerayIndex(x, dimension_size) will throw an error if x is
%    not a valid array index, and if x is greater than dimension_size

%   Copyright 2010, Tim Lin

validateattributes(x,{'numeric'},{'positive','integer','finite','real','nonempty','nonnan','<=',dimension_size})
