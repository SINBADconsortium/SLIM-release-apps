function checkFor_validInteger(x)
% Checks whether x is a valid integer
%
%   Copyright 2010, Tim Lin

validateattributes(x,{'numeric'},{'nonnegative','nonempty','nonnan','finite','real','integer'})
