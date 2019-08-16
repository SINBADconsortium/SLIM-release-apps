function checkFor_validFraction(x)
% Checks whether x is a valid fraction of 1
%
%   Copyright 2010, Tim Lin

validateattributes(x,{'numeric'},{'nonnegative','nonempty','nonnan','finite','real','>=',0,'<=',1})
