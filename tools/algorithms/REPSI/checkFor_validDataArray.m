function checkFor_validDataArray(x)
% Mainly used to check whether endianess is wrong
%
%   Copyright 2010, Tim Lin

validateattributes(x,{'numeric'},{'nonempty','nonnan','finite','real'})
