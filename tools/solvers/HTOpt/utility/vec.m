function [ x ] = vec( x )
%VEC - vectorizes an ND-array
%
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca

    
%x = x(:);
    x = reshape(x,numel(x),1);
end

