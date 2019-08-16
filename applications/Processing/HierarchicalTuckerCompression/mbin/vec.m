function [ x ] = vec( x )
%VEC - vectorizes an ND-array
    
%x = x(:);
    x = reshape(x,numel(x),1);
end