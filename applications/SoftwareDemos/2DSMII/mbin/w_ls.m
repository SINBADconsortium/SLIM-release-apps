function [f,g] = w_ls(r)
% least-squares penalty and gradient
%
% use:
%   [val,grad] = w_ls(r)
%
% input:
%   r  - vector
%
% output:
%   val  - two-norm of residual (.5*norm(r).^2)
%   grad - gradient

% two-norm of residual
f = .5*norm(r)^2; 

% gradient
g = r;
