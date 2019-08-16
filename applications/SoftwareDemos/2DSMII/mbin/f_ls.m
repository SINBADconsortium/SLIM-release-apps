function [val,grad] = f_ls(m,Q,dt,model)
% least-squares misfit for waveform inversion
%
% use:
%   [val,grad] = f_ls(m,Q,dt,model)
%
% input:
%   m  - model (see F.m)
%   Q  - source (see F.m)
%   dt - data (consistent with output of F.m)
%   model - struct with model parameters (see F.m)
%
% output:
%   val  - misfit
%   grad - gradient

% predicted data 
[d,DF] = F(m,Q,model); 

% residual 
r = d - dt;

% penalty
[val,df] = w_ls(r); 

% gradient
if nargout>1
    grad = DF'*df; 
end