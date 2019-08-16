function [f2c,c2f,n_sub] = fine2coarse(n,d,d_sub,type)
% FINE2COARSE - fine grid to coarse grid interpolation. 
%
%   type - one of
%          'cubic' - cubic lagrange interpolation
%          'linear' - linear interpolation
%
% Author: Curt Da Silva, 2015
%
% Usage:
%   [f2c,c2f,n_sub] = fine2coarse(n,d,d_sub,type);   
% 
% Input:
%   n       - fine grid model size (1, 2, or 3 vector, depending on the number of dimensions)
%   d       - fine grid spacing, same length as n
%   d_sub   - coarse grid spacing, must be elementwise greater than d
%  
% Output:
%   f2c     - fine to coarse grid interpolation with cubic spline interpolation
%   c2f     - coarse to fine grid interpolation with cubic spline interpolation
%   n_sub   - coarse grid model size
%
%
%
% Usage:
%   [f2c,c2f] = fine2coarse(n,n_sub,type);
%   
% Input:
%   n        - fine grid model size
%   n_sub    - coarse grid model size
% 
% Output:
%   [f2c,c2f] - same as above
    
    if nargin==3
        if(isa(d_sub,'char'))
            type = d_sub;
            n_sub = d;
        else
            assert( length(n)==1 || length(n)==2 || length(n)==3,'n must have 1, 2, or 3 elements');
            assert( length(d)==length(n), 'd must be the same size as n');
            assert( max( d./d_sub ) <= 1, 'd_sub must be greater than d, elementwise');                        
            n_sub = ceil(n .* d ./d_sub);
        end
    elseif nargin==2
        n_sub = d;        
    else
        n_sub = ceil(n.* d./ d_sub);
    end
    assert( max(n ./n_sub) >= 1, 'n_sub must be greater than n, elementwise');
    
    if length(n)==3 && n(3)==1
        ndims = 2;
    else
        ndims = length(n);
    end
    
    if exist('type','var')==0
        type = 'cubic';
    end
    switch type
      case 'cubic'
        interp_basis = @(xin,xout) opLInterp1D(xin,xout);
      case 'linear'
        interp_basis = @(xin,xout) LinInterp1D(xin,xout);
      otherwise 
        error('Unrecognized type');
    end
        
    
    switch ndims
      case 1
        f2c = interp_basis(linspace(0,1,n),linspace(0,1,n_sub));
        c2f = interp_basis(linspace(0,1,n_sub),linspace(0,1,n));        
      case 2
        f2c = opKron(interp_basis(linspace(0,1,n(2)),linspace(0,1,n_sub(2))),...
                     interp_basis(linspace(0,1,n(1)),linspace(0,1,n_sub(1))));
        c2f = opKron(interp_basis(linspace(0,1,n_sub(2)),linspace(0,1,n(2))),...
                     interp_basis(linspace(0,1,n_sub(1)),linspace(0,1,n(1))));
      case 3        
        c2f = opKron(interp_basis(linspace(0,1,n_sub(3)),linspace(0,1,n(3))), ...
                         interp_basis(linspace(0,1,n_sub(2)),linspace(0,1,n(2))),...
                         interp_basis(linspace(0,1,n_sub(1)),linspace(0,1,n(1))));
        if strcmp(type,'linear')
            f2c = c2f';
        else
            f2c = opKron(interp_basis(linspace(0,1,n(3)),linspace(0,1,n_sub(3))),...
                     interp_basis(linspace(0,1,n(2)),linspace(0,1,n_sub(2))),...
                     interp_basis(linspace(0,1,n(1)),linspace(0,1,n_sub(1))));
        end
    end
        
end