function A = opInterp(type,xin,xout,varargin)
% opInterp - grid to grid interpolation
%
%   type - one of
%          'cubic' - cubic lagrange interpolation
%          'linear' - linear interpolation
%
% Author: Curt Da Silva, 2016
%
% Usage:
%   A = opInterp(type,xin,xout,{yin},{yout},{zin},{zout});   
% 
% Input:
%   xin, xout - x-input, x-output grids, respectively
%   yin, yout - (optional) y-input, y-output grids
%   zin, zout - (optional) z-input, z-output grids
%
% Output:
%   A    - interpolation operator
%
    switch length(varargin)
      case 0 
        ndims = 1;
      case 2
        ndims = 2;
        yin = varargin{1}; yout = varargin{2};
      case 4
        ndims = 3;
        yin = varargin{1}; yout = varargin{2}; zin = varargin{3}; zout = varargin{4};
      otherwise
        error('Must have 0, 2, or 4 additional grid parameter arguments');
    end    
    
    sinc_window = 4;
   
    
    if exist('type','var')==0
        type = 'cubic';
    end
    switch type
      case 'sinc'
        interp_basis = @(xin,xout) opSincInterp(xin,xout,sinc_window);
      case 'cubic'
        interp_basis = @(xin,xout) opLInterp1D(xin,xout);
      case 'linear'
        interp_basis = @(xin,xout) LinInterp1D(xin,xout);
      otherwise 
        error('Unrecognized type');
    end
            
    switch ndims
      case 1
        A = interp_basis(xin,xout);        
      case 2
        A = opKron(interp_basis(yin,yout),interp_basis(xin,xout));
      case 3
        A = opKron(interp_basis(zin,zout),interp_basis(yin,yout),interp_basis(xin,xout));
 
    end
        
end