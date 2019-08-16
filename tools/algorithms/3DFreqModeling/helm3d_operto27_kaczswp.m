function x = helm3d_operto27_kaczswp(wn,h,nt,npml,w,x,b,nsweeps,mode)
% Wrapper function for stencil-based kaczmarz sweeps
% 
% Curt Da Silva, 2016
% 
% Usage:
%   y = helm3d_operto27_kaczswp(wn,h,nt,npml,w,x,b,nsweeps,mode);
%
% Input:
%   wn      - wavenumber
%   h       - grid spacing, 3x1 vector
%   nt      - number of grid points, 3x1 vector
%   npml    - number of pml points, 2x3 matrix (see )
%   w       - overrelaxation parameter
%   x       - current solution estimate
%   b       - right hand side
%   nsweeps - number of sweeps
%   mode    - if mode==1, sweep on helmholtz matrix, otherwise its adjoint
%    
    for k=1:nsweeps            
        if mode==1
            if isreal(wn)
                x = Helm3d_kacz_sweep_mex(wn,h,nt,npml,w,x,b);
            else
                x = Helm3d_kacz_sweep_wnc_mex(wn,h,nt,npml,w,x,b);
            end
        else
            x = Helm3d_kacz_sweep_adj_mex(wn,h,nt,npml,w,x,b);
        end           
    end            

end