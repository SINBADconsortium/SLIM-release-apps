function y = helm3d_operto27_mvp(wn,h,nt,npml,x,n_threads,mode,deriv_mode)
% Simple interface to mex files for multiplying a 3D Helmholtz matrix with a vector
% 
% Curt Da Silva, 2016
% 
% Usage:
%   y = helm3d_operto27_mvp(wn,h,nt,npml,x,n_threads,mode,deriv_mode);
%
% Input:
%   wn         - wavenumber squared (rad^2 s^2/m^2) of size nt
%   h          - grid spacing (m) in [x,y,z]
%   nt         - number of grid points in [x,y,z]
%   npml       - 2 x 3 matrix of 
%   x          - input vector
%   n_threads  - number of threads to use
%   deriv_mode - if true, perform derivative of Helmholtz times vector
%                if false, perform Helmholtz times vector
% 
% Output:
%   y          - Helmholtz matrix (or its derivative) times vector
    
    if mode==1
        if deriv_mode
            y = Helm3dmvp_forw_deriv_mex(wn,h,nt,npml,x,n_threads);
        else
            if isreal(wn)
                y = Helm3dmvp_forw_mex(wn,h,nt,npml,x,n_threads);
            else
                y = Helm3dmvp_forw_wnc_mex(wn,h,nt,npml,x,n_threads);
            end
        end
    else
        if deriv_mode
            y = Helm3dmvp_adj_deriv_mex(wn,h,nt,npml,x,n_threads);            
        else
            if isreal(wn)
                y = Helm3dmvp_adj_mex(wn,h,nt,npml,x,n_threads);
            else
                y = Helm3dmvp_adj_wnc_mex(wn,h,nt,npml,x,n_threads);
            end
        end
    end
   
        
end
