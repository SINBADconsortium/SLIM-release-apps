% [Calandra et al. 2012. Numerical Linear Algebra With Applications 20]
% "An improved two-grid preconditioner for the solution of three-dimensional 
% Helmholtz problems in heterogeneous media"
% 
% This is the geometric multigrid preconditioner explained in the 
% aforementioned paper, Algorithm 1. This routine merely sets it up and returns
% a pointer to the proper function, as it requires quite a lot of parameters 
% setting and it might be a bit too much to expect the user to set them all 
% manually; to be used by (flexible) GMRES,  BiCGStab or any non-Hermitian 
% Krylov solver.
% 
% IMPORTANT:
% Be aware that you should use the 7 points stencil (helmholtz_3d_7p) to use
% this preconditioner. The relaxation parameters rely on a rigorous Fourier
% analysis on the 7 points discretization and thus are very particular for this 
% problem. If you change that, the behaviour of the preconditioner might change
% drastically.
% 
% 
% EXAMPLE:  
%   par.precon = GMG_CGPV12_T(Af,idxf,Ac,idxc,ntf,ntc,flog)
%            u = FGMRES(Af,idxf,b,x0,par,flog);
% 
% INPUT:
%   {Af,idxf}  - matrix in band storage format (no normalization needed) to be 
%                used for the fine level.
%   {Ac,idxc}  - matrix in band storage format (no normalization needed) to be 
%                used for the coarse level.
%   ntf        - grid size for the fine level, including the pml layer (used for
%                the prolongator/retriction)
%   ntc        - grid size for the coarse level, including the pml layer (used 
%                for the prolongator/retriction)
%   flog       - file identifier for the (already open) log file. 
%                Let it unset or set it to 0 if you are targeting 
%                PERFORMANCE or simply do not want any message printed.
% 
% 
% OUTPUT:
% 
% 
% EXAMPLE:
%   par.precon = GMG_CGPV12_T(Af,idxf,Ac,idxc,ntf,ntc,flog)
%            u = FGMRES(Af,idxf,b,x0,par,flog);
%   will compute an approximate solution for Af*x = b using flexible GMRES
%   preconditioned by the multigrid V cycle introduced in the Algorithm 1
%   of Calandra et al.
% 
% AUTHOR: Rafael Lago
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: August, 2014
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
% 
%------------------------------------------------------------------------------
function precon = GMG_CGPV12_T(Af,idxf,Ac,idxc,ntf,ntc,flog)


% Jacobi parameters - the same for all levels
%---------------------------------------------
par_jacobi.maxit = 2;                   % nu
par_jacobi.omega = .8;                  % omega_h
par_jacobi.NNN = ceil(length(idxf)/2);  % Jacobi will use Af(par.NNN,:) as the 
                                        % diagonal of the matrix   
                        

% Smoother parameters
%---------------------
par_smoother.maxcy  =  1;  % varphi
par_smoother.maxit  =  2;  % ms
par_smoother.tol    = -1;  %
par_smoother.precon = @(b)(Jacobi_relax(Af,idxf,b,0*b,par_jacobi));

% Coarse parameters
%---------------------
par_coarse.maxcy    = 10; % varphi_c
par_coarse.maxit    = 10; % mc
par_coarse.tol      = -1;
par_coarse.precon   = @(b)(Jacobi_relax(Ac,idxc,b,0*b,par_jacobi));


% Prolongation and Restriction
%-----------------------------
lxf       = 1:ntf(1);
lyf       = 1:ntf(2);
lzf       = 1:ntf(3);
lxc       = 1:ntf(1)/ntc(1):ntf(1);
lyc       = 1:ntf(2)/ntc(2):ntf(2);
lzc       = 1:ntf(3)/ntc(3):ntf(3);
Lc2f      = opKron(opLInterp1D(lzc,lzf),opLInterp1D(lyc,lyf),opLInterp1D(lxc,lxf));
Lf2c      = opKron(opLInterp1D(lzf,lzc),opLInterp1D(lyf,lyc),opLInterp1D(lxf,lxc));
par.prolongation = @(u)(Lc2f*u);
par.restriction  = @(u)(Lf2c*u);

% Assemble everything
%-----------------------------
par.maxit = 1;
par.pre_smoother = @(b,x0)(FGMRES(Af,idxf,b,x0,par_smoother,0));
par.coarse       = @(b,x0)(FGMRES(Ac,idxc,b,x0,par_coarse,0));
par.pst_smoother = par.pre_smoother;

precon = @(b)(GMG_V(Af,idxf,b,0*b,par,flog));

end
