% [Calandra et al. 2012. Numerical Linear Algebra With Applications 20]
% "An improved two-grid preconditioner for the solution of three-dimensional 
% Helmholtz problems in heterogeneous media"
% 
% This is the geometric multigrid preconditioner explained in the 
% aforementioned paper, Algorithm 3. This routine merely sets it up and returns
% a pointer to the proper function, as it requires quite a lot of parameters 
% setting and it might be a bit too much to expect the user to set them all 
% manually; to be used by (flexible) GMRES,  BiCGStab or any non-Hermitian 
% Krylov solver.
% 
% The preconditioner is composed of a three grid V cycle, but it requires 4
% operators - two Helmholtz and two shifted Laplacian:
% 
%                           Ω_h     Ω_h
%    (V cycle on             \      /     
%     Helmholtz operator)     \    /
%                              \  /          
%                              Ω_2h  ==>  Ω_2h   Ω_2h       
%                                          \     /     
%                                           \   /      (V cycle on shifted 
%                                            \ /        Laplacian operator) 
%                                           Ω_4h       
%
% IMPORTANT:
% Be aware that you should use the 7 points stencil (helmholtz_3d_7p) to use
% this preconditioner with 10 points per wavelength. The relaxation parameters 
% rely on a rigorous Fourier analysis on this particular discretization. If you
% change that, the behaviour of the preconditioner may change drastically.
% 
% EXAMPLE:  
%   par.precon = GMG_CGPV12_T2V(helm_fine,model)
%            u = FGMRES(helm_fine.coef,helm_fine.idx,b,x0,par);
% 
% INPUT:
%   helm_fine  - output of discrete_helmholtz routine. It is expected that
%                you used discretize=true when you called discrete_helmholtz,
%                and that both helm_fine.coef and helm_fine.idx are present.
%                This is the operator used by the finest level Ω_h
%   model      - the same structure you passed to discrete_helmholtz routine
%                in order to obtain a discrete operator.
% 
% AUTHOR: Rafael Lago
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: November, 2014
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
% 
%------------------------------------------------------------------------------
function precon = GMG_CGPV12_T2V(helm_fine,model,adj,flog)

% Parse parameters
%------------------------------
if ~exist('flog','var')
   flog = 0;
end

beta1 =  1;
beta2 = .5*sqrt(-1);

optsm.d      = 2*helm_fine.d;
optsm.pml    = ceil(helm_fine.pml.x/2);
optsm.scheme = helm_fine.scheme;
optsc.d      = 4*helm_fine.d;
optsc.pml    = ceil(helm_fine.pml.x/4);
optsc.scheme = helm_fine.scheme;

%----------------------------------------------------------------------------------------------------
% Get Matrices
%----------------------------------------------------------------------------------------------------
ntf  = helm_fine.nt;

% Discretize level 2 Helmholtz
%--------------------------------------
plog(flog,'* Computing level 2 Helmholtz operator...  ');
tic; Hm = discrete_helmholtz(model,helm_fine.f,optsm,flog); T = toc;
plog(flog,' done in ', T , ' seconds \n');
if adj
  Hm   = conj(Htransp(Hm.coef,Hm.idx));
else
  Hm   = Hm.coef;
end

%--------------------------------------
% Discretize level 2 shifted Laplacian
%--------------------------------------
plog(flog,'* Computing level 2 shifted Laplacian operator...  ');
tic; SLm = discrete_helmholtz(model,sqrt(beta1+beta2)*helm_fine.f,optsm,flog); T = toc;
plog(flog,' done in ', T , ' seconds \n');
ntm  = SLm.nt;
if adj
  idxm = -SLm.idx;
  SLm  = conj(Htransp(SLm.coef,SLm.idx));
else
  idxm = SLm.idx;
  SLm  = SLm.coef;
end

%------------------------------------------------
% Discretize level 1 shifted Laplacian (coarsest)
%-------------------------------------------------
plog(flog,'* Computing coarsest shifted Laplacian operator 27points stencil...  ');
tic; SLc = discrete_helmholtz(model,sqrt(beta1+beta2)*helm_fine.f,optsc,flog); T = toc;
plog(flog,' done in ', T , ' seconds \n');
ntc  = SLc.nt;
if adj
  idxc = -SLc.idx;
  SLc  = conj(Htransp(SLc.coef,SLc.idx));
else
  idxc = SLc.idx;
  SLc  = SLc.coef;
end


%----------------------------------------------------------------------------------------------------
% Build up preconditioner
%----------------------------------------------------------------------------------------------------
% Level 3 - finer level      Ω_h
% Level 2 - in between (lol) Ω_2h
% Level 1 - coarsest level   Ω_4h
%--------------------------------------

% Prolongation and Restriction for level 3
%-------------------------------------------
lxf       = 1:ntf(1);
lyf       = 1:ntf(2);
lzf       = 1:ntf(3);
lxm       = 1:(ntf(1)-1)/(ntm(1)-1):ntf(1);
lym       = 1:(ntf(2)-1)/(ntm(2)-1):ntf(2);
lzm       = 1:(ntf(3)-1)/(ntm(3)-1):ntf(3);
Lm2f      = opKron(opLInterp1D(lzm,lzf),opLInterp1D(lym,lyf),opLInterp1D(lxm,lxf));
Lf2m      = opKron(opLInterp1D(lzf,lzm),opLInterp1D(lyf,lym),opLInterp1D(lxf,lxm));
par_f.prolongation = @(u)(Lm2f*u);
par_f.restriction  = @(u)(Lf2m*u); % More efficient than using pseudoinverse

% Prolongation and Restriction for level 2
%-------------------------------------------
lxm       = 1:ntm(1);
lym       = 1:ntm(2);
lzm       = 1:ntm(3);
lxc       = 1:(ntm(1)-1)/(ntc(1)-1):ntm(1);
lyc       = 1:(ntm(2)-1)/(ntc(2)-1):ntm(2);
lzc       = 1:(ntm(3)-1)/(ntc(3)-1):ntm(3);
Lc2m      = opKron(opLInterp1D(lzc,lzm),opLInterp1D(lyc,lym),opLInterp1D(lxc,lxm));
Lm2c      = opKron(opLInterp1D(lzm,lzc),opLInterp1D(lym,lyc),opLInterp1D(lxm,lxc));
par_m.prolongation = @(u)(Lc2m*u);
par_m.restriction  = @(u)(Lm2c*u); % More efficient than using pseudoinverse

%-----------------------------------------
% 
% Build up level 2
%
%-----------------------------------------

% Jacobi parameters
%---------------------------------------------
par_jacobi_m.maxit =  2;                  % nu_beta
par_jacobi_m.omega = .8;                  % omega_l, with l = 2
par_jacobi_m.NNN = ceil(length(idxm)/2);  % Jacobi will use Af(par.NNN,:) as the
                                          % diagonal of the matrix 

par_jacobi_c.maxit =  2;                  % nu_beta
par_jacobi_c.omega = .2;                  % omega_l, with l = 1
par_jacobi_c.NNN = ceil(length(idxc)/2);  


% Coarse parameters
%-----------------------------------------
par_coarse_m.maxcy    =  1; % varphi_beta
par_coarse_m.maxit    = 10; % m_beta
par_coarse_m.tol      = -1;
par_coarse_m.precon   = @(b)(Jacobi_relax(SLc,idxc,b,0*b,par_jacobi_c));

% Smoother parameters
%---------------------
par_smoother_m.maxcy  =  1;  % varphi
par_smoother_m.maxit  =  2;  % m_beta
par_smoother_m.tol    = -1;  %
par_smoother_m.precon = @(b)(Jacobi_relax(SLm,idxm,b,0*b,par_jacobi_m));

% Assemble everything - Level 2
%-----------------------------
par_m.maxit = 1;
par_m.pre_smoother = @(b,x0)(Jacobi_relax(SLm,idxm,b,x0,par_jacobi_m));
par_m.coarse       = @(b,x0)(FGMRES(SLc,idxc,b,x0,par_coarse_m));
par_m.pst_smoother = par_m.pre_smoother;
level2_precon      = @(b)(GMG_V(SLm,idxm,b,0*b,par_m));

%-------------------------------------------------------------------------------
% 
% Build up level 3
%
%-------------------------------------------------------------------------------

% Jacobi parameters
%---------------------------------------------
par_jacobi_f.maxit = 2;                   % nu
par_jacobi_f.omega = .8;                  % omega_h
par_jacobi_f.NNN = ceil(length(helm_fine.idx)/2);   
                                          

% Smoother parameters
%---------------------
par_smoother_f.maxcy  =  1;  % varphi
par_smoother_f.maxit  =  2;  % ms
par_smoother_f.tol    = -1;  %
par_smoother_f.precon = @(b)(Jacobi_relax(helm_fine.coef,helm_fine.idx,b,0*b,par_jacobi_f));

% Coarse parameters
%-----------------------------------------
par_coarse_f.maxcy    =  2; % varphi_c
par_coarse_f.maxit    = 10; % mc
par_coarse_f.tol      = -1;
par_coarse_f.precon   = level2_precon;


% Assemble everything - Level 3
%-----------------------------
par_f.maxit = 1;
par_f.pre_smoother = @(b,x0)(FGMRES(helm_fine.coef,helm_fine.idx,b,x0,par_smoother_f));
par_f.coarse       = @(b,x0)(FGMRES(Hm,idxm,b,x0,par_coarse_f));
par_f.pst_smoother = par_f.pre_smoother;

precon = @(b)(GMG_V(helm_fine.coef,helm_fine.idx,b,0*b,par_f));

end
