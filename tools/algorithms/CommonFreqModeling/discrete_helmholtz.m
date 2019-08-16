%-------------------------------------------------------------------------------
% Parameters choosing for discrete Helmholtz operator.
%
% To understand better the output of this function, it is important to
% distinguish the computational grid from the physical grid.
%

% This function deals with all the stability issues and outputs the dimension
% and function handlers for the interpolators and extrapolators needed to
% convert a model from (and to) the physical grid.
%
% It may as well optionally output the discrete operator itself (enabled by
% default).
%
% USE:
%   [H,comp_grid,dH_forw,T_adj,ddH] = discrete_helmholtz(v,model,freq,opts,flog)
%
% INPUT:
%    v             - gridded velocity in (m/s) or gridded slowness squared (s^2/m^2) on the COMPUTATIONAL DOMAIN
%    model.{o,d,t} - o,d,n parameters of the model domain
%
%    model.unit    - string set either to "m/s" (meters per second) or "s2/m2" (slowness squared)
%               
%    freq          - frequency in Hertz (can be complex)
%    opts          - Helmholtz parameter struct (details below)
%
% OUTPUT:
%    H             - SPOT operator of the Helmholtz matrix in band-storage format
%    comp_grid     - struct of parameters corresponding to the computational grid, pml, etc., details below
%
%
% Structure of "opts":
% ---------------------------------------
% It may contain the following fields:
%  scheme       - one of
%                  - PDEopts.HELM_OPERTO27 for the 27 points staggered grid (default)
%                  - PDEopts.HELM_STD7 for the 7 points discretization. 
%  mat_free     - if true, use the matrix-free stencil-based codes (default: true)
%  nlam         - number of points per wavelength to be used. Default is 6 for
%                  PDEopts.HELM_OPERTO27
%  dt           - Can be provided as a single scalar or 3 scalars. This is the
%                 grid spacing to be used in the computational grid. 
%  pml          - Can be provided as a single scalar or 3 scalars. Number of points
%                 to be used in the PML layer. 
%  cut_pml      - if true, returning from the PML domain to the computational domain will be
%                 performed by truncation, if false, by the adjoint of extension (default: true)
%  pml_max      - Maximum number of points in the PML layer. This should be a scalar
%                 and it overrides all the previous choices concerning PML. (default: inf).
%  free_surface - If true, the number of PML points at the top of the boundary == 0 (default: false)
%  
%
%  ------------------ for solving the Helmholtz equation  ------------
%  mode         - specifies what operations are needed from the Helmholtz operator in order to store
%                 the least amount of information possible, one of
%
%                  opBandStorage.FORW_MULT_ONLY     - forward multiplication only
%                  opBandStorage.FORW_DIV_ONLY      - forward division only
%                  opBandStorage.MULT_ONLY          - forward/adjoint multiplication only
%                  opBandStorage.DIV_ONLY           - forward/adjoint division only
%                  opBandStorage.MULT_DIV           - multiplication and division
%  solver       - one of 
%                  - LinSolveOpts.CGMN     - conjugate gradient method
%                  - LinSolveOpts.CRMN     - conjugate residual method (default)
%                  - LinSolveOpts.FGMRES   - flexible gmres 
%  maxit        - maximum number of iterations for each PDE solve
%  precond      - one of
%                  - LinSolveOpts.PREC_KACZSWP    - uses Kaczmarz sweeps (default)
%                  - LinSolveOpts.PREC_MLCR       - Multigrid CR preconditioner
%                  - LinSolveOpts.PREC_IDENTITY   - identity, not recommended
%                  - LinSolveOpts.PREC_MLGMRES    - multilevel GMRES
%                  - LinSolveOpts.PREC_SHIFTLAP   - shifted laplacian preconditioner
%                  - LinSolveOpts.PREC_CALANDRA12 - multi-level shifted laplacian preconditioner of Calandra '12
% ---------------------------------------------------------------------
% 
%  
%
% Structure of "comp_grid" (output):
% ---------------------------
%   .{ot,dt,nt}      - o,d,n parameters corresponding to the computational grid
%   .pml_info        - struct detailing pml parameters
%      .x,y             - # PML points in the x, y directions
%      .top             - # PML points at the top of the model (0 if there is a free surface)
%      .bottom          - # PML points at the bottom of the model
%   .phys2comp       - function that interpolates a vector from the physical domain to the computational
%                      domain (w/ PML region)
%   .comp2phys       - function that interpolates a vector from the computational domain (w/ PML region)
%                      to the physical domain
%   .npts_wavelength - number of points per wavelength 
%
% AUTHOR: Curt Da Silva (2015)
%         Rafael Lago (2014)
%         
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: September, 2014
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%
%-------------------------------------------------------------------------------
function [H,comp_grid,T,DT_adj] = discrete_helmholtz(v,model,freq,opts)


% Parse parameters
%------------------------------
mat_free     = check_field(opts,'mat_free',true);
if(~isfield(opts,'solve_opts') || ~isa(opts.solve_opts,'LinSolveOpts') )        
        error('Need to specify solve_opts of type LinSolveOpts');
end
precond      = check_field(opts.solve_opts,'precond',LinSolveOpts.PREC_IDENTITY);

cut_pml      = check_field(opts,'cut_pml',true);
if ~isfield(opts,'scheme')
    error('Need helmholtz scheme');
end
scheme       = opts.scheme;
free_surface = check_field(opts,'free_surface',false);
pml_max      = check_field(opts,'pml_max',inf);
disp_output  = check_field(opts,'disp_output',false);
n_threads    = check_field(opts,'n_threads',5);
wn_offset  = check_field(opts,'wn_offset',0);

comp_grid    = struct;

if numel(v)~=prod(model.n)
    assert(min(size(v)) > 1 && (length(size(v))==2||length(size(v))==3),'v must be a 2D or 3D array or have prod(model.n) elements');
else
    v = reshape(v,model.n);
end

% For 2D default values for the third coordinate
if length(model.n)==2
    model.n(3) = 1; model.o(3) = 0; model.d(3) = 1;    
    opts.pml = model.nb;        
end
mode_3d =  model.n(3)>1;
if ~mode_3d, mat_free = false; end
d_ones = ones(1,3);

if isfield(opts,'dt') && ~isempty(opts.dt)
    if length(opts.dt) == 1
        dt = d_ones*opts.dt;
    else
        dt = opts.dt;              
    end
    if ~mode_3d, dt(3) = 1; end
else
    if numel(v)==prod(model.n)
        dt = model.d;
    else
        error('Specify dt');
    end
end

switch model.unit 
  case 'm/s', vmax = max(v(:)); 
  case 's2/m2',vmax = max(vec(v.^(-1/2)));
  case 's2/km2',vmax = max(vec(1e3*v.^(-1/2))); 
  otherwise, error(['Unrecognized unit ' model.unit]);
end

%------------------------------------------------------------------------------
% PML computation
%------------------------------------------------------------------------------
if isfield(opts,'pml') && ~isempty(opts.pml)
    if length(opts.pml) == 1
        npml = opts.pml*ones(1,3);
    else
        npml = opts.pml;
        if (mode_3d && length(npml)~=3) || ...
            (~mode_3d && (length(npml)~=2 && npml(3)>0))
            error('npml must have length 2 or 3');
        end
    end    
else
    npml = ceil(vmax./(real(freq)*dt)); % Maximum wavelength of the model
end
npml = min(npml,pml_max);
if length(vec(npml))==3
    npml = repmat(npml,2,1);
end
% npml is a 2x3 matrix corresponding to
% in 2D
% [ #pml pts lower z  | # pml pts lower x  | 0 ]
% [ #pml pts upper z  | # pml pts upper x  | 0 ]
% 
% in 3D
% [ #pml pts lower x  | # pml pts lower y  | # pml pts lower z ]
% [ #pml pts upper x  | # pml pts upper y  | # pml pts upper z ]
% 
if free_surface
    if mode_3d
        npml(1,3) = 0; 
    else
        npml(1,1) = 0;
    end
end
if ~mode_3d, npml(:,3) = 0; end

%------------------------------------------------------------------------------
% Define computational grid
%------------------------------------------------------------------------------
nt   = size(v); if ~mode_3d, nt(3) = 1; end
nt_nopml   = nt;
nt         = nt  + sum(npml,1);
ot         = model.o - npml(1,:) .* dt;

% Create padding and extension operators
%-------------------------------------------------------------------------------
if mode_3d
    Pext = opKron(opExtension(nt_nopml(3), [npml(1,3) npml(2,3)]),  ...
                  opExtension(nt_nopml(2), [npml(1,2) npml(2,2)]), ...
                  opExtension(nt_nopml(1), [npml(1,1) npml(2,1)]));
    Ppad = opKron(opExtension(nt_nopml(3), [npml(1,3) npml(2,3)],0),  ...
                  opExtension(nt_nopml(2), [npml(1,2) npml(2,2)],0), ...
                  opExtension(nt_nopml(1), [npml(1,1) npml(2,1)],0));
else
    Pext = opKron(opExtension(nt_nopml(2), [npml(1,2) npml(2,2)]), ...
                  opExtension(nt_nopml(1), [npml(1,1) npml(2,1)]));
    Ppad = opKron(opExtension(nt_nopml(2), [npml(1,2) npml(2,2)],0), ...
                  opExtension(nt_nopml(1), [npml(1,1) npml(2,1)],0));
end

comp_grid.ext           = Pext;
comp_grid.pad           = Ppad;
comp_grid.phys2comp     = Pext;
if cut_pml
    comp_grid.comp2phys = Ppad';
else
    comp_grid.comp2phys = Pext';
end

%------------------------------------------------------------------------------
% Obtain Coefficients 
%------------------------------------------------------------------------------
v_pml = reshape(comp_grid.phys2comp*vec(v),nt);

if mat_free         
    mat_mode = 'implicit';
    wn = param2wavenum(v_pml,freq,model.unit);    
    switch scheme
      case PDEopts.HELM3D_STD7
        Hmvp = FuncObj(@helm3d_mvp_7p_mex,{vec(wn),vec(dt),vec(nt),npml,freq,[],[]});
        jacobi = FuncObj(@helm3d_7pt_jacobi_mex,{vec(wn),vec(dt),vec(nt),npml,freq,[],[],[],[]});
        kacz_sweep = [];
      case PDEopts.HELM3D_OPERTO27              
        Hmvp = FuncObj(@helm3d_operto27_mvp,{wn,dt,nt,npml,[],n_threads,[],false});        
        %Hmvp = FuncObj(@Helm3dmvp,{wn,dt,nt,npml,[],[]});
        jacobi = [];
        kacz_sweep = FuncObj(@helm3d_operto27_kaczswp,{wn,dt,nt,npml,[],[],[],[],[]});        
        if nargout >= 3
            [~,wn] = param2wavenum(v_pml,freq,model.unit);
            dHmvp = FuncObj(@helm3d_operto27_mvp,{wn,dt,nt,npml,[],n_threads,[],true});
            [~,~,wn] = param2wavenum(v_pml,freq,model.unit);
            ddHmvp = FuncObj(@helm3d_operto27_mvp,{wn,dt,nt,npml,[],n_threads,[],true});
        end
      case PDEopts.HELM3D_CHEN2012
        Hmvp = FuncObj(@helm3d_mvp_chen2012_mex,{vec(v_pml),freq,dt,nt,npml,[],[]});
        jacobi = [];
        kacz_sweep = [];
      case PDEopts.HELM3D_CHEN27
        Hmvp = FuncObj(@helm3d_27pt_chen_mvp_mex,{vec(v_pml),freq,dt,nt,npml,[],[]});
        jacobi = [];
        kacz_sweep = [];
    end
    helm_params = struct;
    helm_params.multiply = Hmvp;
    helm_params.jacobi = jacobi;
    helm_params.kacz_sweep = kacz_sweep;
    helm_params.N = prod(nt);
    helm_params.iscomplex = true;
    solve_opts = copy(opts.solve_opts);
    if disp_output && strcmp(solve_opts.solver,LinSolveOpts.SOLVE_FGMRES) 
        addprop(solve_opts,'output_freq');
        solve_opts.output_freq = 1;
    end
    H = opAbstractMatrix(mat_mode,helm_params,solve_opts);
    if nargout >=3
        helm_params = rmfield(helm_params,'jacobi');
        helm_params = rmfield(helm_params,'kacz_sweep');
        helm_params.multiply = dHmvp;
        dH = opAbstractMatrix(mat_mode,helm_params,[]);
        T_forw = @(dm,u) dH*opDiag_swp(dm)*u;
        T_adj = @(z,u) conj(u) .* (dH'*z);
        T = @(u) opFunction_swp(size(u,1),size(u,1),{@(dm)T_forw(dm,u),@(z)T_adj(z,u)});        
        helm_params.multiply = ddHmvp;
        ddH = opAbstractMatrix(mat_mode,helm_params,[]);
        DT_adj_func = @(u,dm,du,z) conj(u) .* (opDiag_swp(dm)*ddH'*z) + conj(du).*(dH'*z);
        DT_adj = @(u,dm,du) opFunction_swp(size(u,1),size(u,1),{@(z) DT_adj_func(u,dm,du,z),@(z) z});
    end
else
    mat_mode = 'explicit';
    switch scheme
      case PDEopts.HELM3D_STD7
        [H,dH,ddH] = helmholtz_3d_7p(v_pml,dt,npml,freq,model.unit);
        helm_params = struct;
        helm_params.coef = H;
        solve_opts = copy(opts.solve_opts);
        if disp_output && strcmp(solve_opts.solver,LinSolveOpts.SOLVE_FGMRES)
            addprop(solve_opts,'output_freq');
            solve_opts.output_freq = 1;
        end
        H = opAbstractMatrix(mat_mode,helm_params,solve_opts);        
        dH_forw = @(dm,u) dH*opDiag_swp(dm)*u;
        T_adj = @(z,u) conj(u) .* (dH'*z);
      case PDEopts.HELM3D_OPERTO27
        [coef,idx] = helmholtz_3d(v_pml, dt,npml,freq,model.unit); 
        H = opBandStorage(coef,idx,opts);    
      case PDEopts.HELM2D_CHEN9P           
        [H,dH,ddH] = Helm2D_opt(v_pml,dt,nt,npml(:,1:2),model.unit,freq,model.f0,wn_offset);  
        helm_params = struct;
        helm_params.coef = H;
        H = opHelmholtz('explicit',helm_params,opts.solve_opts,false,disp_output);        
        
        T_forw = @(dm,u) dH*opDiag_swp(dm)*u;
        T_adj = @(z,u) conj(u) .* (dH'*z);
        T = @(u) opFunction_swp(size(u,1),size(u,1),{@(dm)T_forw(dm,u),@(z)T_adj(z,u)});
        DT_adj_func = @(u,dm,du,z) conj(u) .* (opDiag_swp(dm)*opMatrix(ddH')*z) + conj(du).*(opMatrix(dH')*z);
        DT_adj = @(u,dm,du) opFunction_swp(size(u,1),size(u,1),{@(z) DT_adj_func(u,dm,du,z),@(z) z});
      case PDEopts.POISS2D_FV
        [H,dH_forw,dH_adj,DdH_adj] = Poiss2D(v_pml,model.n);
        T = @(u) opFunction_swp(size(u,1),size(u,1),{@(de) dH_forw(de,u),@(z) dH_adj(z,u)});
        DT_adj = @(u,dm,du) opFunction_swp(size(u,1),size(u,1),{@(z) DdH_adj(u,du,dm,z),@(z) z});
    end
end

%------------------------------------------------------------------------------
% Return all the relevant information of the computational grid
%------------------------------------------------------------------------------
comp_grid.n  = model.n;
comp_grid.d  = model.d;
comp_grid.o  = model.o;
comp_grid.dt = dt;
comp_grid.nt_nopml = nt_nopml;
comp_grid.nt = nt;
comp_grid.ot = ot;
comp_grid.scheme = scheme;
comp_grid.pml = npml;

% Compute + assign the preconditioner to the spot operator
if isa(precond,'char')
switch precond
  case LinSolveOpts.PREC_MLCR
    precond = ML_CRMN(H,v,comp_grid,model,freq,opts);
  case LinSolveOpts.PREC_MLGMRES
    precond = ML_GMRES(H,v,comp_grid,model,freq,opts);    
  case LinSolveOpts.PREC_SHIFTLAP
    precond = SLap_multigrid(v,comp_grid,model,freq,opts);
  case LinSolveOpts.PREC_CALANDRA12
    precond = ML_CALANDRA12(H,v,comp_grid,model,freq,opts);
  case LinSolveOpts.PREC_VGMRES
    precond = V_GMRES(H,v,comp_grid,model,freq,opts);
  otherwise 
    precond = [];
end
end
if ~isempty(precond)
    solve_opts = H.solve_opts;
    solve_opts.precond = precond;
    H.solve_opts = solve_opts;
end

end
