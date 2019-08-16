% Returns a function handler to the appropriate solve (both forward and 
% adjoint).
% 
% This is *NOT* suitable for blocks of right hand sides (yet); you have to loop
% the call for the solve over the number of right hand sides. Future versions 
% may support that.
% 
% This also lacks the use of parallelism.
%
% USE:
%   [H, A] = helmholtz_solver(model,freq,par_helm,par_solver,flog)
%
% INPUT:
%  model.v    - velocity model in meters per second or slowness squared; this
%               is a *3D object* with dimension (nx,ny,nz), NOT A VECTOR!
%  model.dv   - array with the distance between each point in each direction
%               for this velocity model. This is used to compute the physical
%               size of the model - which should not change under any 
%               circunstance.
%  model.unit - string set either to "m/s" (meters per second) or "s2/m2" 
%               (slowness squared)
%  freq       - frequency in Hertz
%  par_helm   - [OPTIONAL] allows the user to enforce some parameters.
%               Read further for more details.
%  par_solver - [OPTIONAL] Contains parameters and options to be passed to the 
%               chosen solve. Read further for more details.
%
% OUTPUT:
%           H - A structure containing at least H.coef (the discrete operator
%               in bandstorage format, e.g. the output of helmholtz_3d) and
%               H.idx, the respective coefficients for the matrix in bandstorage
%               format and a pointer to the solver itself. When used as
%               u = H.solve(b,x0),
%               u will contain the approximation of the solution of the Helmholtz
%               equation, according to the parameters chosen during the
%               creation of this function.
%               Additionally, H contains all the fields from the output of the 
%               function discrete_helmholtz; please refer to "help discrete_helmholtz"
%               for a detailed documentation.
%               If the chosen solver is Matlab's backslash, this structure
%               will additionally contain H.sparse which is a sparse representation
%               of the Helmholtz operator.
%           A - Same as H, but for the adjoint of the wave equation. All the 
%               entries are duplicated (including A.coef) which requires 
%               considerable memory, but might be far more computationally
%               efficient than transposing the matrix for every PDE solve
%               (to be investigated! - Lago)
%               This is an optional output.
%
%                
% Structure of "par_helm" (input, optional):
% -------------------------------------------
% It may contain the following fields:
%  scheme  - set to either "operto27" for the 27 points staggered grid 
%            discretization or to "std7pts" for the 7 points discretization.
%            Only "operto27" is supported at the moment, as "std7pts" is still
%            experimental.
%  discretize - set to true or false. If set to true, the output structure H 
%            will also contain the discrete Helmholtz operator. Otherwise, it
%            will just contain the proper set of parameters for a stable
%            discretization. Default: true.
%  nlam    - number of points per wavelength to be used. Default is 6 for
%            "operto27".
%  d       - Can be provided as a single scalar or 3 scalars. This is the 
%            grid spacing to be used in the computational grid. It is higly
%            adviseable to never use this option and let discrete_helmholtz
%            to automatically compute it.
%  pml     - Can be provided as a single scalar or 3 scalars. Number of points
%            to be used in the PML layer. It is higly adviseable to never use 
%            this option and let discrete_helmholtz to automatically compute it.
%  pml_max - Maximum number of points in the PML layer. This should be a scalar
%            and it overrides all the previous choices concerning PML. Default
%            is "inf".
%  free_surface - set to true or false. If this is true, the number of points
%            in the PML layer is set to zero on the top of the domain only.
%            Default is false.
%  ext_first-First extending the domain (that is, add the PML layer) and then 
%           interpolate it to the proper dimension satisfying stability 
%           conditions. This affects the accuracy: ext_first = true is more 
%           accurate, but ext_first = false uses less memory. Default is 
%           "false".
%           
%           
% Structure of "par_solver" (input, optional):
% --------------------------------------- -----
% par_solver  - Contains parameters and options to be passed to the chosen 
%               solver.
%            In addition to any field present in par_solver, you can define
%            which solver itself to be used by setting the proper flag.
%            For instance
%               par_solver.name = 'carpcg';
%               par_solver.tol = 1e-4;
%               par_solver.maxit = 200;
%            will use at most 200 iterations of CARPCG to try to achieve
%            a tolerance of 1e-4 for the relative preconditioned residual. 
%            Conversely, the following
%               par_solver.tol = 1e-4;
%               par_solver.maxit = 200;
%            will use the maxit and tol fields only if an iterative solver 
%            is chosen by "helmholtz_solver", otherwise those parameters will
%            be ignored.
% 
% AUTHOR: Rafael Lago
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: January, 2014
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
% 
% TODO: Allow a different solver for forward and for the adjoint solve - Lago
% TODO: Allow ONLY adjoint - Lago 
%
%------------------------------------------------------------------------------- 
function [H, A] = helmholtz_solver(model,freq,par_helm,par_solver,flog)

% These are the available solvers
% (well, sort of; this needs to be properly implemented)
%-------------------------------------------------------------------------------

if ~exist('par_helm','var')
  par_helm = [];
end
if ~exist('par_solver','var')
  par_solver = [];
end
if ~exist('flog','var')
  flog = 0;
end

solve_name = check_field(par_solver,'name','');

use_carpcg = 0;
use_carpcr = 0;
use_backslash = 0;
use_mlcrmn = 0;
use_mlcgmn = 0;
use_cgpv12_t = 0;
use_cgpv12_t2v = 0;
use_auto = 0;
switch solve_name
  case 'carpcg'
    use_carpcg = 1;
  case 'carpcr'
    use_carpcr = 1;
  case {'backslash','\'}
    use_backslash = 1;
  case {'mlcrmn'}
    use_mlcrmn = 1;
  case {'mlcgmn'}
    use_mlcgmn = 1;
    error('Solver not implemented yet');
  case {'cgpv12_t'}
    use_cgpv12_t = 1;
    error('Solver not implemented yet');
  case {'cgpv12_t2v'}
    use_cgpv12_t2v = 1;
  otherwise
    use_auto = 1;
end


% We get the smallest stable grid size. This will be used, amongst other things,
% to determine which solve is to be used (e.g. if the grid size is small, 
% we might choose a direct solve).
% 
% This does not get the discrete operator yet, because if multigrid is to be
% used, things need to be rounded up.
stable = discrete_helmholtz(model,freq,...
                            setfield(par_helm,'discretize',false));

%-------------------------------------------------------------------------------
% AUTO
% This is what will select the solve automatically when the user doesn't - 
% based on the number of right hand sides, the size of the problem, memory 
% available, frequency, number of processors, all that jazz.
% 
% In the future. Some future.
% 
% TODO: Implement the future
%-------------------------------------------------------------------------------
if use_auto
  
  use_carpcr = 1; % Default!
  
  % For loose accuracy, forget fancy solvers and stick to the simple CRMN
  if isfield(par_solver,'tol') & par_solver.tol <1e-4 
    use_carpcr = 1;
  
  % For very large cases where 7 pts is being used, go for the economic
  % Multigrid T2V
  elseif sum(stable.nt) > 1800 && strcmp(stable.scheme,'std7pts')
    use_cgpv12_t2v = 1;
  
  % TODO: This is a very wild guess...
  elseif sum(stable.nt) > 100
    use_mlcrmn = 1;
  
  % Finally, if you don't know what to do, go for CRMN
  else
    use_carpcr = 1;
  end
end

%-------------------------------------------------------------------------------
% BACKSLASH, CARPCG and CARPCR
%-------------------------------------------------------------------------------
if use_carpcr || use_carpcg || use_backslash

  H = discrete_helmholtz(model,freq,...
                               setfield(par_helm,'discretize',true),flog);
  
  %----------------------------------------------------------------------------
  % BACKSLASH...
  % Does nothing with x0, uses Matlab's backslash for the rest.
  % 
  % FIXME:TODO: We should be able to STORE the approximate inverse generated by
  % backslash in order to solve for each right hand side separately rather
  % than for the whole block. 
  %----------------------------------------------------------------------------
  if use_backslash
    % TODO: Storing Sparse AND Bandstorage format! This is not good, fix this asap - Lago
    H.sparse = H2sparse(H.coef,H.idx); 
    H.solve = @(b,x0)(H.sparse\b(:));
    H.solver_name = 'backslash';
    
    warning(['MATLAB''s backslash has been chosen in setup_solve, but in its ' ...
            'current implementation, it is not optimal. Each PDE solved will  '...
            'recompute the approximate inverse for every right hand side. \n'...
            'This should be optimized in a future release, ask Bas Peters '...
            'and Curt da Silva about this issue.']);
    
    if nargout>1
      % TODO: Stores the sparse AND its adjoint! Might be too much, investigate this asap - Lago
      A = H;
      A.sparse = A.sparse';
      A.solve = @(b,x0)(A.sparse\b(:));
      A.solver_name = 'backslash';
    end
  end

  %----------------------------------------------------------------------------
  % CARPCG
  % We use the same setup for CARPCR
  %----------------------------------------------------------------------------
  if use_carpcg || use_carpcr
    Nt = prod(H.nt);
    HScale = 1./sqrt(sum(abs(H.coef).^2,1));
    HScale = spdiags(transpose(HScale),0,Nt,Nt);

    % FIXME: This line might force matlab to store H.coef twice - the existing
    % one and the H.coef*HScale one. Check this further - Lago
    par_solver.size = H.nt;
    H.solve = @(b,x0)(CGMN(H.coef*HScale,H.idx,HScale*b(:),...
                           x0,par_solver));
    H.solver_name = 'carpcg';

    if nargout>1
        % FIXME: This will force matlab to store H.coef twice - the existing
        % one and the transposed one. This is not necessarily a bad thing, 
        % considering that it might take too long to transpose it back and forth.
        % 
        % This is something that deserves further investigation.
        % - Lago
        A = H;
        A.coef = conj(Htransp(H.coef,H.idx));
        A.idx  = -H.idx;
        AScale = 1./sqrt(sum(abs(A.coef).^2,1));
        AScale = spdiags(transpose(AScale),0,Nt,Nt);
        A.solve = @(b,x0)(CGMN(A.coef*AScale,A.idx,AScale*b(:),...
                               x0,par_solver));
        A.solver_name = 'carpcg';
    end
  end

  %----------------------------------------------------------------------------
  % CARPCR
  %----------------------------------------------------------------------------
  if use_carpcr
        H.solve = @(b,x0)(CRMN(H.coef*HScale,H.idx,HScale*b(:),...
                                  x0,par_solver));
        H.solver_name = 'carpcr';
    if nargout>1
        A.solve = @(b,x0)(CRMN(A.coef*AScale,A.idx,AScale*b(:),...
                                  x0,par_solver));
        A.solver_name = 'carpcr';
    end
  end
end

%-------------------------------------------------------------------------------
% MLCRMN, MLCGMN and CGPV12_T2V
%-------------------------------------------------------------------------------
if use_mlcrmn || use_mlcgmn || use_cgpv12_t2v
  
  % Divides points by four, round it up, multiply by four;
  % This ensures that the coarser grid has a number of point which is
  % a multiple of four.
  par_helm.d   = ceil((model.nv-1).*model.dv)./(4*ceil((stable.n-1)/4));
  par_helm.discretize = 1;
  
  % Now we get the discrete operator
  H = discrete_helmholtz(model,freq,par_helm,flog);

  %----------------------------------------------------------------------------
  % ML-CRMN
  %----------------------------------------------------------------------------
  if use_mlcrmn

    par_solver.precon = ML_CRMN(H,model,0);
    H.solve = @(b,x0)(FGMRES(H.coef,H.idx,b(:),x0(:),par_solver));
    H.solver_name = 'mlcrmn';
    
    if nargout>1
      % FIXME: This will force matlab to store H.coef twice - the existing
      % one and the transposed one. This is not necessarily a bad thing, 
      % considering that it might take too long to transpose it back and forth.
      % 
      % This is something that deserves further investigation.
      % - Lago
      A = H;
      A.coef = conj(Htransp(H.coef,H.idx));
      A.idx  = -H.idx;
      par_solver.precon = ML_CRMN(A,model,1);
      A.solve = @(b,x0)(FGMRES(A.coef,A.idx,b(:),x0(:),par_solver));
      A.solver_name = 'mlcrmn';
    end
  end

  %----------------------------------------------------------------------------
  % GMG_CGPV12_T2V
  % Note that for this to work we need to have the optional argument "model" 
  % passed.
  %----------------------------------------------------------------------------
  if use_cgpv12_t2v
    if ~strcmp(H.scheme,'std7pts')
      warning(['GMG_CGPV12_T2V was either choosen or selected as solve by ' ...
              ' setup_solve, but the chosen discretization scheme is "' ...
              H.scheme '". This preconditioner is only guaranteed to converge'...
              'for std7pts.']);
    end
    
    par_solver.precon = GMG_CGPV12_T2V(H,model);
    H.solve = @(b,x0)(FGMRES(H.coef,H.idx,b,x0,par_solver));
    H.solver_name = 'cgpv12_t2v';
    
    if nargout>1
      % FIXME: This will force matlab to store H.coef twice - the existing
      % one and the transposed one. This is not necessarily a bad thing, 
      % considering that it might take too long to transpose it back and forth.
      % 
      % This is something that deserves further investigation.
      % - Lago
      A = H;
      A.coef = conj(Htransp(H.coef,H.idx));
      A.idx  = -H.idx;
      par_solver.precon = ML_CRMN(A,model,1);
      A.solve = @(b,x0)(FGMRES(A.coef,A.idx,b(:),x0(:),par_solver));
      A.solver_name = 'mlcrmn';
    end
  end
end

clear stable;

end
