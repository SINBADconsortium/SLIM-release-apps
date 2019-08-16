function [x,r,g,info] = spgl1( A, b, tau, sigma, x, options )
%SPGL1  Solve basis pursuit, basis pursuit denoise, and LASSO
%
% [x, r, g, info] = spgl1(A, b, tau, sigma, x0, options)
%
% ---------------------------------------------------------------------
% Solve the basis pursuit denoise (BPDN) problem
%
% (BPDN)   minimize  ||x||_1  subj to  ||Ax-b||_2 <= sigma,
%
% or the l1-regularized least-squares problem
%
% (LASSO)  minimize  ||Ax-b||_2  subj to  ||x||_1 <= tau.
% ---------------------------------------------------------------------
%
% INPUTS
% ======
% A        is an m-by-n matrix, explicit or an operator.
%          If A is a function, then it must have the signature
%
%          y = A(x,mode)   if mode == 1 then y = A x  (y is m-by-1);
%                          if mode == 2 then y = A'x  (y is n-by-1).
%
% b        is an m-vector.
% tau      is a nonnegative scalar; see (LASSO).
% sigma    if sigma != inf or != [], then spgl1 will launch into a
%          root-finding mode to find the tau above that solves (BPDN).
%          In this case, it's STRONGLY recommended that tau = 0.
% x0       is an n-vector estimate of the solution (possibly all
%          zeros). If x0 = [], then SPGL1 determines the length n via
%          n = length( A'b ) and sets  x0 = zeros(n,1).
% options  is a structure of options from spgSetParms. Any unset options
%          are set to their default value; set options=[] to use all
%          default values.
%
% OUTPUTS
% =======
% x        is a solution of the problem
% r        is the residual, r = b - Ax
% g        is the gradient, g = -A'r
% info     is a structure with the following information:
%          .tau     final value of tau (see sigma above)
%          .rNorm   two-norm of the optimal residual
%          .rGap    relative duality gap (an optimality measure)
%          .gNorm   Lagrange multiplier of (LASSO)
%          .stat    = 1 found a BPDN solution
%                   = 2 found a BP sol'n; exit based on small gradient
%                   = 3 found a BP sol'n; exit based on small residual
%                   = 4 found a LASSO solution
%                   = 5 error: too many iterations
%                   = 6 error: linesearch failed
%                   = 7 error: found suboptimal BP solution
%                   = 8 error: too many matrix-vector products
%          .time    total solution time (seconds)
%          .nProdA  number of multiplications with A
%          .nProdAt number of multiplications with A'
%
% OPTIONS
% =======
% Use the options structure to control various aspects of the algorithm:
%
% options.fid         File ID to direct log output
%        .verbosity   0=quiet, 1=some output, 2=more output.
%        .iterations  Max. number of iterations (default if 10*m).
%        .bpTol       Tolerance for identifying a basis pursuit solution.
%        .optTol      Optimality tolerance (default is 1e-4).
%        .decTol      Larger decTol means more frequent Newton updates.
%        .subspaceMin 0=no subspace minimization, 1=subspace minimization.
%        .quitPareto  0=normal execution, 1=forces an exit when the pareto curve is reached
%        .minPareto   Minimum number of spgl1 iterations before checking for quitPareto
%
% EXAMPLE
% =======
%   m = 120; n = 512; k = 20; % m rows, n cols, k nonzeros.
%   p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
%   A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';
%   b  = A*x0 + 0.005 * randn(m,1);
%   opts = spgSetParms('optTol',1e-4);
%   [x,r,g,info] = spgl1(A, b, 0, 1e-3, [], opts); % Find BP sol'n.
%
% AUTHORS
% =======
%  Ewout van den Berg (ewout78@cs.ubc.ca)
%  Michael P. Friedlander (mpf@cs.ubc.ca)
%    Scientific Computing Laboratory (SCL)
%    University of British Columbia, Canada.
%
% BUGS
% ====
% Please send bug reports or comments to
%            Michael P. Friedlander (mpf@cs.ubc.ca)
%            Ewout van den Berg (ewout78@cs.ubc.ca)

% 15 Apr 07: First version derived from spg.m.
%            Michael P. Friedlander (mpf@cs.ubc.ca).
%            Ewout van den Berg (ewout78@cs.ubc.ca).
% 17 Apr 07: Added root-finding code.
% 18 Apr 07: sigma was being compared to 1/2 r'r, rather than
%            norm(r), as advertised.  Now immediately change sigma to
%            (1/2)sigma^2, and changed log output accordingly.
% 24 Apr 07: Added quadratic root-finding code as an option.
% 24 Apr 07: Exit conditions need to guard against small ||r||
%            (ie, a BP solution).  Added test1,test2,test3 below.
% 15 May 07: Trigger to update tau is now based on relative difference
%            in objective between consecutive iterations.
% 15 Jul 07: Added code to allow a limited number of line-search
%            errors. 
% 23 Feb 08: Fixed bug in one-norm projection using weights. Thanks
%            to Xiangrui Meng for reporting this bug.
% 26 May 08: The simple call spgl1(A,b) now solves (BPDN) with sigma=0.
%
% 18 Feb 10: Code branched to SPGL1-SLIM to include custom features:
%               -added option to force exit when Pareto curve is reached
%               -included ability to handle Kaczmarz operators.
%               -hacked out the linesearch fail conditions using a large limit. 
%            Tim Lin (timtylin@gmail.com).
% 03 May 10: Added capabilities to work on distributed vectors. Added options for
%            parallelized capabilities, max line-search iterations.
% 20 May 10: Made minPareto a user-configurable condition
% 20 Jul 10: Improved performance of spgLine for costly Aprod
% 17 Nov 10: Made default allowance for feasible line-search artificially large to allow for badly scaled problems
% 20 Dec 10: Added option to allow L1-projection failures to issue warning istead of error (ignorePErr)
% 27 Apr 11: Some optimizations via the following:
%               -L1 projection on distributed arrays no longer uses sorting
%               -disabled restoring to xBest to save memory space
%               -disabled subspace minimizaiton permanately to avoid computing active set to save space
%               -re-written some expressions to avoid calculating imtermediate results
% 03 May 11: Further memory optimizations, calculated many quantites in-place to avoid temporary allocation of memory, made oneProjector nested

%
%   ----------------------------------------------------------------------
%   This file is part of SPGL1 (Spectral Projected-Gradient for L1).
%
%   Copyright (C) 2007 Ewout van den Berg and Michael P. Friedlander,
%   Department of Computer Science, University of British Columbia, Canada.
%   All rights reserved. E-mail: <{ewout78,mpf}@cs.ubc.ca>.
%
%   SPGL1 is free software; you can redistribute it and/or modify it
%   under the terms of the GNU Lesser General Public License as
%   published by the Free Software Foundation; either version 2.1 of the
%   License, or (at your option) any later version.
%
%   SPGL1 is distributed in the hope that it will be useful, but WITHOUT
%   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
%   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General
%   Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with SPGL1; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
%   USA
%   ----------------------------------------------------------------------
REVISION = '$Rev: 46 $';
DATE     = '$Date: 2011-06-14 19:35:37 -0700 (Tue, 14 Jun 2011) $';
REVISION = REVISION(6:end-1);
DATE     = DATE(35:50);

% Set to true to display debug flags
PRINT_DEBUG_FLAGS = false;

% A check to verify that this is the SLIM version of SPGl1
if isstr(A) && strcmp(A,'is_SLIM_version')
    x = 1;
    return
end

tic;              % Start your watches!
m = length(b);

%----------------------------------------------------------------------
% Check arguments. 
%----------------------------------------------------------------------
if ~exist('options','var'), options = []; end
if ~exist('x','var'), x = []; end
if ~exist('sigma','var'), sigma = []; end
if ~exist('tau','var'), tau = []; end

if nargin < 2 || isempty(b) || isempty(A)
   error('At least two arguments are required');
elseif isempty(tau) && isempty(sigma)
   tau = 0;
   sigma = 0;
   singleTau = false;
elseif isempty(sigma) % && ~isempty(tau)  <-- implied
   singleTau = true;
else
   if isempty(tau)
      tau = 0;
   end
   singleTau = false;
end

%----------------------------------------------------------------------
% Grab input options and set defaults where needed. 
%----------------------------------------------------------------------
defaultopts = spgSetParms(...
'fid'        ,      1 , ... % File ID for output
'verbosity'  ,      2 , ... % Verbosity level
'iterations' ,   10*m , ... % Max number of iterations
'nPrevVals'  ,      3 , ... % Number previous func values for linesearch
'bpTol'      ,  1e-06 , ... % Tolerance for basis pursuit solution 
'optTol'     ,  1e-04 , ... % Optimality tolerance
'decTol'     ,  1e-04 , ... % Req'd rel. change in primal obj. for Newton
'stepMin'    ,  1e-16 , ... % Minimum spectral step
'stepMax'    ,  1e+05 , ... % Maximum spectral step
'rootMethod' ,      2 , ... % Root finding method: 2=quad,1=linear (not used).
'activeSetIt',    Inf , ... % Exit with EXIT_ACTIVE_SET if nnz same for # its.
'subspaceMin',      0 , ... % Use subspace minimization
'iscomplex'  ,    NaN , ... % Flag set to indicate complex problem
'maxMatvec'  ,    Inf , ... % Maximum matrix-vector multiplies allowed
'weights'    ,      1 , ... % Weights W in ||Wx||_1
'Kaczmarz'   ,      0 , ... % Toggles whether Kaczmarz mode is on (experimental)
'KaczScale'  ,      1 , ... % Scaling factor for Tau when using Kaczmarz-type submatrices
'quitPareto' ,      0 , ... % Exits when pareto curve is reached
'minPareto'  ,      3 , ... % If quitPareto is on, the minimum number of iterations before checking for quitPareto conditions
'lineSrchIt' ,      1 , ... % Maximum number of line search iterations for spgLineCurvy, originally 10 ...
'feasSrchIt' ,  10000 , ... % Maximum number of feasible direction line search iteraitons, originally 10 ...
'ignorePErr' ,      0 , ... % Ignores projections error by issuing a warning instead of an error ...
'project'    , @NormL1_project , ...
'primal_norm', @NormL1_primal  , ...
'dual_norm'  , @NormL1_dual      ...
   );
options = spgSetParms(defaultopts, options);

fid           = options.fid;
logLevel      = options.verbosity;
maxIts        = options.iterations;
nPrevVals     = options.nPrevVals;
bpTol         = options.bpTol;
optTol        = options.optTol;
decTol        = options.decTol;
stepMin       = options.stepMin;
stepMax       = options.stepMax;
activeSetIt   = options.activeSetIt;
subspaceMin   = options.subspaceMin;
maxMatvec     = max(3,options.maxMatvec);
weights       = options.weights;
quitPareto    = options.quitPareto;
minPareto     = options.minPareto;
Kaczmarz      = options.Kaczmarz;
lineSrchIt    = options.lineSrchIt;
feasSrchIt    = options.feasSrchIt;
ignorePErr    = options.ignorePErr;

% maxLineErrors TEMPORARILY DISABLED to prevent very large scaling issue to throw off algorithm
maxLineErrors = Inf;     % Maximum number of line-search failures.
pivTol        = 1e-12;  % Threshold for significant Newton step.

%----------------------------------------------------------------------
% Initialize local variables.
%----------------------------------------------------------------------
iter          = 0;  itnTotLSQR = 0; % Total SPGL1 and LSQR iterations.
nProdA        = 0;  nProdAt    = 0;
lastFv        = -inf(nPrevVals,1);  % Last m function values.
nLineTot      = 0;                  % Total no. of linesearch steps.
printTau      = false;
nNewton       = 0;
bNorm         = norm(b,2);
stat          = false;
timeProject   = 0;
timeMatProd   = 0;
nnzIter       = 0;                  % No. of its with fixed pattern.
nnzIdx        = [];                 % Active-set indicator.
subspace      = false;              % Flag if did subspace min in current itn.
stepG         = 1;                  % Step length for projected gradient.
testUpdateTau = 0;                  % Previous step did not update tau

% Determine initial x, vector length n, and see if problem is complex
explicit = ~(isa(A,'function_handle'));
if isempty(x)
   if isnumeric(A)
      n = size(A,2);
      realx = isreal(A) && isreal(b);
   else
      x = Aprod(b,2);
      n = length(x);
      realx = isreal(x) && isreal(b);
   end
   % NOTE: branch at the initialization below to deal with distributed input of b
   % x = zeros(n,1);
else
   n     = length(x);
   realx = isreal(x) && isreal(b);
end
if isnumeric(A), realx = realx && isreal(A); end;

% Override options when options.iscomplex flag is set
if (~isnan(options.iscomplex)), realx = (options.iscomplex == 0); end

% Check if all weights (if any) are strictly positive. In previous
% versions we also checked if the number of weights was equal to
% n. In the case of multiple measurement vectors, this no longer
% needs to apply, so the check was removed.
if ~isempty(weights)
  if any(~isfinite(weights))
     error('Entries in options.weights must be finite');
  end
  if any(weights <= 0)
     error('Entries in options.weights must be strictly positive');
  end
else
  weights = 1;
end

% Quick exit if sigma >= ||b||.  Set tau = 0 to short-circuit the loop.
if bNorm <= sigma
   printf('W: sigma >= ||b||.  Exact solution is x = 0.\n');
   tau = 0;  singleTau = true;
end 
  
% Don't do subspace minimization if x is complex.
if ~realx && subspaceMin
   printf('W: Subspace minimization disabled when variables are complex.\n');
   subspaceMin = false;
end

% Pre-allocate iteration info vectors
xNorm1 = zeros(min(maxIts,10000),1);
rNorm2 = zeros(min(maxIts,10000),1);
lambda = zeros(min(maxIts,10000),1);

% Exit conditions (constants).
EXIT_ROOT_FOUND    = 1;
EXIT_BPSOL1_FOUND  = 2;
EXIT_BPSOL2_FOUND  = 3;
EXIT_OPTIMAL       = 4;
EXIT_ITERATIONS    = 5;
EXIT_LINE_ERROR    = 6;
EXIT_SUBOPTIMAL_BP = 7;
EXIT_MATVEC_LIMIT  = 8;
EXIT_ACTIVE_SET    = 9; % [sic]
EXIT_AT_PARETO     = 10;

%----------------------------------------------------------------------
% Log header.
%----------------------------------------------------------------------
printf('\n');
printf(' %s\n',repmat('=',1,82));
printf(' SPGL1_SLIM v.%s (%s) based on v.1017\n', REVISION, DATE);
printf(' %s\n',repmat('=',1,82));
printf(' %-22s: %8i %4s'   ,'No. rows'          ,m       ,'');
printf(' %-22s: %8i\n'     ,'No. columns'       ,n          );
printf(' %-22s: %8.2e %4s' ,'Initial tau'       ,tau     ,'');
printf(' %-22s: %8.2e\n'   ,'Two-norm of b'     ,bNorm      );
printf(' %-22s: %8.2e %4s' ,'Optimality tol'    ,optTol  ,'');
if singleTau
   printf(' %-22s: %8.2e\n'  ,'Target one-norm of x'  ,tau       );
else
   printf(' %-22s: %8.2e\n','Target objective'  ,sigma      );
end
printf(' %-22s: %8.2e %4s' ,'Basis pursuit tol' ,bpTol   ,'');
printf(' %-22s: %8i\n'     ,'Maximum iterations',maxIts     );
printf('\n');
if singleTau
   logB = ' %5i  %13.7e  %13.7e  %9.2e  %6.1f  %6i  %6i';
   logH = ' %5s  %13s  %13s  %9s  %6s  %6s  %6s\n';
   printf(logH,'Iter','Objective','Relative Gap','gNorm','stepG','nnzX','nnzG');
else
   logB = ' %5i  %13.7e  %13.7e  %9.2e  %9.3e  %6.1f  %6i  %6i';
   logH = ' %5s  %13s  %13s  %9s  %9s  %6s  %6s  %6s  %13s\n';
   printf(logH,'Iter','Objective','Relative Gap','Rel Error',...
          'gNorm','stepG','nnzX','nnzG','tau');
end    
    
% Project the starting point and evaluate function and gradient.
if isempty(x)
    r         = b;  % r = b - Ax
    g         = -Aprod(r,2);  % g = -A'r
    f         = norm(r)^2 / 2;
    dx        = project(-g, tau);
else
    x         = project(x,tau);
    r         = b - Aprod(x,1);  % r = b - Ax
    g         =   - Aprod(r,2);  % g = -A'r
    f         = norm(r)^2 / 2;
    dx        = project(x - g, tau) - x;
end

% Compute initial steplength.
dxNorm = norm(dx,inf);
if dxNorm < (1 / stepMax)
   gStep = stepMax;
else
   gStep = min( stepMax, max(stepMin, 1/dxNorm) );
end

% Required for nonmonotone strategy.
lastFv(1) = f;
fBest     = f;
% xBest     = x;
fOld      = f;

dispFlag('fin Init')

clear dx

%----------------------------------------------------------------------
% MAIN LOOP.
%----------------------------------------------------------------------
while 1
    
    %------------------------------------------------------------------
    % Test exit conditions.
    %------------------------------------------------------------------

    % Compute quantities needed for log and exit conditions.
    gNorm   = undist(options.dual_norm(g,weights)); % originally options.dual_norm(-g,weights), but for true norms the sign should not matter
    rNorm   = norm(r, 2);
    gap     = dot(r,(r-b)) + tau*gNorm;
    rGap    = abs(gap) / max(1,f);
    aError1 = rNorm - sigma;
    aError2 = f - sigma^2 / 2;
    rError1 = abs(aError1) / max(1,rNorm);
    rError2 = abs(aError2) / max(1,f);
    
    dispFlag('fin CompConditions')
    
    if isempty(x)
        nnzX    = 0;
    else
        nnzX    = sum(abs(x) >= min(.1,10*options.optTol));
    end
    nnzG    = 0; % this is a temporary measure to make sure we don't destroy the log formatting
    
    % Single tau: Check if we're optimal.
    % The 2nd condition is there to guard against large tau.
    if singleTau
       if rGap <= optTol || rNorm < optTol*bNorm
          stat  = EXIT_OPTIMAL;
       end
 
    % Multiple tau: Check if found root and/or if tau needs updating.
    else
       
       if rGap <= max(optTol, rError2) || rError1 <= optTol
          % The problem is nearly optimal for the current tau.
          % Check optimality of the current root.
          test1 = rNorm       <=   bpTol * bNorm;
          test2 = gNorm       <=   bpTol * rNorm;
          test3 = rError1     <=  optTol;
          test4 = rNorm       <=  sigma;
          
          if test4, stat=EXIT_SUBOPTIMAL_BP;end % Found suboptimal BP sol.
          if test3, stat=EXIT_ROOT_FOUND;   end % Found approx root.
          if test2, stat=EXIT_BPSOL2_FOUND; end % Gradient zero -> BP sol.
          if test1, stat=EXIT_BPSOL1_FOUND; end % Resid minim'zd -> BP sol.
       end

       testRelChange1 = (abs(f - fOld) <= decTol * f);
       testRelChange2 = (abs(f - fOld) <= 1e-1 * f * (abs(rNorm - sigma)));
       testUpdateTau  = ((testRelChange1 && rNorm >  2 * sigma) || ...
                         (testRelChange2 && rNorm <= 2 * sigma)) && ...
                         ~stat && ~testUpdateTau;
       
       if testUpdateTau
          if quitPareto && iter >= minPareto, stat=EXIT_AT_PARETO;end % Chose to exit out of SPGL1 when pareto is reached 
          % Update tau.
          tauOld   = tau;
          tau      = max(0,tau + (rNorm * aError1) / (gNorm * options.KaczScale));
          nNewton  = nNewton + 1;
          printTau = abs(tauOld - tau) >= 1e-6 * tau; % For log only.
          if tau < tauOld
             % The one-norm ball has decreased.  Need to make sure that the
             % next iterate if feasible, which we do by projecting it.
             if ~isempty(x)
                 x = project(x,tau);
             end
          end
          % See if new rows need to be chosen for Kaczmarz
          if Kaczmarz
              A([],'new_rows');
          end
       end
    end

    % Too many its and not converged.
    if ~stat  &&  iter >= maxIts
        stat = EXIT_ITERATIONS;
    end
    dispFlag('fin CheckConverge')

    %------------------------------------------------------------------
    % Print log, update history and act on exit conditions.
    %------------------------------------------------------------------
    if logLevel >= 2 || singleTau || printTau || iter == 0 || stat
       tauFlag = '              '; subFlag = '';
       if printTau, tauFlag = sprintf(' %13.7e',tau);   end
       if subspace, subFlag = sprintf(' S %2i',itnLSQR); end
       if singleTau
          printf(logB,undist(iter),undist(rNorm),undist(rGap),undist(gNorm),log10(undist(stepG)),undist(nnzX),undist(nnzG));
          if subspace
             printf('  %s',subFlag);
          end
       else
          printf(logB,undist(iter),undist(rNorm),undist(rGap),undist(rError1),undist(gNorm),log10(undist(stepG)),undist(nnzX),undist(nnzG));
          if printTau || subspace
             printf(' %s',[tauFlag subFlag]);
          end
       end
       printf('\n');
    end
    printTau = false;
    subspace = false;
    
    % Update history info
    if isempty(x)
        xNorm1(iter+1) = 0;
    else
        xNorm1(iter+1) = options.primal_norm(x,weights);
    end
    rNorm2(iter+1) = rNorm;
    lambda(iter+1) = gNorm;
    
    if stat, break; end % Act on exit conditions.
        
    %==================================================================
    % Iterations begin here.
    %==================================================================
    iter = iter + 1;
    xOld = x;  fOld = f; rOld = r; % gOld update moved down to coincide with gradient update
    try
       %---------------------------------------------------------------
       % Projected gradient step and linesearch.
       %---------------------------------------------------------------
       dispFlag('begin LineSrch')
       [nLine,stepG,lnErr] = spgLineCurvy(max(lastFv));
       dispFlag('fin LineSrch')
       nLineTot = nLineTot + nLine;
       if lnErr
          %  Projected backtrack failed. Retry with feasible dir'n linesearch.
          dispFlag('begin FeasLineSrch')
          clear x
          clear r
          x    = xOld;
          
          % In-place scaling of gradient and updating of x to save memory
          if ~isempty(xOld)
              dx = project(xOld - gStep.*g, tau) - xOld;
          else
              dx = project(-gStep .* g, tau);
          end
          
          gtd  = dot(g,dx);
          
          [f,step,r,nLine,lnErr] = spgLine(f,dx,gtd,rOld,max(lastFv),@Aprod,b,feasSrchIt);
          dispFlag('fin FeasLineSrch')
          
          if isempty(xOld)
              x = step*dx;
          else
              x = xOld + step*dx;
          end
          clear dx
          
          x = project(x, tau);
          
          nLineTot = nLineTot + nLine;
       end
       if lnErr
       %  Failed again.  Revert to previous iterates and damp max BB step.
          if maxLineErrors <= 0
             stat = EXIT_LINE_ERROR;
          else
             stepMax = stepMax / 10;
             printf(['W: Linesearch failed with error %i. '...
                     'Damping max BB scaling to %6.1e.\n'],lnErr,stepMax);
             maxLineErrors = maxLineErrors - 1;
          end
       end
       
       primNorm_x = undist(options.primal_norm(x,weights));
       targetNorm = tau+optTol;
       
       if ignorePErr
            if primNorm_x > targetNorm
                warning('Primal norm of projected x is larger than expected, project again to be safe')
                primNorm_x
                targetNorm
            end
       else
            ensure(primNorm_x <= targetNorm);
       end
       
       dispFlag('fin UpdateX')
       
       %---------------------------------------------------------------
       % Update gradient and compute new Barzilai-Borwein scaling.
       %---------------------------------------------------------------
       if isempty(xOld)
           s    = x;
       else
           xOld    = x - xOld; % in-place calculating of s
           s       = xOld;
           clear xOld
       end
       
       gOld = g;
       g    = - Aprod(r,2);
       y    = g - gOld;
       clear gOld
       
       sts  = norm(s)^2;
       sty  = dot(s,y);
       if   sty <= 0,  gStep = stepMax;
       else            gStep = min( stepMax, max(stepMin, sts/sty) );
       end
       
       dispFlag('fin CompScaling')
       
       clear s
       clear y
       
    catch % Detect matrix-vector multiply limit error WARNING: the latest round of optimizations may have broke this, do testing to veriify
       err = lasterror;
       if strcmp(err.identifier,'SPGL1:MaximumMatvec')
         stat = EXIT_MATVEC_LIMIT;
         iter = iter - 1;
         x = xOld;  f = fOld;  r = rOld;
         break;
       else
         rethrow(err);
       end
    end

    %------------------------------------------------------------------
    % Update function history.
    %------------------------------------------------------------------
    if singleTau || f > sigma^2 / 2 % Don't update if superoptimal.
       lastFv(mod(iter,nPrevVals)+1) = undist(f);
       if fBest > f
          fBest = f;
          % xBest = x;
       end
    end
    
laabel = ['x_',num2str(iter) '.mat'];
save(laabel,'x');

end % while 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Restore best solution (only if solving single problem).
if singleTau && f > fBest
    %% Restoring of xBest is disabled to save memory for large problems
    printf('NOTE: solution not actually optimal, best objective value is %13.7e',sqrt(2*fBest))
    % rNorm = sqrt(2*fBest);
    % printf('\n Restoring best iterate to objective %13.7e\n',rNorm);
    % x = xBest;
    % r = b - Aprod(x,1);
    % g =   - Aprod(r,2);
    % gNorm = options.dual_norm(g,weights);
    % rNorm = norm(r,  2);
end

% Final cleanup before exit.
info.tau         = tau;
info.rNorm       = rNorm;
info.rGap        = rGap;
info.gNorm       = gNorm;
info.rGap        = rGap;
info.stat        = stat;
info.iter        = iter;
info.nProdA      = nProdA;
info.nProdAt     = nProdAt;
info.nNewton     = nNewton;
info.timeProject = timeProject;
info.timeMatProd = timeMatProd;
info.itnLSQR     = itnTotLSQR;
info.options     = options;
info.timeTotal   = toc;

info.xNorm1      = xNorm1(1:iter);
info.rNorm2      = rNorm2(1:iter);
info.lambda      = lambda(1:iter);

% Print final output.
switch (stat)
   case EXIT_OPTIMAL
      printf('\n EXIT -- Optimal solution found\n')
   case EXIT_ITERATIONS
      printf('\n ERROR EXIT -- Too many iterations\n');
   case EXIT_ROOT_FOUND
      printf('\n EXIT -- Found a root\n');
   case {EXIT_BPSOL1_FOUND, EXIT_BPSOL2_FOUND}
      printf('\n EXIT -- Found a BP solution\n');
   case EXIT_LINE_ERROR
      printf('\n ERROR EXIT -- Linesearch error (%i)\n',lnErr);
   case EXIT_SUBOPTIMAL_BP
      printf('\n EXIT -- Found a suboptimal BP solution\n');
   case EXIT_MATVEC_LIMIT
      printf('\n EXIT -- Maximum matrix-vector operations reached\n');
   case EXIT_ACTIVE_SET
      printf('\n EXIT -- Found a possible active set\n');
   case EXIT_AT_PARETO
      printf('\n EXIT -- Reached the pareto curve\n');
   otherwise
      error('Unknown termination condition\n');
end
printf('\n');
printf(' %-20s:  %6i %6s %-20s:  %6.1f\n',...
   'Products with A',nProdA,'','Total time   (secs)',info.timeTotal);
printf(' %-20s:  %6i %6s %-20s:  %6.1f\n',...
   'Products with A''',nProdAt,'','Project time (secs)',timeProject);
printf(' %-20s:  %6i %6s %-20s:  %6.1f\n',...
   'Newton iterations',nNewton,'','Mat-vec time (secs)',timeMatProd);
printf(' %-20s:  %6i %6s %-20s:  %6i\n', ...
   'Line search its',nLineTot,'','Subspace iterations',itnTotLSQR);
printf('\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED FUNCTIONS.  These share some vars with workspace above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function z = Aprod(x,mode)
   if (nProdA + nProdAt >= maxMatvec)
     error('SPGL1:MaximumMatvec','');
   end
     
   tStart = toc;
   if mode == 1
      nProdA = nProdA + 1;
      dispFlag('begin Aprod')
      if   explicit, z = A*x;
      else           z = A(x,1);
      end
      dispFlag('end Aprod')
   elseif mode == 2
      nProdAt = nProdAt + 1;
      dispFlag('begin Atprod')
      if   explicit, z = A'*x;
      else           z = A(x,2);
      end
      dispFlag('end Atprod')
   else
      error('Wrong mode!');
   end
   timeMatProd = timeMatProd + (toc - tStart);
end % function Aprod

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function printf(varargin)
  if logLevel > 0
     fprintf(fid,varargin{:});
  end
end % function printf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = project(x, tau)
   dispFlag('begin Project')
    
   tStart      = toc;
   x = options.project(x,weights,tau); 
   timeProject = timeProject + (toc - tStart);
   
   dispFlag('fin Project')
end % function project

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [iterG,stepG,errG] = ...
    spgLineCurvy(fMax)
% Projected backtracking linesearch.
% On entry,
% g  is the (possibly scaled) steepest descent direction.

EXIT_PBLINE_CONVERGED  = 0;
EXIT_PBLINE_ITERATIONS = 1;
EXIT_PBLINE_NODESCENT  = 2;
gamma  = 1e-4;
stepG   =  1;
scale  =  1;      % Safeguard scaling.  (See below.)
nSafe  =  0;      % No. of safeguarding steps.
iterG   =  0;
debug  =  false;  % Set to true to enable log.
% n      =  length(x);

if debug
   fprintf(' %5s  %13s  %13s  %13s  %8s\n',...
           'LSits','fNew','step','gts','scale');  
end

if isempty(xOld)
    gtxOld = 0;
else
    gtxOld = gStep * dot(g,xOld);
end

while 1
    
    % Calculate scaling of gradient
    g_scale = -stepG*scale*gStep;
    
    % Evaluate trial point
    clear x
    if isempty(xOld)
        x     = project(g_scale .* g, tau);
    else
        x     = project(xOld + g_scale .* g, tau);
    end
    
    % Evaluate function value
    clear r
    r     = b - Aprod(x,1);
    f     = norm(r)^2 / 2;
    gtx   = gStep * dot(g,x);
    gts   = scale * (gtx - gtxOld);
    if gts >= 0 % Should we check real and complex parts individually?
       errG = EXIT_PBLINE_NODESCENT;
       break
    end

    if debug
        fprintf(' LS %2i  %13.7e  %13.7e  %13.6e  %8.1e\n',...
                iterG,f,stepG,undist(gts),scale);
    end
    
    % 03 Aug 07: If gts is complex, then should be looking at -abs(gts).
    if f < fMax - gamma*stepG*abs(gts)  % Sufficient descent condition.
       errG = EXIT_PBLINE_CONVERGED;
       break
    elseif iterG >= lineSrchIt                 % Too many linesearch iterations.
       errG = EXIT_PBLINE_ITERATIONS;
       break
    end
    
    % New linesearch iteration.
    iterG = iterG + 1;
    stepG = stepG / 2;
    
end % while 1
    
end % function spgLineCurvy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dispFlag(flagMsg)
    
    if PRINT_DEBUG_FLAGS
        disp(flagMsg)
        pause(5)
    end
    
end % function dispFlag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of nested functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % function spg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRIVATE FUNCTIONS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nnzX,nnzG,nnzIdx,nnzDiff] = activeVars(x,g,nnzIdx,options)
% Find the current active set.
% nnzX    is the number of nonzero x.
% nnzG    is the number of elements in nnzIdx.
% nnzIdx  is a vector of primal/dual indicators.
% nnzDiff is the no. of elements that changed in the support.
  xTol    = min(.1,10*options.optTol);
  gTol    = min(.1,10*options.optTol);
  gNorm   = options.dual_norm(g,options.weights);
  nnzOld  = nnzIdx;

  % Reduced costs for postive & negative parts of x.
  z1 = gNorm + g;
  z2 = gNorm - g;

  % Primal/dual based indicators.
  xPos    = x >  xTol  &  z1 < gTol; %g < gTol;%
  xNeg    = x < -xTol  &  z2 < gTol; %g > gTol;%
  nnzIdx  = xPos | xNeg;

  % Count is based on simple primal indicator.
  nnzX    = sum(abs(x) >= xTol);
  nnzG    = sum(nnzIdx);
  
  if isempty(nnzOld)
     nnzDiff = inf;
  else
     nnzDiff = sum(nnzIdx ~= nnzOld);
  end
  
end % function spgActiveVars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = LSQRprod(Aprod,nnzIdx,ebar,n,dx,mode)
% Matrix multiplication for subspace minimization.
% Only called by LSQR.
  nbar = length(ebar);
   if mode == 1
      y = zeros(n,1);
      y(nnzIdx) = dx - (1/nbar)*(ebar'*dx)*ebar; % y(nnzIdx) = Z*dx
      z = Aprod(y,1);                            % z = S Z dx
   else
      y = Aprod(dx,2);
      z = y(nnzIdx) - (1/nbar)*(ebar'*y(nnzIdx))*ebar;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fNew,step,rNew,iter,err] = spgLine(f,d,gtd,r,fMax,Aprod,b,feasSrchIt)
% Nonmonotone linesearch.

EXIT_CONVERGED  = 0;
EXIT_LINEITERATIONS = 1;
maxIts = feasSrchIt;
step   = 1;
iter   = 0;
gamma  = 1e-4;
gtd    = -abs(undist(gtd)); % 03 Aug 07: If gtd is complex,
                    % then should be looking at -abs(gtd).
                    
Ad = Aprod(d,1);

while 1

    % Evaluate trial point and function value.
    rNew = r - step*Ad;
    fNew = norm(rNew)^2 / 2;

    % Check exit conditions.
    if fNew < fMax + gamma*step*gtd  % Sufficient descent condition.
       err = EXIT_CONVERGED;
       break
    elseif  iter >= maxIts           % Too many linesearch iterations.
       err = EXIT_LINEITERATIONS;
       break
    end
    
    % New linesearch iteration.
    iter = iter + 1;
    
    % Safeguarded quadratic interpolation.
    if step <= 0.1
       step  = step / 2;
    else
       tmp = (-gtd*step^2) / (2*(fNew-f-step*gtd));
       if tmp < 0.1 || tmp > 0.9*step || isnan(tmp)
          tmp = step / 2;
       end
       step = tmp;
    end
    
end % while 1

end % function spgLine

