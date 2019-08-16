function [y,res] = linearsolve(op,b,x0,mode,lsopts)    
% LINEARSOLVE - Iterative solver setup + solve for a spot operator
%  
%  Curt Da Silva, 2015
%
%  Usage:
%     y = linearsolve(op,b,x0,mode,lsopts);
%   
%  Input:
%     op     - SPOT operator (forward operator, not op')
%     b      - right hand side
%     mode   - if mode==1, forward solve, otherwise adjoint solve
%     x0     - initial guess (default: 0)
%     lsopts - LinSolveOpts object specifying solver parameters
%              (default: op.solve_opts)
%
%  Output:
%     y      - solution to op*y = b if mode==1 or op'*y = b if mode ~=1
%
    isprop = @(obj,x) ~isempty(find(ismember(properties(obj),x),1));
    if exist('x0','var')==0 || isempty(x0)
        x0 = zeros(size(b));
    end
    
    assert(exist('lsopts','var') && isa(lsopts,'LinSolveOpts'),'lsopts must be a LinSolveOpts object');
    assert(size(op,1)==size(op,2),'op must be square');
    assert(size(b,1)==size(op,1),'Rhs dimension mismatch');
    assert(size(x0,1)==size(op,1),'x0 dimension mismatch');
    
    % Create a copy of lsopts, 
    if isprop(lsopts,'output_freq')
        output_freq = lsopts.output_freq;
    else output_freq = 0; 
    end
    lsopts = copy(lsopts);
    if output_freq > 0, 
        lsopts.addprop('output_freq'); 
        lsopts.output_freq = output_freq; 
    end 
    
    is_solver = @(x) strcmp(lsopts.solver,x);
    is_prec = @(x) strcmp(lsopts.precond,x);
            
    if mode == 1
        A = op;
    else
        A = op';
    end
        
    % Construct preconditioner
    % Construct CRMN/CGMN solver
    if is_solver(LinSolveOpts.SOLVE_CRMN) || is_solver(LinSolveOpts.SOLVE_CGMN)
        if isa(lsopts.precond,'char') && is_prec(LinSolveOpts.PREC_KACZSWP)
            w = 1.5; nsweeps = 1; zb = 0*b;
            b = kaczswp(op,w,zb,b,nsweeps,mode);
            A = opDirac(op.n) - opFunction_swp(op.n,op.n,@(y,m) kaczswp(op,w,y,zb,nsweeps,mode),1);
            if is_solver(LinSolveOpts.SOLVE_CRMN)
                solve = @(A,b) CRMN(A,b,x0,lsopts);
            else
                solve = @(A,b) CGMN(A,b,x0,lsopts);
            end
        else
            if is_solver(LinSolveOpts.SOLVE_CRMN)
                solve = @(A,b) CRMN(A'*A,A'*b,x0,lsopts);
            else
                solve = @(A,b) CGMN(A'*A,A'*b,x0,lsopts);
            end                
        end
    elseif isa(lsopts.precond,'LinSolveOpts')
        par = struct;
        par.tol = lsopts.precond.tol;
        par.maxit = lsopts.precond.maxit;
        par.maxinnerit = lsopts.precond.maxinnerit;
        par.precond = @(x) x;
        switch lsopts.precond.solver
          case LinSolveOpts.SOLVE_BICGSTAB
            precond = @(x) BICGSTAB(A,x,0*x,par);   
          case LinSolveOpts.SOLVE_JACOBI
            alpha = lsopts.precond.params.alpha;
            precond = @(x) jacobi(op,alpha,x,0*x,mode,par.maxit);
          case LinSolveOpts.SOLVE_FGMRES
            precond = @(x) FGMRES(A,x,0*x,par);
          case LinSolveOpts.SOLVE_CRMN
            if isa(lsopts.precond.precond,'char') && strcmp(lsopts.precond.precond,LinSolveOpts.PREC_KACZSWP)
                w = 1.5; nsweeps = 1; zb = 0*b;
                A1 = opDirac(op.n) - opFunction_swp(op.n,op.n,@(y,m) kaczswp(op,w,y,zb,nsweeps,mode),1);                
                precond = @(x) CRMN(A1,kaczswp(op,w,zb,x,nsweeps,mode),zb,lsopts);
            else
                precond = @(x) CRMN(A'*A,A'*x,0*x,lsopts);
            end
          otherwise 
            error('Unrecognized preconditioner');
        end        
    elseif isa(lsopts.precond,'char')
        switch lsopts.precond
          case LinSolveOpts.PREC_IDENTITY
            precond = opDirac(size(A,1));          
          otherwise
            error('Unrecognized preconditioner');
        end
    elseif isa(lsopts.precond,'opSpot')
        if mode==1
            precond = lsopts.precond;
        else
            precond = lsopts.precond';
        end
        assert(size(precond,1)==size(op,2),'Preconditioner size mismatch');
    elseif isa(lsopts.precond,'function_handle')
        T = lsopts.precond;
        if mode==1            
            precond = @(x) T(x,1);
        else
            precond = @(x) T(x,0);
        end
        
    elseif isa(lsopts.precond,'FuncObj')
        precond = lsopts.precond;
        if precond.hasarg('mode')
            precond.setarg('mode',mode);
            precond.setarg('x',x0);
        end
    end
    
    % Actual solver setup
    switch lsopts.solver
      case LinSolveOpts.SOLVE_JACOBI
        if isfield(lsopts.params,'alpha')
            y = jacobi(op,lsopts.params.alpha,b,x0,mode,lsopts.maxit);
        else
            y = jacobi(op,[],b,x0,mode,lsopts.maxit);
        end
        res = 0;
      case LinSolveOpts.SOLVE_BICGSTAB        
        lsopts.precond = precond;        
        y = BICGSTAB(A,b,x0,lsopts);
        res = 0;
      case LinSolveOpts.SOLVE_GMRES
        [y,~,~,~,res] = gmres(A,b,lsopts.maxinnerit,lsopts.tol,lsopts.maxit,opDirac(length(b)),precond);
      case LinSolveOpts.SOLVE_FGMRES
        lsopts.precond = precond;        
        [y,res] = FGMRES(A,b,x0,lsopts);      
      case LinSolveOpts.SOLVE_LU
        res = 0;
        if isfield(op.params,'is_full') && op.params.is_full
            L = op.params.L; U = op.params.U; P = op.params.P;
            if mode==1
                y = U\(L\(P*b));
            else
                y = P'*(L'\(U'\b));
            end
        else
            L = op.params.L; U = op.params.U; P = op.params.P; Q = op.params.Q; R = op.params.R;
            if mode == 1
                y = Q*(U\(L\(P*(R\(b))))); 
            else
                y = R'\(P'*(L'\(U'\(Q'*b))));
            end
        end
      case LinSolveOpts.SOLVE_BACKSLASH
        coef = op.params.coef;
        res = 0;
        if mode == 1
            y = coef\b; 
        else
            y = coef'\b;
        end
      case {LinSolveOpts.SOLVE_CRMN,LinSolveOpts.SOLVE_CGMN}
        [y,res] = solve(A,b);
        
      otherwise
        error(['Unrecognized solver ' lsopts.solver]);
    end
  
    
   
end