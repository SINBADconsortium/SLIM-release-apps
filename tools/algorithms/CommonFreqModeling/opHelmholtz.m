classdef opHelmholtz < opSpot & handle
% opHelmholtz - SPOT wrapper for the Helmholtz operator
%
%  Curt Da Silva, 2015
%
%
    properties( SetAccess = protected )
       mode,params,deriv_mode;
    end
    properties
        solve_opts,disp_output;
    end
    methods
        function op = opHelmholtz(mode,params,solve_opts,deriv_mode,disp_output)
            if exist('deriv_mode','var')==0
                deriv_mode = false;
            end
            if strcmp(mode,'explicit')
                Nt = size(params.coef,2);
                if strcmp(solve_opts.solver,LinSolveOpts.SOLVE_LU)
                    [L,U,P,Q,R] = lu(params.coef);
                    params.L = L; 
                    params.U = U;
                    params.P = P;
                    params.Q = Q;
                    params.R = R;
                end
            elseif strcmp(mode,'operto27');
                Nt = prod(params.nt);
                if deriv_mode==0
                    wn = param2wavenum(params.v,params.freq,params.unit);
                elseif deriv_mode == 1
                    [~,wn] = param2wavenum(params.v,params.freq,params.unit);
                elseif deriv_mode == 2
                    [~,~,wn] = param2wavenum(params.v,params.freq,params.unit);
                end                
                params.wn = wn;
                params = rmfield(params,'v');
            else
                error('Unrecognized mode');
            end
                        
            op = op@opSpot('Helmholtz operator',Nt,Nt);
            op.mode = mode;
            op.params = params;            
            op.deriv_mode = deriv_mode;
            op.solve_opts = solve_opts;
            op.cflag = true;
            op.sweepflag = true;   
            if exist('disp_output','var')
                op.disp_output = disp_output;
            else
                op.disp_output = false;
            end
        end
    end
    
    methods        
        function x = kacz_sweep(op,w,x,b,nsweeps,mode)                        
            p = op.params;
            for k=1:nsweeps            
                if mode==1
                    if isreal(p.wn)
                        x = Helm3d_kacz_sweep_mex(p.wn,p.h,p.nt,p.npml,w,x,b);
                    else
                        x = Helm3d_kacz_sweep_wnc_mex(p.wn,p.h,p.nt,p.npml,w,x,b);
                    end
                else
                    x = Helm3d_kacz_sweep_adj_mex(p.wn,p.h,p.nt,p.npml,w,x,b);
                end           
            end            
        end        
   end        
    
    methods ( Access = protected )
        function y = multiply(op,x,mode)
            switch op.mode
              case 'explicit'
                if mode == 1
                    y = op.params.coef*x;
                else
                    y = op.params.coef'*x;
                end
              case 'operto27'     
                p = op.params;                
                if mode==1
                    if op.deriv_mode > 0
                        y = Helm3dmvp_forw_deriv_mex(p.wn,p.h,p.nt,p.npml,x,p.n_threads);
                    else
                        if isreal(p.wn)                                                
                            y = Helm3dmvp_forw_mex(p.wn,p.h,p.nt,p.npml,x,p.n_threads);
                        else
                            y = Helm3dmvp_forw_wnc_mex(p.wn,p.h,p.nt,p.npml,x,p.n_threads);
                        end
                    end
                else
                    if op.deriv_mode > 0
                        y = Helm3dmvp_adj_deriv_mex(p.wn,p.h,p.nt,p.npml,x,p.n_threads);
                    else                    
                        if isreal(p.wn)
                            y = Helm3dmvp_adj_mex(p.wn,p.h,p.nt,p.npml,x,p.n_threads);
                        else
                            y = Helm3dmvp_adj_wnc_mex(p.wn,p.h,p.nt,p.npml,x,p.n_threads);
                        end
                    end
                end
            end
        end
        
        function y = divide(op,x,mode)
            if norm(x) < 1e-12, y = zeros(size(x));
            else
                y = linearsolve(op,x,[],mode,op.solve_opts);
            end
        end
        
    end
end