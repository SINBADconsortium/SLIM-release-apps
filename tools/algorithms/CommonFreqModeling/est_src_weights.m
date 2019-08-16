function w = est_src_weights(w,Dpred,Dobs,src_freq_idx,misfit_func)
% EST_SRC_WEIGHTS - Estimate per-frequency source weights
%
% Curt Da Silva, 2016
%
% Usage:
%   w = est_src_weights(w,Dpred,Dobs,src_freq_idx,misfit_func);
%
% Input:
%   w            - initial estimate of weights, nfreq x 1 vector
%   Dpred        - predicted data (serial or distributed matrix, nrec x npdes)
%   Dobs         - observed data
%   src_freq_idx - N x 2 matrix of [source_idx, freq_idx] corresponding to 
%                  each column of Dpred/Dobs
%   misfit_func  - pointwise misfit function handle, [f,g,h] = misfit(x,y);
%  
% Output:
%   w           - estimated weights, that minimize, per-frequency, misfit(w*Dpred,Dobs)
%
    f = src_freq_idx(:,2);
    nrec = size(Dpred,1);

    if isdistributed(Dpred)
        W = @(t) oppKron2Lo(opDiag_swp(t(f)),opDirac(nrec));
        A = oppKron2Lo(opPartialSum(src_freq_idx,2),opOnes(1,nrec));
    else
        W = @(t) opKron(opDiag_swp(t(f)),opDirac(nrec));
        A = opKron(opPartialSum(src_freq_idx,2),opOnes(1,nrec));
    end
    Dpred = vec(Dpred);
    Dobs = vec(Dobs);

    function [f,g,h] = srcw(t)                
        [f,g,h] = misfit_func(W(t)*Dpred,Dobs);
        g = A*(g.*conj(Dpred));
        h = A*((h*conj(Dpred)).*Dpred);
        if isdistributed(g)
            g = gather(g); h = gather(h); 
        end
    end
    
    w = newton(@srcw,w,1e-6);    
end

function [x] = newton(fh, x0, tol)
% Newton's method for problems with a strictly diagonal Hessian.
% 
% function [x] = newton(func, x0, tol)
%
% Input:  func  - function handle returning [value, gradient, hessian]
%                 where hessian is a vector containing the diagonal of the
%                 Hessian.
%         x0    - initial value of parameter
%         tol   - tolerance

% initialization
maxIter = 10;
mult    = 0.8;
c       = 0.01;
x       = x0;
iter    = 0;

[f, g, h] = fh(x);

while(norm(g) > tol)&&(iter < maxIter)
    
    % search direction
    d = -g./h;
    
    % Armijo linesearch
    step   = 1;
    lsiter = 0;
    while( fh(x + step*d) > f - c*step*abs(g'*d) )
        step = step*mult;
        lsiter = lsiter + 1;
        if lsiter > 20
            return;
        end
    end
    
    % update
    x = x + step*d;
    iter = iter + 1;
    [f, g, h] = fh(x);
 
end  % while
end