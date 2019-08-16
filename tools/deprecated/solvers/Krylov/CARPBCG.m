function [X,res,niter] = CARPBCG(R,idx,Q,X,n,options)
% serial CARP- BlockCG algorithm for band-storage matrices
%
% use:
%   [x,res] = CARPBCG(R,idx,Q,x0,n,options)
%
% input:
%   {R,idx} - normalized matrix in band storage format, see mat2R
%   Q       - right hand side matrix
%   X0      - initial guess matrix
%   n       - size of domain
%   options.w     - relaxation parameter, default = 1.5
%   options.tol   - CG tolerance, default = 1e-6
%   options.maxit - max. iterations, default = length(q)
%   options.nb    - number of blocks, default = 1.
%
% output:
%   X   - solution matrix
%   res - residual vector    
% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2013
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

% parse options
w     = 1.5;
tol   = 1e-6;
maxit = size(Q,1);

if isfield(options,'w')
    w = options.w;
end
if isfield(options,'tol')
    tol = options.tol;
end
if isfield(options,'maxit')
    maxit = options.maxit;
end
if isfield(options,'nb')
    nb = options.nb;
end

res = zeros(maxit,nb);

% Make sure diagonals of matrix are stored as rows
if size(R,1)==length(idx)
	% do nothing
elseif size(R,2)==length(idx)
	R = R.';
else
	error('Dimensions of R and idx do not match.');
end

% construct r.h.s. and operator
B    = dsweep(R,idx,w,X,Q);
Afun = @(X)(X - dsweep(R,idx,w,X,zeros(size(B))));

I = [1:floor(size(B,2)/nb):size(B,2) size(B,2)+1];
% start CG 
for i = 1:nb
    Ii    = I(i):I(i+1)-1;
    normB = norm(B(:,Ii),'fro');
    
    k = 1;

    R = B(:,Ii) - Afun(X(:,Ii));
    P = R;
    RR = R'*R;
    res(k,i) = trace(RR);

    while (k<maxit)&&(res(k,i)>(tol*normB)^2)    

        Q  = Afun(P);

        alpha = (P'*Q)\(RR);

        X(:,Ii) = X(:,Ii) + P*alpha;
        R = R - Q*alpha;

        RRt  = R'*R;
        beta = RR\RRt;
        RR   = RRt;

        P = R + P*beta;

        k = k + 1;
        res(k,i) = trace(RR);
    end

    res(:,i) = sqrt(res(:,i))/normB;
    niter(i) = k;
end
niter = mean(niter);

end


function X = dsweep(R,idx,w,X,B)
    for k=1:size(X,2);
        % forward sweep
        X(:,k) = sweepR_mex(R,idx,X(:,k),B(:,k),w,1);

        % backward sweep
        X(:,k) = sweepR_mex(R,idx,X(:,k),B(:,k),w,-1);
    end
end