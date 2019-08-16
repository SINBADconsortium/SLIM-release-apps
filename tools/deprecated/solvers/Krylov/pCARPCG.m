function [x,res,niter] = pCARPCG(R,idx,q,x,n,options)
% parallel CARPCG algorithm for band-storage matrices
%
% use:
%   x = pCARPCG(R,idx,q,x0,options)
%
% input:
%   {R,idx} - normalized, distributed, matrix in band storage format, see mat2R
%   q       - distributed right hand side
%   x0      - distributed initial guess
%   n       - size of domain
%   options.w     - relaxation parameter, default = 1.5
%   options.tol   - CG tolerance, default = 1e-6
%   options.maxit - max. iterations, default = length(q)
%   options.n_threads - number of CARP blocks to use for an additional
%   parallelization over threads. Each MATLAB worker will split up its block
%   into options.n_threads sub-blocks
%
% output:
%   x   - solution  
%   res - residual vector    
%
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


% get number of workers
spmd
	np = numlabs;
end
np = np{1};
if np==1
    error('Only one worker available, use serial version');
end

% parse options
w     = 1.5;
tol   = 1e-6;
maxit = length(q);

if isfield(options,'w')
    w = options.w;
end
if isfield(options,'tol')
    tol = options.tol;
end
if isfield(options,'maxit')
    maxit = options.maxit;
end

%%% Note: the following section also appears in pCARPCG.m in this directory.

% If supplied in the call to the function, that value of n_threads takes
% precedence.
if isfield(options, 'n_threads')
	n_threads = options.n_threads;
% Otherwise use the environment variable.
else
	n_threads = getenv('OMP_NUM_THREADS');
end
% If after all that n_threads is still not a valid value, fall back to the
% default.
if not(n_threads >= 1)
	n_threads = 1;
end
%%% End of copy-and-pasted section.

N = prod(n(1:end-1));

spmd
	% Make sure diagonals of matrix are stored as rows
	if size(R,1)==length(idx)
		% do nothing
	elseif size(R,2)==length(idx)
		R = R.';
	else
		error('Dimensions of R and idx do not match.');
	end

	x = pSPOT.pWindow.funWindowLast1HaloMakeCodist(x,n(end),np,1);
	q = pSPOT.pWindow.funWindowLast1HaloMakeCodist(q,n(end),np,1);

	b = dsweep(R,idx,0*x,q,w,n,n_threads);  
	assert(not(any(isnan(getLocalPart(b)))), 'pCARPCG: an element of b, the right hand side, is NaN after the initial sweep.');

	% start CG 
	bloc = getLocalPart(b);
	normb = norm(bloc((labindex~=1)*2*N+1 : end));
	normb = sqrt(gplus(normb.^2));
	tol  = (tol*normb)^2;
	res  = [];
	k    = 0;

	r = b - x + dsweep(R,idx,x,0*b,w,n,n_threads);
	p = r;
	rloc = getLocalPart(r);
	normr = norm(rloc((labindex~=1)*2*N+1 : end))^2;
	normr = gplus(normr);
	res   = [res normr];

	while (k<=maxit)&&(normr>tol)    

		q = p - dsweep(R,idx,p,0*b,w,n,n_threads);
		
		ploc = getLocalPart(p);
		qloc = getLocalPart(q);
		pdotq = ploc((labindex~=1)*2*N+1 : end)'*qloc((labindex~=1)*2*N+1 : end);
		pdotq = gplus(pdotq);
		% alpha needs to be codistributed to make x codistributed also (as
		% opposed to a composite)
		alpha = codistributed(normr/pdotq);
		
		x = x + alpha*p;
		r = r - alpha*q;
		
		rloc = getLocalPart(r);
		normr = norm(rloc((labindex~=1)*2*N+1 : end))^2;
		normr = gplus(normr);
		beta  = codistributed(normr/res(end));
		
		p = r + beta*p;
		
		k = k + 1;
		res = [res normr];
	end
	x = pSPOT.pWindow.funWindowLast1HaloDropCodist(x,n(end),np,1);
end % spmd

res   = sqrt(res{1})/normb{1};
niter = k{1};

end % function pCARPCG


function x = dsweep(R,idx,x,b,w,n,n_threads)

    N  = prod(n(1:end-1));
    
    % setup
    dist = getCodistributor(x);
    Rloc = getLocalPart(R);
    xloc = getLocalPart(x);
    bloc = getLocalPart(b);
    idxloc = idx + (labindex~=1)*N;

    % forward sweep
    xloc = sweepR_mex(Rloc,idxloc,xloc,bloc,w,1,n_threads);  
	assert(not(any(isnan(xloc))),'dsweep: an element of xloc, the local part of x, is NaN after the forward sweep.');
    
    % get halos
    if labindex<numlabs; labTo = labindex + 1; else labTo = []; end;
    if labindex>1; labFrom = labindex - 1; else labFrom = []; end;
    halol = labSendReceive(labTo, labFrom, xloc(end-2*N+1:end));
    
    if labindex>1; labTo = labindex - 1; else labTo = []; end;
    if labindex<numlabs; labFrom = labindex + 1; else labFrom = []; end;
    halor = labSendReceive(labTo, labFrom, xloc(1:2*N));
    
    % update halos
    if labindex > 1; xloc(1:2*N)=(xloc(1:2*N)+halol)/2; end;
    if labindex < numlabs; xloc(end-2*N+1:end)=(xloc(end-2*N+1:end)+halor)/2; end;
    
    % backward sweep
    xloc = sweepR_mex(Rloc,idxloc,xloc,bloc,w,-1,n_threads);
	assert(not(any(isnan(xloc))),'dsweep: an element of xloc, the local part of x, is NaN after the backward sweep.');
    
    % get halos
    if labindex<numlabs; labTo = labindex + 1; else labTo = []; end;
    if labindex>1; labFrom = labindex - 1; else labFrom = []; end;
    halol = labSendReceive(labTo, labFrom, xloc(end-2*N+1:end));
    
    if labindex>1; labTo = labindex - 1; else labTo = []; end;
    if labindex<numlabs; labFrom = labindex + 1; else labFrom = []; end;
    halor = labSendReceive(labTo, labFrom, xloc(1:2*N));
    
    % update halos
    if labindex > 1; xloc(1:2*N)=(xloc(1:2*N)+halol)/2; end;
    if labindex < numlabs; xloc(end-2*N+1:end)=(xloc(end-2*N+1:end)+halor)/2; end;
    
    % wrap up
    x = codistributed.build(xloc,dist,'noCommunication');

end
