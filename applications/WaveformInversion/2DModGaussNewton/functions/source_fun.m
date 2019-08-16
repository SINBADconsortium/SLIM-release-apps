function  [b,simp,G] = source_fun(gl,nx,nz,nx_border,nz_border,ns,sp,sdep,stype,p) 
% This function define b and G for helmholtz equation.
% Frequency domain wavefield simulation needs a linear equation to be solved.
% for example
% 	u = A\b
% 	where u is wavefield, A is helmholtz matrix and b is source signature.
% 	b is usually a diagnal matrix, each columne of it represents one sequential
% 	shot. We also can do simultaneous shots by multiply b from the right hand 
% 	side with a matrix G.
% 	
% Inputs
% 	gl: grid length
% 	nx: number of grids in x dimension
% 	nz: number of grids in z dimension
% 	nx_border: number of grids in the boundary in x dimension 
% 	nz_border: number of grids in the boundary in z dimension
% 	ns: number of shots
% 	sp: source position 
% 
% Outputs	
% 	b: source
% 	simp: source position on the grid
% 
% Author: Xiang Li
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: 02, 2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

Q = (speye(ns)./(gl^2));

if stype == 1
	G   = opGaussian(p,ns);
	b   = Q * G';
	simp  = [];
elseif stype == 2
	simp= randperm(ns);
	simp= sort(simp(1:p));
	idx = sub2ind([p,ns],1:p,simp);
	G   = spalloc(p,ns,p);
	G(idx) = 1;
	b   = Q * G';
end


return

% 
% % if ~exist('israndn'),israndn = 0;end
% if ~exist('p'),p = ns;end
% % hx = xd/(nx-1);
% % hz = zd/(nz-1);
% hx   = gl;
% hz   = gl;
% nx_tot = nx + 2 * nx_border;
% nz_tot = nz + 2 * nz_border;
% % Q  = zeros(nz,nx,ns);
% % Q(1,1:nx,1:ns) = eye(ns)/(hx*hz);
% Q  = spalloc(nz_tot,nx_tot*ns,ns);
% 
% idx = sp + nx_border;
% idx = sub2ind([nx_tot,ns],idx,1:ns);
% 
% Q(sdep+nz_border,idx) = 1/(hx*hz);
% Q  = Q(:);
% if stype == 1
% 	G   = opGaussian(p,ns);
% 	RM  = opKron(G,speye(nx_tot),speye(nz_tot));
% 	b   = RM * Q;
% 	sname = ['Sim'];
% 	simp  = [];
% elseif stype == 2
% 	simp= randperm(ns);
% 	simp= sort(simp(1:p));
% 	idx = sub2ind([p,ns],1:p,simp);
% 	G   = spalloc(p,ns,p);
% 	G(idx) = 1;
% 	RM  = opKron(G,speye(nx_tot),speye(nz_tot));
% 	b   = RM * Q;
% 	b   = sparse(b);
% 	sname = ['Seq'];
% end
% b = reshape(b,nx_tot*nz_tot,p);
% 
% (speye(model.nshot)./(model.gl^2)) 


