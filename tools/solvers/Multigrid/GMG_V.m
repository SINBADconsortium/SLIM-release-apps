%------------------------------------------------------------------------------
% Generic geometric multigrid routine; it merely uses function handles to do
% everything. This is meant to be used for band-storage matrices as well.
%
% USE:
%   x = GMG_V(A,b,x,par,mode)
%
% INPUT:
%   A         - explicit matrix or spot operator
%   b         - right hand side
%   x         - initial guess (default: [], zero vector)
%   par.maxit - Maximum number of iterations to be performed. Usually this is
%               set to one or two
%
% OUTPUT:
%   x   - approximate solution 
% 
% AUTHOR: Rafael Lago
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: August, 2014
%
% Updated by
%   Curt Da Silva, 2015
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
% 
%------------------------------------------------------------------------------
function x = GMG_V(A,b,x,par,mode)

    if mode==0
        A = A';
    end
    
    % Get initial residual (if needed)
    if isempty(x) || nnz(x)==0
        x = zeros(size(b)); 
    end
    
    if ~isfield(par,'pst_smoother')
        par.pst_smoother = par.pre_smoother;
    end
    
    R = par.restriction;
    coarse = par.coarse;
    P = par.prolongation;
    S_pre = par.pre_smoother;
    S_pst = par.pst_smoother;
    l=1;
    while (l<=par.maxit)
        x = S_pre(b,x,mode);
        r = b - A*x;
        xc = R*r;   
        xc = coarse(xc,0*xc,mode);
        x = x + P*xc;  
        x = S_pst(b,x,mode);      
        l = l + 1;
    end
    
end
