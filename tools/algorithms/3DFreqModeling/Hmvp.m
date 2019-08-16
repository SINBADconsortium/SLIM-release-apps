%  Matrix vector multiplication with the matrix coming from 
%  helmholtz_3d_operto discretization. It calls a mex function for 
%  performing the matrix-vector product itself. Before using this function, 
%  make sure to run.
%  
%  $ mex Hmvp_mex.c
%  
%  in order to obtain the executable file Hmvp_mex.mexa64.
%  
%  USE:
%    b = Hmvp(H,idx,x)
%  
%  INPUT:
%    H   - coming from helmholtz_3d_operto discretization.
%    idx - offsets of diagonals stored in H
%    x   - vector to multiply
%
% OUTPUT:
%    b   - result of the product H*x.
%
% AUTHOR: Art Petrenko (edited by Rafael Lago, Curt Da Silva)
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
%-----------------------------------------------------------------------------
function output = Hmvp(H,idx,x,nthreads)
   assert(~isdistributed(x) && ~iscodistributed(x));
   if exist('nthreads','var')==0 || isempty(nthreads) || nthreads==0, nthreads = size(x,2); end
   output = Hmvp_MT_mex(H,idx,x,nthreads);
         
end
