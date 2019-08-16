% Compute the norm of a distributed vector
%
% USE:
%   out = pnorm(u,N)
%
%   Computes out = ||u||, where u is a distributed vectors with local
%   dimension N.
%
% AUTHOR: Rafael Lago
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2014
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
% 
%-----------------------------------------------------------------------------
function out = pnorm(u,N) 
   
   uloc = getLocalPart(u);
   out = norm(uloc((labindex~=1)*2*N+1 : end));
   out = sqrt(gplus(out.^2));
   
end