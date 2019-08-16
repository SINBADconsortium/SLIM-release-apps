% Compute the dot product between two vectors in parallel
%
% USE:
%   out = pdot(u,v,N)
%
%   Computes out = u'*v, where u and v are distributed vectors with local
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
function out = pdot(u,v,N) 
   
   uloc = getLocalPart(u);
   vloc = getLocalPart(v);
   out = uloc((labindex~=1)*2*N+1 : end)'*vloc((labindex~=1)*2*N+1 : end);
   out = gplus(out);
   
end