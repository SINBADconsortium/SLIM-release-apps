function y = project_tangent(x, dx, dimTree)   
% PROJECT_TANGENT - Project a tangent vector onto the horizontal space along the
% vertical space.
% 
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
% Usage:
%    y = project_tangent(x,dx,dimTree)
%
%  
% Input:
%    x     - current point on the manifold (assumed to be horizontal)
%   dx     - tangent vector
%
% Output:
%    y     - horizontal vector
    y = dx - project_vertical(x,dx,dimTree);
end
