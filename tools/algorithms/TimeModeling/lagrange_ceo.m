function [c] = lagrange_ceo(x0,x);
% function [c] = lagrange_ceo(x0,x);
% this function computes lagrange coefficients of a single point (cubic interpatation)
% if you want to comput x0 with x1-x4 as follows
%                    x0
%                    |
% x1 --------- x2 -------- x3 ------- x4
% 
% x0 = c1 * x1 + c2 * x2 + c3 * x3 + c4 * x4
% then output of this function with give you
% c = [c1 c2 c3 c4]
% 
% Author: Xiang Li, OCT 30, 2013

if x0  <  x(2) | x0 > x(3)
	disp(['x0 need in between x(2:3)'])
end

if isvector(x)
	x =x(:);
end

c = zeros(size(x));

c(1,:) = (x0 - x(2,:)) .* (x0 - x(3,:)) .* (x0 - x(4,:)) ./ (x(1,:) - x(2,:)) ./ (x(1,:) - x(3,:)) ./ (x(1,:) - x(4,:));
c(2,:) = (x0 - x(1,:)) .* (x0 - x(3,:)) .* (x0 - x(4,:)) ./ (x(2,:) - x(1,:)) ./ (x(2,:) - x(3,:)) ./ (x(2,:) - x(4,:));
c(3,:) = (x0 - x(1,:)) .* (x0 - x(2,:)) .* (x0 - x(4,:)) ./ (x(3,:) - x(1,:)) ./ (x(3,:) - x(2,:)) ./ (x(3,:) - x(4,:));
c(4,:) = (x0 - x(1,:)) .* (x0 - x(2,:)) .* (x0 - x(3,:)) ./ (x(4,:) - x(1,:)) ./ (x(4,:) - x(2,:)) ./ (x(4,:) - x(3,:));

	
	
	
	
	
	
	
	
	