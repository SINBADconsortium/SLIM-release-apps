function [c] = lagrange_coe_special(x0);
% function [c] = lagrange_ceo(x0,x);
%
%
%	x0 must between [0,1]. 
%
%	if you want to compute value at off-grid point 5.3 from [4 5 6 7]. just input demical part
%	of 5.3 (which 0.3). this function will give you coefficient [c1 c2 c3 c4] for [4 5 6 7].
%	Then x(5.3) = c1 * x(4) + c2 * x(5) + c3 * x(6) + c4 * x(7);
%	
%
% 
% Author: Xiang Li, Aug 30, 2014



x = [-1 0 1 2];


c = zeros(length(x0),4);


c(:,1) = (x0       ) .* (x0 -  1  ) .* (x0 -  2  ) ./ (-1) ./ (-2) ./ (-3);
c(:,2) = (x0 +  1  ) .* (x0 -  1  ) .* (x0 -  2  ) ./ ( 1) ./ (-1) ./ (-2);
c(:,3) = (x0 +  1  ) .* (x0       ) .* (x0 -  2  ) ./ ( 2) ./ ( 1) ./ (-1);
c(:,4) = (x0 +  1  ) .* (x0       ) .* (x0 -  1  ) ./ ( 3) ./ ( 2) ./ ( 1);


% c = zeros(size(x));
%
% c(1,:) = (x0 - x(2,:)) .* (x0 - x(3,:)) .* (x0 - x(4,:)) ./ (x(1,:) - x(2,:)) ./ (x(1,:) - x(3,:)) ./ (x(1,:) - x(4,:));
% c(2,:) = (x0 - x(1,:)) .* (x0 - x(3,:)) .* (x0 - x(4,:)) ./ (x(2,:) - x(1,:)) ./ (x(2,:) - x(3,:)) ./ (x(2,:) - x(4,:));
% c(3,:) = (x0 - x(1,:)) .* (x0 - x(2,:)) .* (x0 - x(4,:)) ./ (x(3,:) - x(1,:)) ./ (x(3,:) - x(2,:)) ./ (x(3,:) - x(4,:));
% c(4,:) = (x0 - x(1,:)) .* (x0 - x(2,:)) .* (x0 - x(3,:)) ./ (x(4,:) - x(1,:)) ./ (x(4,:) - x(2,:)) ./ (x(4,:) - x(3,:));

%
%
%
%
%
	
	
	
	