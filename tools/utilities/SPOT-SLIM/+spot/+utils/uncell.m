function y = uncell(x)
%UNCELL  Uncellification utility
%
%   y = uncell(x) uncells cell array x to its most un-nested state and
%   return it as a vector.

y = [];

for i = 1:length(x)
    if iscell(x(i))
        y = [y spot.utils.uncell(x{i})];
    else
        y = [y x(i)];
    end
end