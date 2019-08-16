function y = compositeDef(x)
%COMPOSITEDEF Default compositioning of pspot operators
%   
%   y = compositeDef(x) returns a Composite array with the list of ops, x
%   stored in the default distribution configuration given the current
%   parallel poll conditions.
%
%   x must be a 1 by n array or cell

if isa(x,'Composite')
    x = [x{:}];
end

y       = Composite();
glodist = pSPOT.utils.defaultDistribution(length(x));

ind = 0;
for i=1:length(y)
    y{i} = x(ind+1:ind+glodist(i));
    ind  = ind + glodist(i);
end
