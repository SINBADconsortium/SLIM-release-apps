function y = defGlobInd(x)
%DEFGLOBIND  Default Global Indices corresponding to codistributor
%
%   y = defGlobInd(x) gives you the default global indices in a cell array
%   that "tells you the relationship between indices on a local
%   part and the corresponding indices on the distributed array.  The
%   global indices method on a codistributor object allows you to get this
%   relationship without actually creating the array."

defdist   = pSPOT.utils.defaultDistribution(x);
total     = 0;
y         = cell(1,length(defdist));
for i=1:length(defdist)
    y{i}  = total+1:total+defdist(i);
    total = total + defdist(i);
end