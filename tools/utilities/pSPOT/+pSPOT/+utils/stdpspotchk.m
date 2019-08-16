function [opList,m,n,cflag,linear] = stdpspotchk(varargin)
%STDCHK    Performs the standard checking and data extraction routine
%          in pSpot.
%
%   [opList,m,n,cflag,linear] = stdpspotchk(LIST) takes in a cell array of
%   Spot operators or numerical matrices and returns:
%       opList: List of Spot operators in a cell array
%       [m,n] : The array of the sizes of the operators
%       cflag : 1 if complex, 0 if not.
%       linear: 1 if linear, 0 if not.

% Check for empty operators and remove them
ops = ~cellfun(@isempty,varargin);
assert(any(ops),'At least one operator must be specified.');
arrayfun(@(ind) warning('pSpot:NoInput','input "%d" is empty',ind), find(~ops));
opList = varargin(ops);

% Check for pSpot operators
ops = cellfun(@(p) isa(p,'oppSpot'), opList);
assert(~any(ops),' oppSpot operators are not supported');

% Convert all non-Spot arguments to operators
ops = cellfun(@(p) ~isa(p,'opSpot'), opList);
opList(ops) = cellfun(@(p) {opMatrix(p)}, opList(ops));

% Check complexity and setup sizes
[m,n] = cellfun(@size,opList);
real = cellfun(@isreal,opList);
cflag = ~all(real);
linear = cellfun(@(p) logical(p.linear), opList);
linear = all(linear);