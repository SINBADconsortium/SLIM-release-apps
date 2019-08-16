function y = global_sum(varargin)
%GLOBAL SUM    (SPMD only) Global Summation
%   y = global_sum(x) returns the addition of the x's from each lab.
%   The result is replicated on all labs.
%
%   y = global_sum(x, labtarget) places all of the results on lab labtarget. Y
%   will be equal to [] on all other labs.
%
%   This is actually a wrapper to the PCT funcion gplus which is more
%   performance-y compared to our old implementation of the concept.
%
%   type "help gplus" for more information and examples

y = gplus(varargin{:});