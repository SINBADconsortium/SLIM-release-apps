function test_suite = test_ToolsSolversQuasiNewton
% test for lbfgs.m
%

% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

%
initTestSuite;

function setup


function teardown
delete('iter.log');

function testLBFGS
% test lbfgs on simple quadratic function
% added June 15, 2011

N  = 100;

A  = randn(N);
x0 = randn(N,1);
fh = @(x)qm(x,A,1:N);

M = 5000;

options.itermax = 5000;
options.M=M;
options.tol = 1e-15;
options.fid = fopen('iter.log','w');
xn = lbfgs_wWolfe(fh,x0,options);

assertVectorsAlmostEqual(xn,zeros(N,1),'absolute',1e-14);

function [f,g,H] = qm(x,A,I)
% Quadratic model x'*A'*A*x for testing of optimization codes
%
% Tristan van Leeuwen, 2011
% tleeuwen@eos.ubc.ca
%
% use:
%   [f,g,H] = qm(x,A,I);
%
% input:
%   x - vector 
%   A - matrix
%   I - index set determining which rows of A to use.
%
% output:
%   f - value
%   g - gradient
%   H - Hessian

H = A(I,:)'*A(I,:);

f = .5*x'*H*x;
g = H*x;
