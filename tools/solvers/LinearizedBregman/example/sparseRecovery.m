%% Sparse recovery example for the linearized Bregman method
%
% Author: Philipp Witte
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%         
% Date: May, 2015
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.



clear; clc; close all;

% dimensions of underdetermined system
m = 2000;
n = 4000;
A = randn(m,n);

% sparse solution solution with 1% non-zero entries
x = diag(rand(n,1)<=0.01)*randn(n,1); 

% measurements
bTrue = A*x;
% add noise
b = bTrue + 0.5*rand(size(bTrue));
sigma = norm(bTrue-b,2);


%% Set-Up Optimization
options.verbosity = 1;
options.iterations = 100;
options.xtrue = x;
x0 = [];
tau = 0;
lambda = 1;

% spgl1 solution
tic
[x1, res1, g1, info1] = spgl1(A,b,tau,sigma,x0,options);
t1=toc;
 
% linearized Bregman with dynamic step size (default)
tic
[x2, res2, info2] = lbm(A,b,lambda,sigma,x0,options);
t2=toc;

% linearized Bregman with exact step size
options.step = 'e';
tic
[x3, res3, info3] = lbm(A,b,lambda,sigma,x0,options);
t3=toc;


%% Plots

ymin = 1.3*min(x);
ymax = 1.3*max(x);
xax = 1:length(x);
FigHandle1 = figure('Position', [100, 100, 1049, 695]);

% SPGL1 solution
subplot(3,1,1);
plot(x); title(['SPGl1 (SNR = ',num2str(SNR(x,x1),'%4.2f'),', t = ',...
    num2str(t1,'%4.2f'),' s)']);
axis([0 m ymin ymax])
hold on
for i=1:length(x)
if (x1(i) ~= 0)
plot(xax(i),x1(i),'ro')
end
end
legend('True x', 'SPGl1','location','SouthWest')

% Linearized Bregman dynamic step solution
subplot(3,1,2);
plot(x); title(['LB (SNR = ',num2str(SNR(x,x2),...
    '%4.2f\n'),', t = ',num2str(t2,'%4.2f'),' s)']);
axis([0 m ymin ymax])
hold on
for i=1:length(x)
if (x2(i) ~= 0)
plot(xax(i),x2(i),'ro')
end
end
legend('True x', 'LB','location','SouthWest')

% Linearized Bregman exact step solution
subplot(3,1,3);
plot(x); title(['LB exact (SNR = ',num2str(SNR(x,x3),...
    '%4.2f\n'),', t = ',num2str(t3,'%4.2f'),' s)']);
axis([0 m ymin ymax])
hold on
for i=1:length(x)
if (x3(i) ~= 0)
plot(xax(i),x3(i),'ro')
end
end
legend('True x', 'LB exact','location','SouthWest')
