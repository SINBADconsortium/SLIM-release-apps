% Seismic trace interpolation in common receiver records with curvelets

% Author: Philipp Witte
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%         
% Date: May, 2015

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.


clear; 
clc;

%% Load seismic data and geometry
load ../data/receiver_gather.mat
load ../data/SuTraceHeaders.mat
load ../data/SuHeader.mat

% Set source/receiver coordinates
for k=1:length(SuTraceHeaders)
    s(k)=SuTraceHeaders(k).SourceX*1e-3;
    r(k)=SuTraceHeaders(k).GroupX*1e-3;
end

% Time vector
t = (0:SuHeader.ns-1)'*SuHeader.dt*1e-6;

% Set trace informations
h = (s - r)/2;
m = (s + r)/2;
nt = length(t);
nr = length(unique(r));
ns = length(unique(s));
Ir = find(unique(r)==2850);


%% Seismic image recovery with curvelets

% Sub-sampling mask
p = 2/3;    % subsampling ration
I = opDirac(nt);
idx = randperm(ns);
R1 = opRestriction(ns,idx(1:round(ns*p)));
R = kron(R1,I);

% Curvelet operator
C = opCurvelet(nt,ns,4,16,1,'WRAP');

% CS matrix
M = R*C';


%% Optimization

% Set optimization parameters
options.iterations = 200;
options.xtrue = C*CR(:);
x0 = [];
tau = 0;
sigma = 1e-8;

% SPGl1 solution
tic
[x1, res1, grad1, info1] = spgl1(M,R*CR(:),tau,sigma,x0,options);
t1=toc;

% Standard linearized Bregman with dynamic step length
options.method = 's';
options.step = 'd';
lambda = [];    % default lambda
tic
[x2, res2, info2] = lbm(M,R*CR(:),lambda,sigma,x0,options);
t2=toc;


%% Plot results

FigHandle1 = figure('Position', [400, 200, 1200, 1000]);
% Input data
subplot(1,3,1); 
imagesc(s(Ir),t,reshape(R'*R*CR(:),nt,ns),[-8e-1, 8e-1]); colormap(gray); title('Original Data')
xlabel('Source No.'); ylabel('TWT [s]');
% SPGl1 solution
subplot(1,3,2)
imagesc(s(Ir),t,reshape(C'*x1,nt,ns),[-8e-1, 8e-1]); colormap(gray);
title(['SPGL1 (SNR = ', ...
       num2str(SNR(CR,reshape(C'*x1,nt,ns)),'%4.2f'),', t = ',num2str(t1,'%4.2f'),' s)']);
   xlabel('Source No.'); ylabel('TWT [s]');
% Linearized Bregman solution
subplot(1,3,3)
imagesc(s(Ir),t,reshape(C'*x2,nt,ns),[-8e-1, 8e-1]); colormap(gray);
title(['Linearized Bregman (SNR = ', ...
       num2str(SNR(CR,reshape(C'*x2,nt,ns)),'%4.2f'),', t = ',num2str(t2,'%4.2f'),' s)']);
xlabel('Source No.'); ylabel('TWT [s]');

suptitle(['Trace interpolation, \lambda=',num2str(lambda,'%d'),...
    ' (30 % missing traces), ',num2str(options.iterations,'%d'),...
    ' iterations'])
   
% Residual plot
FigHandle2 = figure('Position', [400, 200, 1049, 895]);
iter = 1:1:options.iterations;
plot(iter,info1.rNorm2,iter,info2.residual)
axis([1 options.iterations 0 1.2*max(info1.rNorm2)])
legend('SPGl1','Linearized Bregman')
xlabel('Iteration No.')
ylabel('2-norm of residual')
title('Residual as a function of iteration number')


