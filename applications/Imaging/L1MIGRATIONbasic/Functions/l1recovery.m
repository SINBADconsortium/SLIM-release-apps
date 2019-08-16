function [x,rNorm,C,info] = l1recovery(K,b,nx,ny,mode,iter,x0,dy,SubScl)
% function x = l1recovery(K,b,nx,ny,mode)
% This function calculate the L1 recovery of Ax = b by using SPGL1
% 
% Inputs:
% 
%   K: opreator 
%   b: b = Ax
%   if Yogi's operator
%   nx: gridpoints of x in x dimension
%   ny: gridpoints of y in y dimension
%   
%   if Tim's operator, switch nx ny;
%   nx: gridpoints of y in y dimension
%   ny: gridpoints of x in x dimension
%
%   mode: if mode = 1; calculate in physical domain
%         if mode = 2; calculate in wavelet domain
%         if mode = 3; calculate in curvelet domain
%         if mode = 4; calculate in Shearlet domain
%         if mode = 5; calculate in curvelet domain with padding zeros
%         if mode = 6; calculate in curvelet domain by jump out when reaching
%                      the pareto curve.
%   iter: iterations used in SPGL1
%             
%  Author: Xiang Li
% 

if(nargin<9)
    SubScl = 2;
end

if(nargin<8)
    dy   = 15;
end
if(nargin<7)
    x0   = [];
end
if(nargin<6)
    iter = 500;
end
if(nargin<5)
    mode = 1;
end
% quit pareto
quitp = 0;
if mode == 8, quitp = 1;end


opts = spgSetParms('iterations' , iter , ...
    'verbosity' , 1 , ...
    'bpTol' , 1e-3 , ...
    'optTol' , 1e-3, ...
    'decTol', 5e-2, ...
    'quitPareto', quitp);
if     mode == 1
    
    [x] = spgl1(K,b,0,0,x0,opts);
    
elseif mode == 2
    W   = opWavelet(nx,ny);   
    KW  = K*W;
    [xw]= spgl1(KW,b,0,0,x0,opts);
    x   = W*xw;
    
elseif mode == 3
    C   = opCurvelet(ny,nx,4,16,0);
    KC  = K*C';
	[x,r,g,info]= spgl1(KC,b,norm(x0,1),0,x0,opts);
	iter_num    = info.iter 
   
elseif mode == 4
    S   = opShearlet(ny,nx);
    KS  = K*S;
    [xs]= spgl1(KS,b,0,0,x0,opts);
    x   = S*xs;

elseif mode == 5
    C   = opMeCurvelet2d(ny,nx,4,16,0,'ME');
	KC  = K*C;
    [xc]= spgl1(KC,b,0,0,x0,opts);
    x   = C*xc;
elseif mode == 6
    C   = opCurvelet(ny,nx,4,16,0,'ME')';
	KC  = K*C;
	[x,r,g,info]= spgl1(KC,b,norm(x0,1),0,x0,opts);
	rNorm    = info.rNorm2;
elseif mode == 7
    CV  = opCurvelet(ny,nx,4,16,0,'ME',0);
	S   = opCoreScale(ny,nx,4,16,0,SubScl,'ME');
	C   = CV'*S;
	KC  = K*C;
	[x,r,g,info]= spgl1(KC,b,norm(x0,1),0,x0,opts);
	iter_num    = info.iter
elseif mode == 8
    C   = opCurvelet2d(ny,nx,4,16,0,'ME');
	KC  = K*C';
	l1itr= 0
	x   = x0;
	rNorm = [];
	while l1itr < iter
	[x,r,g,info]= spgl1(KC,b,norm(x,1),0,x,opts);
	l1itr = l1itr + info.iter;
	rNorm    = [rNorm;info.rNorm2];
	end
	info.iter   = l1itr;
	info.rNorm2 = rNorm;
elseif mode == 9
%   C   = opMeCurvelet2d(ny,nx,4,16,0,'ME');

	P   = opWaExtend(ny,nx,'extend')
	m = 1;
	while 2^m < max(ny,nx)
		m = m + 1;
		n = 2^m;
	end
	C   = opWaveAtom(n,'ME');
	C   = C * P;
	
	KC  = K*C';
	[x,r,g,info]= spgl1(KC,b,norm(x0,1),0,x0,opts);
	rNorm    = info.rNorm2;
elseif mode == 10
	opts.quitPareto = 1;
	C   = opMeCurvelet2d(ny,nx,4,16,0,'ME');
	KC  = K*C;
	[x,r,g,info]= spgl1(KC,b,norm(x0,1),0,x0,opts);
	rNorm    = info.rNorm2;
elseif mode == 11
	opts.quitPareto = 1;
	C   = opMeCurvelet2d(ny,nx,4,16,0,'ME');
	KC  = K*C;
	l1itr= 0
	x   = x0;
	rNorm = [];
	for m = 1:2
		[x,r,g,info]= spgl1(KC,b,norm(x,1),0,x,opts);
		l1itr = l1itr + info.iter;
		rNorm    = [rNorm;info.rNorm2];
	end
	info.iter   = l1itr;
	info.rNorm2 = rNorm;
elseif mode == 12
    C   = opMeCurvelet2d(ny,nx,4,16,0,'ME');
	KC  = K*C;
	[x,r,g,info]= spgl2(KC,b,norm(x0,1),0,x0,opts);
	rNorm    = info.rNorm2;
elseif mode == 13
	clear opts
	opts = spgSetParms('iterations' , iter);
        C   = opMeCurvelet2d(ny,nx,4,16,0,'ME');
	KC  = K*C;
	% [x,r,g,info]= spgl1(KC,b,norm(x0,1),0,x0,opts);
	% rNorm    = info.rNorm2;
	if x0 == [],x0 = zeros(size(C,2),1);end
	corrections = 1;
	 [x,~,~,info] = pqnl1_2(KC,b, 0, 1e-3, x0, opts,0,corrections);
	rNorm =[];
	% l1 migration
	% opts.iterations = 100;
	% opts.fid = fopen('log_freq_spg.txt', 'w'); 
	% [x_spg,r_spg,g_spg,info_spg] = spgl1(A, b, 0, 1e-3, zeros(size(A,2),1), opts);
	% m_spg = C'*x_spg;
	% % opts.fid = fopen('log_freq_pqn.txt', 'w'); 
	% sigma_ref = info_spg.rNorm;
	% corrections_list = 1;
	% for i = 1:length(corrections_list)
	%     corrections = corrections_list(i);
	%     [x_pqn(:,i),~,~,info_pqn(:,i)] = pqnl1_2(A, b, 0, 1e-3, zeros(size(A,2),1), opts,sigma_ref,corrections);
	%     m_pqn(:,i) = C'*x_pqn(:,i);
	%     nz = length(z);
	%     figure(i); subplot(2,1,1);imagesc(reshape(real(m_spg),nz,nr));
	%     subplot(2,1,2);imagesc(reshape(real(m_pqn(:,i)),nz,nr));
	%     save info info_spg info_pqn m_spg m_pqn
	% end
	
end




% P   = opWaExtend(opts.ny,opts.nx,'extend')
% C   = opWaveAtom(n,'ME');

