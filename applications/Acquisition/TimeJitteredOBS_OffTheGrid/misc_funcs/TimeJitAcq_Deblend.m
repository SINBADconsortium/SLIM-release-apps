function Dest = TimeJitAcq_Deblend(label, jitacq, nt, nr, ns, ds, RM, jitD, opt, options)

% Setup measurement operator [A]
% The operator A is a combination of the sampling operator RM and the sparsifying operator S.

% Non-equispaced positions for opNFDCT
pos = sort([jitacq.sjitb1arr1 jitacq.sjitb1arr2]);
posNC = pos/(ns*ds) - 0.5;
fname_posNC = [label '_posNC.rsf'];
rsf_write_all(fname_posNC, {'out=stdout'}, posNC)

% Non-equispaced curvelet transform
% NOTE - mode 1 ---> synthesis (NC), mode 2 ---> analysis (NC')
NC = opNFDCT(nr, ns, posNC);

% opSplineWavelet(U, V, filt, smooth, levels)
% NOTE - mode 1 ---> synthesis (W), mode 2 ---> analysis (W')
W = opSplineWaveletSPOT(nt, 1, nt, 3, 5);                   

% oppKron2Lo: kronecker tensor product to act on a distributed vector
S = oppKron2Lo(NC', W', 1);

% Measurement operator
A = RM*S'; 


% Setup recovery operator [Sreg]

% Regular curvelet transform
C = opFDCT(nr, ns, max(1,ceil(log2(min(nr,ns)) - 3)), 16);

% Recovery operator
Sreg = oppKron2Lo(C', W', 1);                              


% Solve the one-norm recovery problem
% spgl1: (A, b, tau, sigma, x, options)
% xest: estimated (synthesis) curvelet coefficients
xest = spgl1(A, jitD(:), opt.spgl1_tau, opt.spgl1_sigma, opt.spgl1_x, options);

% Recover data
Dest = reshape(Sreg'*xest, [nt nr ns]);

end  % function

