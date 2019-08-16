function op = opFDCT(m, n, nbscales, nbangles)
%
%   Copyright 2008, Gilles Hennenfent
% 
% July, 2012 : Wrapped all sparco operators into SPOT operators
%

if nargin < 3, nbscales = ceil(log2(min(m,n)) - 3); end;
if nargin < 4, nbangles = 16;                       end;

T = opFDCTcore(m, n, nbscales, nbangles);
Ft = opFFT2C(m, n);
%Ftspot = opTranspose(opFunction(m*n, m*n, Ft));
Ftspot = opFunction(m*n, m*n, Ft);
op = Ftspot'*T;

%op = @(x, mode) opFDCT_intrnl(A, x, mode);


%function y = opFDCT_intrnl(A, x, mode)

%if mode == 0
%    infos = A([],0);
%    y = {infos{1}, infos{2}, infos{3}, {'FDCT'}};
%elseif mode == 1
%    % Synthesis mode  
%    y = A(x,1);
%else
%    % Analysis mode
%    y = A(x,2);
%end
