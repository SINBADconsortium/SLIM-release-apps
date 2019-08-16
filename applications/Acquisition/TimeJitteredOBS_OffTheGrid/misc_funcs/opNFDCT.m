function op = opNFDCT(m, n, pos, nbscales, nbangles, ft_opt)
%
%   Copyright 2008, Gilles Hennenfent
%
% July, 2012 : Wrapped sparco operators into SPOT operators
%              Haneet Wason

pos = pos(:);
if nargin < 4, nbscales = ceil(log2(min(m,n)) - 3); end;
if nargin < 5, nbangles = 16;                       end;
if nargin < 6
    ft_opt.method = 'gaussian';
    ft_opt.m      = 6;
    ft_opt.sigma  = 2;
end

T = opFDCTcore(m, n, nbscales, nbangles);
F = opFFT1C(m);
N = opNFFT(n, pos, ft_opt);

parms = {m, n, T, F, N, pos};
fh    = @(x,mode) opNFDCT_intrnl(parms{:},x,mode);
op    = opFunction(length(pos)*m, size(T,2), fh);

%-------------------------------------------------------------

function y = opNFDCT_intrnl(m, n, T, F, N, pos, x, mode)

if mode == 0
%    infos = T([],0);
%    y = {length(pos)*m,infos{2},[1,1,1,1],{'NFDCT'}};
     y = {length(pos)*m, size(T,2), [1,1,1,1], {'NFDCT'}};
elseif mode == 1
    % Synthesis mode  
%    tmp = reshape(T(x,1),m,n);
    tmp = reshape(T*x, m, n);
    y   = zeros(m, length(pos));
    
    %%% trace-by-trace
    for i = 1:n
%        tmp(:,i) = conj(F(tmp(:,i),1));
        tmp(:,i) = conj(F*tmp(:,i));
    end
    
    %%% time-by-time
    for i = 1:m
%        y(i,:) = -conj(N(tmp(i,:).',1));
        y(i,:) = -conj(N*tmp(i,:).');
    end
    
    y = y(:);
else
    % Analysis mode
    tmp  = reshape(x, m, length(pos));
    tmp1 = zeros(m, n);
    
    %%% time-by-time
    for i = 1:m
%        tmp1(i,:) = N(-tmp(i,:)',2);
        tmp1(i,:) = N'*(-tmp(i,:)');
    end
    
    %%% trace-by-trace
    for i = 1:n
%        tmp1(:,i) = F(conj(tmp1(:,i)),2);
        tmp1(:,i) = F'*conj(tmp1(:,i));
    end
    
%    y = T(tmp1(:),2);
    y = T'*tmp1(:);
end
