function op = opFDCTcore(m, n, nbscales, nbangles)
%
%   Copyright 2008, Gilles Hennenfent
% 
% July, 2012 : Wrapped into a SPOT operator
%              Haneet Wason


if nargin < 3, nbscales = ceil(log2(min(m,n)) - 3); end;
if nargin < 4, nbangles = 16;                       end;

% Compute length of curvelet coefficient vector
C = fdct_wrapping_core(randn(m,n),nbscales,nbangles);

hdr{1}{1} = size(C{1}{1});
cn = prod(hdr{1}{1});  
for i = 2:nbscales
    nw = length(C{i});
    for j = 1:nw
        hdr{i}{j} = size(C{i}{j});
        cn = cn + prod(hdr{i}{j});
    end
end

parms = {m, n, cn, hdr, nbscales, nbangles};
fh = @(x,mode) opFDCTcore_intrnl(parms{:}, x, mode);
op = opFunction(m*n, cn, fh);


function y = opFDCTcore_intrnl(m, n, cn, hdr, nbs, nba, x, mode)

if mode == 0
    y = {m*n, cn, [1,1,1,1], {'FDCTcore'}};
elseif mode == 1
    % Synthesis mode  
    x = fdct_v2c_mat(x,hdr);
    y = ifdct_wrapping_core(x,m,n);
    y = y(:);
else
    % Analysis mode
    y = fdct_wrapping_core(reshape(x,m,n),nbs,nba);
    y = fdct_c2v_mat(y,cn);
end
