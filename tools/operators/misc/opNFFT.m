function op = opNFFT(N, pos, ft_opt)
% OPNFFT 1D nonuniformly sampled fast Fourier transform
%
% Copyright 2008, Gilles Hennenfent
%  
% N      : number of samples (or seismic traces)
% pos    : positions or nodes \subest [-1/2, 1/2)
% ft_opt : transform options
%
% July, 2012 : Wrapped into a SPOT operator
%              Haneet Wason   


% options for NFFT
if nargin < 3
    ft_opt.method = 'gaussian';  % window function
    ft_opt.m      = 6;           % cut-off or Taylor-order     
    ft_opt.sigma  = 2;           % oversampling factor
end

parms = {N, pos, ft_opt};
fh    = @(x, mode) opNFFT_intrnl(parms{:}, x, mode);
op    = opFunction(length(pos), N, fh);

%------------------------------------------------------------

function y = opNFFT_intrnl(N, pos, ft_opt, x, mode)

if mode == 0
    y = {length(pos), N, [1,1,1,1], {'NFFT'}};

elseif mode == 1
    % Synthesis mode  
    y = -nfft(x, pos, N, ft_opt, 'notransp')/sqrt(N);

else
    % Analysis mode
    y = nfft(-x, pos, N, ft_opt, 'transp')/sqrt(N);
end
