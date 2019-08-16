function y = filter_bank(n,dx,filter_name,cutoff_freq, opts)   
% FILTER_BANK - Implementation of a low-pass filter, based on
%   http://www.labbookpages.co.uk/audio/firWindowing.html
%
% Curt Da Silva, 2015
%
% Usage:
%   y = filter_bank(n,dx,filter_name,cutoff_freq,opts);
%
% Input:
%   n           - signal length
%   dx          - signal grid spacing
%   filter_name - one of 
%                   'hamming' - hamming windowed low pass filter
%                   'kaiser' - kaiser windowed low pass filter
%   cutoff_freq - cutoff frequency for the low-pass
%   opts        - struct with options
%       .full        - if true, zero pads the filter to have length
%                      n (default: false)
%       .ripple      - for the 'kaiser' filter, the amount of ripple
%                      around 1 tolerated in the frequency domain
%                      (default: 0.05)
%       .trans_width - for the 'kaiser' filter, transition width of
%                      the passband -> stopband, in Hz (default:
%                      0.1*sampling frequency)
%   
% Output:
%   y - filter
%
    if isfield(opts,'full'),full = opts.full; else full = false; end

    [x,f] = spacefreq_grid(n,dx);
    samp_freq = 2*max(f);
    
    ft = cutoff_freq/samp_freq;
    if strcmp(filter_name,'kaiser')
        if ~isfield(opts,'ripple'),opts.ripple = 0.05; end
        if ~isfield(opts,'trans_width'),opts.trans_width = 0.1*samp_freq; end
        A = -20*log10(opts.ripple);
        
        tw = 2*pi*opts.trans_width/samp_freq;
        if A > 21
            M = ceil((A-7.95)/(2.285*tw));
            if A < 50
                beta = 0.5842*(A-21)^(0.4)+0.07886*(A-21);
            else
                beta = 0.1102*(A-8.7);
            end
        else
            M = ceil(5.79/tw);
            beta = 0;
        end        
    else    
        M = opts.filter_length-1;
    end
    ix = vec(0:M) - M/2;
    y = 2*ft*ones(M+1,1);
    I = ix~=0;
    y(I) = sin(2*pi*ft*ix(I))./(pi*ix(I));    
    w = ones(size(y));
    switch filter_name
      case 'hamming'
        ix = vec(0:filter_length-1);
        w = 0.54 - 0.46*cos(2*pi*ix/M);
      case 'kaiser'
        w = kaiser_window(ix,M/2,beta);
    end
    y = y.*w;
    M = M+1;
    % zero padding
    if full
        if mod(n-M,2)==0
            y = [zeros((n-M)/2,1); y; zeros((n-M)/2,1)];
        else
            y = [zeros((n-M-1)/2,1);y; zeros((n-M+1)/2,1)];
        end
    end
end