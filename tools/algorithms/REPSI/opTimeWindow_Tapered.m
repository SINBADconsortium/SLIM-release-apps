function op = opTimeWindow_Tapered(A, window_start, window_end, taperlength)
% OPTIMEWINDOW_TAPERED      Zeros out everything outside of a time window (1st dimension),
%                           with a smooth cosine taper window inside the time window
%
%    opTimeWindow_Tapered(A, window_start, window_end) zeros out everything in nd-matrix A
%    outside of the time range (the 1st dimension) specified by window_start and window_end,
%    in terms of absolute samples, inclusive.
%    
%    Inside the time range, the signal gets weighted by a cosine taper on each side by
%    taperlength samples. This is meant to ensure smoothness in the time domain and decrease
%    side-lobing effects in the frequency spectrum.

%   Copyright 2010, Tim Lin

sA = size(A);
nt = sA(1);
n = prod(sA);

% Handle distributed A using pSPOT operators. For distributed A this works on the last dimension.
if isa(A,'distributed')
    nt = sA(end);
    n = nt;
end

% Validate window specs
checkFor_validArrayIndex(window_start, nt)
checkFor_validArrayIndex(window_end, nt)
assert(window_end >= window_start, 'window_end must be greater than or equal to window_start')

% Validate taper specs
window_length = window_end - window_start + 1;
checkFor_validArrayIndex(taperlength*2, window_length)

taper=0.5*(1-cos((1:taperlength*2-1)*pi/taperlength)');

subfunc_handle = @(x,mode) opTimeWindow_intrnl(x,mode);

op = opFunction(n,n,subfunc_handle); % return a SPOT operator using constructor opFunction

if isa(A,'distributed') % return a pSPOT operator if A is distributed, that works on the last dimensions
    op = oppKron2Lo(op, opDirac(prod(sA(1:end-1))));
end

    function y = opTimeWindow_intrnl(x,mode)

        switch mode
            case 0
                y = {n,n,[0,0,0,0],{'opTimeWindow'}};

            case 1
                x = reshape(x,nt,[]);
                
                x(1:window_start-1,:) = 0;
                x(window_end+1:end,:) = 0;
                
                x(window_start:window_start+taperlength-1,:) = x(window_start:window_start+taperlength-1,:).*repmat(taper(1:taperlength),[1 size(x,2)]);
                x(window_end:-1:window_end-taperlength+1,:) = x(window_end:-1:window_end-taperlength+1,:).*repmat(taper(1:taperlength),[1 size(x,2)]);
                
                y = x(:);

            case 2
                x = reshape(x,nt,[]);
                
                x(1:window_start-1,:) = 0;
                x(window_end+1:end,:) = 0;
                
                x(window_start:window_start+taperlength-1,:) = x(window_start:window_start+taperlength-1,:).*repmat(taper(1:taperlength),[1 size(x,2)]);
                x(window_end:-1:window_end-taperlength+1,:) = x(window_end:-1:window_end-taperlength+1,:).*repmat(taper(1:taperlength),[1 size(x,2)]);
                
                y = x(:);
                
        end
    end
    
end