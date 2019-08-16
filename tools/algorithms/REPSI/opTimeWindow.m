function op = opTimeWindow(A, window_start, window_end)
% OPTIMEWINDOW      Zeros out everything outside of a time window (1st dimension)
%
%    opTimeWindow(A, window_start, window_end) zeros out everything in nd-matrix A
%    outside of the time range (the 1st dimension) specified by window_start and window_end,
%    in terms of absolute samples, inclusive.

%   Copyright 2010, Tim Lin

sA = size(A);
nt = sA(1);
n = prod(sA);

% Handle distributed A using pSPOT operators. For distributed A this works on the last dimension.
if isa(A,'distributed')
    nt = sA(end);
    n = nt;
end

checkFor_validArrayIndex(window_start, nt)
checkFor_validArrayIndex(window_end, nt)
assert(window_end >= window_start, 'window_end must be greater than or equal to window_start')

subfunc_handle = @(x,mode) opTimeWindow_intrnl(x,mode);

op = opFunction(n,n,subfunc_handle); % return a SPOT operator using constructor opFunction

% if isa(A,'distributed') % return a pSPOT operator if A is distributed, that works on the last dimensions
%     op = oppKron2Lo(opFunction(n,n,subfunc_handle), opDirac(prod(sA(1:end-1))));
% end

function y = opTimeWindow_intrnl(x,mode)

    switch mode
        case 0
            y = {n,n,[0,0,0,0],{'opTimeWindow'}};

        case 1
            x = reshape(x,nt,[]);
            x(1:window_start-1,:) = 0;
            x(window_end+1:end,:) = 0;
            y = x(:);

        case 2
            x = reshape(x,nt,[]);
            x(1:window_start-1,:) = 0;
            x(window_end+1:end,:) = 0;
            y = x(:);
    end
end
    
end