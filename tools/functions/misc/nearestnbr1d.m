function I = nearestnbr1d(xin,xout)
% 1D Nearest neighbour algorithm. For each point in xout, 
% finds the index of the closest point in xin. Uses a 
% sorting + bracketing system to avoid unnecessary distance
% computations/memory overhead.
%
% Curt Da Silva, 2015
%
% Usage:
%   I = nearestnbr1d(xin,xout);
%
% Input:
%   xin     - input grid
%   xout    - output grid
%  
% Output:
%   I       - indices mapping output grid pt -> closest input grid pt
%
    assert(min(xin) <= min(xout) && max(xin) >= max(xout), 'Output grid must be contained within input grid');
    % Sort xin, xout
    [xin,Jin] = sort(xin,'ascend');
    [xout,Jout] = sort(xout,'ascend');

    xend = xin(2);
    istart = 1;
    I = zeros(length(xout),1);
    for i=1:length(xout)
        % Keep incrementing until the current output point
        % is in the interval [xin(istart),xin(istart+1)]
        while xout(i) > xend
            istart = istart+1;
            if istart==length(xin),
                xend = inf;
            else
                xend = xin(istart+1);
            end
        end
        % Find out whether the current point is closer to the
        % start or the end of the interval
        if abs(xout(i)-xin(istart)) <= abs(xout(i)-xend)
            I(i) = istart;
        else
            I(i) = istart+1;
        end
    end
    
    % Convert results back to unsorted xin/xout 
    Jin_inv(1:length(Jin)) = Jin;
    Jout_inv(1:length(Jout)) = Jout;    
    
    I = Jin_inv(I);
    I(Jout_inv) = I;
    
end