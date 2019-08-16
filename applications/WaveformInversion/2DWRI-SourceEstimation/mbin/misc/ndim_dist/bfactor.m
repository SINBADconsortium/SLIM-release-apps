function y = bfactor(x,d)
% BFACTOR - (approximate) balanced integer factorization
% 
%   Writes an integer I as
%    I = q1 * q2 * ... * qd
%   where max(qi) - min(qi) is as small as possible
%   
%   Curt Da Silva
%   curtd@math.ubc.ca
%   Feb 2015
%
%   Input: 
%     x - positive integer
%     d - number of factors
%   
%   Output:
%     y - [q1,q2,...,qd] approximately balanced factors

z = [0,inf];

for i=1:10
    
    y = factor(x);
    
    if length(y) < d
        error('input has less than d factors');
    end
    while length(y) > d
        y = sort(y);
        % Half of the time, randomly swap the 2nd smallest factor with the 3rd smallest factor
        if rand >= 0.5 && length(unique(y)) > 2, 
            [yu,ia,ib] = unique(y);
            t = yu(ib(2)); 
            y(ia(2)) = yu(ib(1)); 
            y(ia(1)) = t;
        
        end
        y = [y(1)*y(2),y(3:end)];
    end
    if max(y)-min(y) <= max(z)-min(z)
        z = y;
    end
end
y = z;
end