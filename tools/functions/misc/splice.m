function y = splice(v,num_chunks)
% Splice - Split up a vector in to a fixed number of chunks
%
% Curt Da Silva, 2016
% 
% Usage:
%   y = splice(v,num_chunks);
%
    N = length(v);
    M = floor(N/num_chunks);
    R = N-M*num_chunks;
    
    y = cell(num_chunks,1);
    tot = 1;
    for i=1:num_chunks
        if i <= R
            y{i} = v(tot:tot+M);
        else
            y{i} = v(tot:tot+M-1);
        end
        tot = tot+length(y{i});
    end
end