function indices = index_block(idx,block_size)
% Chunks up an array in to pieces of a specified size
% 
% Curt Da Silva, 2016
% 
% Usage:
%   indices = index_block(idx,block_size);
%
% Input:
%   idx        - N x 1 vector
%   block_size - size of each block
% 
% Output:
%   indices    - cell array whose union is idx, each element is a block_size x 1 
%                vector
    N = length(idx);
    if N <= block_size
        indices = cell(1,1); indices{1} = idx;
    else
        kup = floor(N/block_size);
        indices = cell(kup,1);
        total = 1;
        for j=1:kup
            indices{j} = idx(total:total+block_size-1);
            total = total + block_size;
        end
        if kup*block_size < N
            indices{end+1} = idx(total:end);
        end
    end
end