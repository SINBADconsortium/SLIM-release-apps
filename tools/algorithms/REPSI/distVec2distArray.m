function x = distVec2distArray(x, sizes)
% function y = distVec2distCube(x, sizes)
% 
% x         is a distributed vector
% sizes     is an array of sizes in different dimensions
% 
% Takes in distributed vector x of length prod(sizes) and do a "reshape"
% to a N-D array of size "sizes" and is distributed over the last dimension
%
% example:
%   x = distributed.randn(1000,1);
%   x = distVec2distArray(x, [10 10 10]);
% 
% x is now a 10x10x10 array that is distributed over the last (3rd) dimension
% using the default partitioning scheme.
%
% Tim Lin, SLIM, EOS-UBC 2011
    
    if not(isvector(x))
        error('distVec2distArray: Input x needs to be a column vector')
    end
    
    sizes = sizes(:).';
    assert(length(x) == prod(sizes), 'Length of input must equal product of sizes')
    
    if ~isa(x,'distributed') % if x is local then simply distribute correctly
        x = reshape(x, sizes);
        x = distributed(x);
        
    else
        % reshape according to sizes, distribute over last dimension
        size_lastDim = sizes(end);
        lastDim = length(sizes);
        
        spmd
            % checks if redistribution is necessary, if the size of each local part of x is evenly 
            % divisible by prod(sizes(1:end-1)) then redistribution is not necessary
            codistr = getCodistributor(x);
            partition_vec = codistr.Partition;
            
            if all(mod(partition_vec,prod(sizes(1:end-1)))) == 0
                
                partition_lastDim = partition_vec ./ prod(sizes(1:end-1));
                new_codistr = codistributor1d(lastDim, partition_lastDim, sizes);
                
                local_x = getLocalPart(x);
                
                sizeLocal_lastDim = length(local_x) / prod(sizes(1:end-1));
                sizeLocal = [sizes(1:end-1) sizeLocal_lastDim];
                local_x = reshape(local_x, sizeLocal);
                x = codistributed.build(local_x, new_codistr,'noCommunication');
                
            else
                % probe the behaviour of default distribution scheme using a test vector
                test_vec = codistributed.zeros(size_lastDim,1);
                
                codistr = getCodistributor(test_vec);
                partition_lastDim = codistr.Partition;
                partition_vec = partition_lastDim .* prod(sizes(1:end-1));
            
                new_codistr_vec = codistributor1d(1, partition_vec, [prod(sizes) 1]);
                new_codistr = codistributor1d(lastDim, partition_lastDim, sizes);
            
                % redistribute as a vec, then reshape locally, then build as a 3d array at the end
                x = redistribute(x, new_codistr_vec);
                local_x = getLocalPart(x);
            
                sizeLocal_lastDim = length(local_x) / prod(sizes(1:end-1));
                sizeLocal = [sizes(1:end-1) sizeLocal_lastDim];
                local_x = reshape(local_x, sizeLocal);
                x = codistributed.build(local_x, new_codistr,'noCommunication');
            end
        end
    end
end
            
