function x = distVectorize(x)
% function y = distVectorize(x)
% 
% Vectorizes a distributed ND-array, provided that x is distributed over the last dimension
%
%
% Tim Lin, SLIM, EOS-UBC 2011
    
    if isvector(x)
        return
    end
    
    sizes = size(x);
    nd = length(sizes);
    
    if ~isa(x,'distributed') % if x is local then simply vectorize
        x = x(:);
    else
        spmd
            codistr = getCodistributor(x);
            assert(codistr.Dimension == nd, 'distVectorize: input must be distributed over the last dimension')
            
            partition_lastDim = codistr.Partition;
            partition_vec = partition_lastDim .* prod(sizes(1:end-1));
            
            new_codistr_vec = codistributor1d(1, partition_vec, [prod(sizes) 1]);

            local_x = getLocalPart(x);
            local_x = local_x(:);
            
            x = codistributed.build(local_x, new_codistr_vec,'noCommunication');
        end
    end
end
            
