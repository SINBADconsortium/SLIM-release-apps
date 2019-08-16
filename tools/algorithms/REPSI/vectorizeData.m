function x = vectorizeData(x)
% function y = vectorizeData(x)
%
%   Similar to MATLAB's vec(), but calls the more optimized "distVectorize" for
%   distributed x if it's distributed over the last dimension.
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
        end
        
        codistr = codistr{1}; % get the first lab
        
        if codistr.Dimension == nd     % distributed over the last dimension, call fast vec
            x = distVectorize(x);
        else
            x = x(:);   % requires the more general MATLAB vectorize
        end
    end
end
            
