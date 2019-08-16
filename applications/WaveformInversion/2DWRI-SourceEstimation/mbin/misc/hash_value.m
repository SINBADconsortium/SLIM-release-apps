function y = hash_value(x)
% HASH_VECTOR - Computes a SHA-256 hash of an array, string, or cell array
% of arrays and strings. 
%
% Author: Curt Da Silva
% 
%
% Usage: 
% 
%  y = hash_value(x);
% 
% Input:
%  x - numerical array, string, or cell array of strings + numerical arrays
% 
% Output:
%  y - SHA-256 hash string
% 
import java.security.MessageDigest;
import java.math.*;    
    function z = hash(w)      
        if isempty(w), z = ''; return; end
        if isfloat(w), w = num2str(w,15); w = w(:); end
        mD = java.security.MessageDigest.getInstance('SHA-256');
        java_string = java.lang.String(char(w));
        mD.update(java_string.getBytes());
        z = mD.digest();
        bi = java.math.BigInteger(1, z);
        z = char(bi.toString(16));
    end

L = 1e5;
if iscell(x)
    y = hash_value(x{1});
    for i=2:length(x)
        y = hash([y, hash_value(x{i})]);
    end    
elseif isfloat(x)    
    if issparse(x), [i,j,s] = find(x); x = [i,j,s]; end    
    x = x(:);        
    if length(x) < L      
        w = num2str(x,15); w = w(:);
        y = hash(w);
    else        
        y = hash(x(1:L));
        nel = L+1;
        for i=2:floor(length(x)/L)           
           y = hash([y,hash(x(nel:nel+L-1))]);
           nel = nel + L;
        end
        y = hash([y,hash(x(nel:end))]);
    end    
elseif ischar(x)
    y = hash(x);
end


end

