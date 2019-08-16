function [f, gradf] = LSLinearMisfit(x,rhs, A)
% LSLINEARMISFIT - Least-squares misfit with optional masking, i.e.
%
%    f(x) = 0.5\|A(x) - rhs \|_2^2         (*)
%
%  or
%
%    f(x) = 0.5\|P_{A} x - rhs \|_2^2      (**)
%   
%
%    P_{A} x = x(i)  if  A(i) == 0
%            = 0     if  A(i) == 1
% 
%  Curt Da Silva
%  HTOpt v0.1
%  curtd@math.ubc.ca
%
% Usage:
%    [f, gradf] = LSLinearMisfit(x,rhs,{A})
%
% Input:
%    x     - input vector of length == length(rhs), can be
%            distributed or nondistributed
%    rhs   - data vector
%    A     - linear operator applied to x, (*) above
%            OR a logical vector with length == length(rhs) with 1s
%            where rhs is zero, (**) above
%
% Output:
%    f     - least-squares misfit value
%    gradf - Euclidean gradient

    isVec = exist('A','var') & not(isa(A,'opSpot') ) & (size(A,2)==1 || size(A,1)==1);
    distMode = isdistributed(rhs) || iscodistributed(rhs);
    if distMode
        isLogical = isVec;
    else
        isLogical = isVec & islogical(A);
    end
    if ~isLogical & isVec
        error('Non-logical vectors unsupported');
    end
    
    if (exist('A','var') == 0 )
        r = x - rhs;
        gradf = r;
    else
        if(isVec )
            if isLogical
                r = x;
                if distMode
                    spmd
                        codist = getCodistributor(r);
                        rloc = getLocalPart(r);
                        Aloc = getLocalPart(A);
                        rloc(Aloc) = 0;
                        r = codistributed.build(rloc,codist,'noCommunication');
                    end
                else
                    r(A) = 0;
                end
            else                
                r = A .* x;
            end
        else
            r = A *x;                        
        end
        if(isdistributed(x))
            clear x;
        end        
        r = r - rhs;
        if nargout >= 2
            if(isVec)
                if isLogical
                    gradf = r;
                else
                    gradf = (conj(A) .* r);
                end
            else
                gradf = A'*r;        
            end
        end
       
    end
         
    n = norm(r);
    f = 0.5*n^2; 
    
    if nargout >=2
        gradf = gradf;
    end
end







