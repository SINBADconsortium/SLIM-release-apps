function y = mtimes(A,B)
%*   Product of two operators.
%
%   A*B  returns an operator that is the product of two operators.
%
%   See also opFoG.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.
   
%   http://www.cs.ubc.ca/labs/scl/spot

% Note that either A or B must be belong to the opSpot class because this
% function gets called for both M*C, C*M, where C is the class and M is a
% matrix or vector. This gives the following options, with s for scalar and
% C for any instance of an opSpot:
%
% 1) M*C, implemented as (C'*M')'
% 2) C*M
% 3) s*C
% 4) C*s
% 5) C*C, either of which can be a foreign class

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SeisDataContainer Proprocessing
% Always strip the data at the highest level possible(mtimes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isa(B,'SeisDataContainer') && isa(A,'opSpot')
    y = dataMultiply(B,A); % Let the datacontainer itself handle the stripping
    return;
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mode: Explicit opFog Saver
%
% This is to save memory and prevent unnecessary multiplies in opFog
% Only applies to in-core serial explicit matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (isnumeric(A) || isa(A, 'opMatrix')) && isa(B, 'opMatrix')

    y = opMatrix(double(A)*double(B));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mode 1: M*C
% Mode 3: s*C - Here we also handle the special case where C is 1-by-M.
%               If so, then we recast this as (C'*s)', which results in
%               a call to the "usual" matrix-vector product.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~isa(A,'opSpot') % A is not spot
    if isscalar(A) && (B.m ~= 1)
       % s*C (mode 3)
       if isa(A,'oppSpot') || isa(B,'oppSpot')
           y = oppFoG(A,B);
       else
           y = opFoG(A,B);
       end
    else
       % M*C (mode 1)
       y = (B' * A')';
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mode 2: C*M
% Mode 4: C*s - Here we also handle the special case where C is N-by-1.
%               If so, then we recast this as (C'*s)', which results in
%               a call to the "usual" matrix-vector product.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~isa(B,'opSpot') % A is spot, B isnt
        p = size(B,1);
        
        if isscalar(B) && A.n ~= 1
            % C*s (mode 4)
            if isa(A,'oppSpot') || isa(B,'oppSpot')
                y = oppFoG(A,B);
            else
                y = opFoG(A,B);
            end
        
        elseif A.n ~= p && ~isscalar(A)
        % Raise an error when the matrices do not commute. We make an
        % exception for 1-by-1 operators.
            sizB = ['Matrix[' num2str(size(B)) ']'];
            error('Matrix dimensions must agree when multiplying by %s and %s.',...
            char(A), sizB);
        
        else % Perform operator*matrix
        % Unsupported data types warning
        if ~isnumeric(B) && ~isa(A,'oppSpot') && ~islogical(B)
            warning(['Data type "' class(B)...
            '" is not officially supported by spot, proceed at own risk']);
        end
        
            if isempty(A)
                y = zeros(A.m,size(B,2), class(B));
            elseif isempty(B)
                y = zeros(A.m,0, class(B));
            else
                y = applyMultiply(A,B,1);
            end
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Both args are Spot ops.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    if isa(A,'oppSpot') || isa(B,'oppSpot')
        y = oppFoG(A,B);
    else
        y = opFoG(A,B);
    end
end
