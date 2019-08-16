function Y = ttm(X,A,dims)
%TTM - Multilinear product of a tensor with matrices/SPOT operators
%  
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
% 
%
%  Usage:
%    Y = TTM(X,A,N) computes the N-mode matrix product Y = A o_N X
%    for a tensor X and a matrix/SPOT operator A that satisfies 
%      
%    Y = TTM(X,A,DIMS) for a cell array A containing matrices/SPOT
%    operators A{i}, computes the N-mode products 
%    Y = A{1} o_DIMS(1) A{2} o_DIMS(2) o_... A{length(A)} o_DIMS(length(A)) X
%   
%    If A{i} is empty, this is treated as an identity operator in dimension i
%  
%    Y = TTM(X,A) for a cell array A is the same as Y = TTM(X,A, 1:length(A))
%
%  
if nargin < 2
    error('Need at least 2 arguments');
end

if ~isfloat(X) 
    error('First argument must be a MATLAB multidimensional array');
end

if ~isa(A,'cell')
    if ~isfloat(A) && ~isa(A,'opSpot')
        error(['Second argument must be either a cell array, a matrix, or a SPOT operator']);
    end
    A = {A};
else
    if ~all(cellfun( @(x)( isfloat(x) || isa(x,'opSpot') ), A))
        error(['All elements of A must be either matrices or SPOT operators']);
    end
end

if ~exist('dims','var')
    dims = 1:length(A);
end

B = cell(max(length(size(X)),max(dims)),1);
for i=1:length(B)
    if isempty(find( i == dims, 1 ))
        B{i} = opDirac(size(X,i)); 
    else
        if isempty(A{find(i == dims,1)})
            B{i} = opDirac(size(X,i)); 
        else
            B{i} = A{find( i == dims,1)};
        end
    end
end

newdims = zeros(1,length(dims));
for i=1:length(B)
    newdims(i) = size(B{i},1);
end

kronstr = '';
for i=length(B):-1:1
   kronstr = [kronstr, 'B{' num2str(i) '},']; 
end
kronstr = kronstr(1:end-1);
eval(['P = opKron(' kronstr ');']);


Y = reshape(P * vec(X), newdims);



end