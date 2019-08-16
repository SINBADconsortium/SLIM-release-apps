function [Rtrain, Rtest, Jtrain, Jtest] = ndimSubsampling(dims, subDims, subPercent, subType, parallel )
% NDIMSUBSAMPLING - Generates SPOT operators for subsampling multidimensional arrays (i.e. the training set) + comparing interpolation results on the subsampled indices ( i.e. the test set ).
%
% Usage:
%   [Rtrain, Rtest, Jtrain, Jtest] = ndimSubsampling(dims,subDims,subPercent,{subType}, {parallel});
%  
% Input:
%   dims        - size of dimensions
%   subDims     - subset of 1:length(dims) corresponding to the dimensions to be subsampled
%   subPercent  - ratio of points to remove along subDims, between 0 and 1
%   subType     - 0 : remove points (default)
%                 1 : Gaussian subsampling
%   parallel    - if true, will use pSPOT instead of SPOT (default: false)
%
% Output:
%   Rtrain      - SPOT operator which restricts to training set
%   Rtest       - SPOT operator which restricts to test set
%   Jtrain      - random indices used for the training set (when subType == 0), empty otherwise
%   Jtest       - random indices used for the test set (when subType == 0), empty otherwise

POINTS = 0;
GAUSSIAN = 1;

if nargin < 3
    error('At least 3 arguments required');
end

assert( 0 < subPercent && subPercent < 1, 'subPercent must be between 0 and 1' );

if exist('subType','var')==0
    subType = POINTS;
end
if exist('parallel','var')==0    
    parallel = false;
end

nRemovePoints = round( prod(dims(subDims))*subPercent );
keepDims = setdiff(1:length(dims),subDims);
nsubDims = prod(dims(subDims));
nkeepDims = prod(dims(keepDims));

switch subType
  case POINTS   
    Jtrain = sort(randperm( nsubDims , nsubDims - nRemovePoints),'ascend');
    Jtest = sort(setdiff(1:nsubDims, Jtrain ),'ascend');    
    Strain = opRestriction_swp( nsubDims, Jtrain );
    Stest = opRestriction_swp( nsubDims, Jtest );
  case GAUSSIAN
    Jtrain = [];
    Jtest = [];
    Strain = opGaussian(nsubDims - nRemovePoints, nsubDims ) / sqrt( nsubDims - nRemovePoints );
    Stest = opDirac(nsubDims);
  otherwise
    error('Unsupported subsampling type');
end

% subsample all dimensions
if length(subDims) == length(dims)
    Rtrain = Strain;
    Rtest = Stest;    
% subsample first max(subDims) dimensions 
elseif ismember(1,subDims) && max(subDims) == length(subDims)
    if ~parallel
        Rtrain = opKron( opDirac(nkeepDims), Strain );
        Rtest = opKron( opDirac(nkeepDims), Stest );
    else
        Rtrain = oppKron2Lo( opDirac(nkeepDims), Strain );
        Rtest = oppKron2Lo( opDirac(nkeepDims), Stest );
    end
    
% subsample last min(subDims) dimensions
elseif ismember(length(dims),subDims) && ((max(subDims) - min(subDims) + 1) == length(subDims))
    if ~parallel
        Rtrain = opKron( Strain, opDirac(nkeepDims) );
        Rtest = opKron( Stest, opDirac(nkeepDims) );
    else
        Rtrain = oppKron2Lo( Strain, opDirac(nkeepDims) );
        Rtest = oppKron2Lo( Stest, opDirac(nkeepDims) );
    end

% need to permute before subsampling
else 
    P = opPermute(dims, [subDims, keepDims] );
    if ~parallel
        Rtrain = opKron( opDirac(nkeepDims), Strain ) * P;
        Rtest = opKron( opDirac(nkeepDims), Stest ) * P;
    else
        Rtrain = oppKron2Lo( opDirac(nkeepDims), Strain ) * P;
        Rtest = oppKron2Lo( opDirac(nkeepDims), Stest ) * P;
    end
    
end
    
    
