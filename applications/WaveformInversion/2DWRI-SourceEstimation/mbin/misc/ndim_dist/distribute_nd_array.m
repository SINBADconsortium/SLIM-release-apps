function [A,loc_size,start_indices,end_indices] = distribute_nd_array(A, dims, ndist_dims,factors)              
% DISTRIBUTE_ND_ARRAY - Rearranges a distributed, n-dimensional array to ensure that 
% local parts of the array correspond to hypercubes of the underlying array. This effectively generalizes 
% codistributed2d to multiple dimensions.
%
% Curt Da Silva
% March 2015
% curtd@math.ubc.ca
%
% Usage:
%   [A,loc_size,start_indices,end_indices] = distribute_nd_array(A,dims,ndist_dims,{factors});
%  
% Input:
%   A             - distributed nd array, vectorized
%   dims          - dimensions of the underlying array
%   ndist_dims    - number of dimensions to distribute. Distributed dimensions will always be the last dimensions of dims
%   factors       - vector with prod(factors)==parpool_size(), specifying the size reduction in each dimension
%                   (default): balanced factorization
%                   note: number of open matlab pools must have a prime factorization with at LEAST length(ndist_dims) factors, otherwise
%
% Output:
%   A             - distributed nd array, vectorized
%   loc_size      - composite, size of the hypercube at the current worker
%   start_indices - composite, starting indices for each dimension of the local array
%   end_indices   - composite, ending indices for each dimension of the local array
%
% Example:
%   dims = [10,10,10,10];
%   X = pSPOT.utils.distVectorize(distributed.randn(prod(dims),1));
%   
%   % distributed over the last 3 dimensions                                               
%   [X,loc_size,start_indices,end_indices] = distributed_nd_array(X,dims,3);
%
%   Xg = reshape(gather(X),dims);
%   spmd,
%     Xloc = reshape(getLocalPart(X),loc_size);  %Each local part is a hypercube
%     i1 = 1:10;                                 %All non-distributed dimensions are always present
%     i2 = start_indices(1):end_indices(1);      %Indices of each distributed dimension
%     i3 = start_indices(2):end_indices(2);
%     i4 = start_indices(3):end_indices(3);
%     assert(norm(vec(Xloc) - vec(Xg(i1,i2,i3,i4)))==0); 
%   end
%


nlabs = parpool_size();
assert(nlabs > 0, 'Need open matlab workers');
assert(ndist_dims <= length(dims), 'Cant distribute over more dimensions than there are');

% Distributed dimensions are always the last ones, by convention
dist_dims = dims(end-ndist_dims+1:end);
non_dist_dims = dims(1:length(dims)-ndist_dims);
d = length(dist_dims);

if d == 1,  
    A = reshape(A,prod(dims(1:end-1)),dims(end));
    spmd,
        A = redistribute(A,codistributor1d(2, [],[prod(dims(1:end-1)),dims(end)]));
        loc_size = [non_dist_dims, size(getLocalPart(A),2)];
        [start_indices,end_indices] = globalIndices(A,2);
    end
    return;
end

% Balanced factors of #labs
if exist('factors','var')==0 || isempty(factors) || prod(factors)~= nlabs || min(factors) < 2
    %Sort factors so that the largest dimensions get the biggest reduction in size
    if length(factor(nlabs)) < d, 
        error('Number of workers must have at least #(distributed dims) prime factors ');
    end
    if exist('factors','var')==0 || isempty(factors) 
        factors = sort(bfactor(nlabs,d),'ascend');
    end        
    
    [~,I] = sort(dist_dims,'ascend');
    I_unsort = 1:length(I); 
    I_unsort(I) = I_unsort;
    factors = factors(I_unsort);    
end

% Each dimension gets shrunk by the values in factors
new_size = floor(dist_dims./factors);            


Nloc = prod(non_dist_dims);
Ndist = prod(dist_dims);

part_size = prod(new_size)*Nloc*ones(1,nlabs);

is = cell(1,d); IS = cell(1,d);

loc_size = Composite();
start_indices = Composite(); 
end_indices = Composite();

block_to_lab = distributed.zeros(Ndist,1);

% Compute the indices that each worker wants
% If this part is too slow, you might have to work this in to an explicit spmd block
for i=1:nlabs
    l_size = new_size;
    
    [is{:}] = ind2sub(factors,i);
    
    s_indices = zeros(1,d); 
    e_indices = zeros(1,d);  
    idx_range = cell(1,d); 
    for j=1:d
        s_indices(j) = (is{j}-1)*l_size(j)+1; 
        e_indices(j) = (is{j}) * l_size(j);
        % Check if we're at the boundary of the jth dimension, if so adjust block sizes in that dimension
        s = is{j} == factors(j);                
        if s && mod(dist_dims(j),factors(j))>0,                         
            l_size(j) = l_size(j)+ mod(dist_dims(j),factors(j)); 
            e_indices(j) = e_indices(j) + mod(dist_dims(j),factors(j));
        end        
        idx_range{j} = s_indices(j):e_indices(j);
    end
    
    [IS{:}] = ndgrid(idx_range{:});
    
    % Convert d-dimensional indices to 1D indices   
    I_1d = sub2ind(dist_dims, IS{:} );
    
    % Place indices in the correct locations, so each lab knows where to send which chunk of data it currently owns
    block_to_lab(vec(I_1d)) = i;
    
    start_indices{i} = s_indices; end_indices{i} = e_indices;
    
    loc_size{i} = l_size;
    part_size(i) = prod(l_size)*Nloc;
end

spmd,
    codist = codistributor1d(1,part_size,[prod(dims),1]);
    A = redistribute(A,codist);         
    Aloc = getLocalPart(A); 
    % Aloc : nondistributed dimensions x distributed dimensions
    Aloc = reshape(Aloc, Nloc,prod(loc_size));
    loc_size = [non_dist_dims,loc_size];
    block_to_lab = redistribute(block_to_lab,codistributor1d(1,part_size/Nloc,[Ndist,1]));
    block_to_lab_loc = getLocalPart(block_to_lab);

    [imin,imax] = globalIndices(A,1);
    
    % imin, imax range over 1:Ndist
    imin = (imin-1)/Nloc+1; imax = imax/Nloc;    
    iloc = imin:imax;                 
    
    ind_self = find(block_to_lab_loc==labindex);
    Anewloc = zeros(size(Aloc));  
    Inewloc = zeros(size(Aloc,2),d);
    index = 1;
    
    % For each lab, trade off the part of the array that we have with the part of the array it needs, if any             
    for i=1:nlabs
        % If we're at our own lab, keep track of the data we currently have/need
        if i==labindex            
            Anewloc(:,index:index+length(ind_self)-1) = Aloc(:,ind_self);
            idx_nd = cell(1,d);
            [idx_nd{:}] = ind2sub(dist_dims,vec(iloc(ind_self)));               
            Inewloc(index:index+length(ind_self)-1,:) = cell2mat(idx_nd);
            index = index + length(ind_self);
        else
            I_to_send = find(block_to_lab_loc == i);                   
            Arecv = labSendReceive(i,i, Aloc(:, I_to_send));        
            idxrecv = labSendReceive(i,i,iloc(I_to_send));            
            Anewloc(:,index:index+size(Arecv,2)-1) = Arecv;
            
            % Convert received 1D indices to ND indices
            idx_nd = cell(1,d);
            [idx_nd{:}] = ind2sub(dist_dims,vec(idxrecv));                        
            Inewloc(index:index+size(Arecv,2)-1,:) = cell2mat(idx_nd);
            index = index + size(Arecv,2);            
        end
    end    
    % Ensure indices are increasing, from outermost dimension inwards
    [~,I] = sortrows(Inewloc,d:-1:1);
    Anewloc = Anewloc(:,I);
    
    Aloc = vec(Anewloc);
    A = codistributed.build(Aloc,codist,'noCommunication');                    
end            
end
