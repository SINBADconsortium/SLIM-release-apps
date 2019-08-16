d = 4; n = 40; ndist_dims = 2;
dims = randi(n,1,d) + 10;

X = pSPOT.utils.distVectorize(distributed.randn(prod(dims),1));

Xg = reshape(gather(X),dims);
[Xdist,loc_size,start_indices,end_indices] = distribute_nd_array(X,dims,ndist_dims);

spmd,
    i1 = start_indices(1):end_indices(1);
    i2 = start_indices(2):end_indices(2);
    Xloc = getLocalPart(Xdist);
    e = norm(vec(Xloc) - vec(Xg(:,:,i1,i2)));    
    e = pSPOT.utils.global_sum(e);
end
assert( e{1} == 0);

%%
d = 4; n = 40;
dims = randi(n,1,d) + 10;

X = pSPOT.utils.distVectorize(distributed.randn(dims));

Xg = reshape(gather(X),dims);
[X,loc_size,start_indices,end_indices] = distribute_nd_array(X,dims,ndist_dims);

spmd,
    i1 = start_indices(1):end_indices(1);
    i2 = start_indices(2):end_indices(2);
    i3 = start_indices(3):end_indices(3);
    Xloc = getLocalPart(X);
    e = norm(vec(Xloc) - vec(Xg(:,i1,i2,i3)));   
    e = pSPOT.utils.global_sum(e);    
end
assert( e{1} == 0);

%%
d = 4; n = 40;
dims = randi(n,1,d) + 10;

X = pSPOT.utils.distVectorize(distributed.randn(dims));

Xg = reshape(gather(X),dims);
[X,loc_size,start_indices,end_indices] = distribute_nd_array(X,dims,ndist_dims);

spmd,
    i1 = start_indices(1):end_indices(1);
    i2 = start_indices(2):end_indices(2);
    i3 = start_indices(3):end_indices(3);
    i4 = start_indices(4):end_indices(4);
    Xloc = getLocalPart(X);
    e = norm(vec(Xloc) - vec(Xg(i1,i2,i3,i4))); 
    e = pSPOT.utils.global_sum(e);   
end
assert( e{1} == 0);


%% 
d = 5; 
dims = randi(20,1,d) + 10;

X = pSPOT.utils.distVectorize(distributed.randn(dims));

Xg = reshape(gather(X),dims);
[X,loc_size,start_indices,end_indices] = distribute_nd_array(X,dims,ndist_dims);

spmd,
    i1 = start_indices(1):end_indices(1);
    i2 = start_indices(2):end_indices(2);
    i3 = start_indices(3):end_indices(3);
    Xloc = getLocalPart(X);
    e = norm(vec(Xloc) - vec(Xg(:,:,i1,i2,i3)));   
    e = pSPOT.utils.global_sum(e);    
end
assert( e{1} == 0);


%% 
d = 5; n = 20;
dims = randi(n,1,5) + 10;

X = pSPOT.utils.distVectorize(distributed.randn(dims));

Xg = reshape(gather(X),dims);
[X,loc_size,start_indices,end_indices] = distribute_nd_array(X,dims,ndist_dims);

spmd,
    i1 = start_indices(1):end_indices(1);
    i2 = start_indices(2):end_indices(2);
    i3 = start_indices(3):end_indices(3);
    i4 = start_indices(4):end_indices(4);
    Xloc = getLocalPart(X);
    e = norm(vec(Xloc) - vec(Xg(:,i1,i2,i3,i4)));   
    e = pSPOT.utils.global_sum(e);   
end
assert( e{1} == 0);

%%
