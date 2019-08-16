function idx = randjet(n,k)

idx(k) = 0;

if n/k < 2
   idx = randperm(n,k);
else
   j = floor(n/k);
   ke = j * (k-1);
   kk = 0;
   for i = 1:k-1
       ii = randperm(j,1);
       idx(i) = kk + ii;
       kk     = kk + j;
   end
   ii = randperm(n-ke,1);
   idx(k) = ii + ke;
end

idx = sort(idx, 'ascend');
