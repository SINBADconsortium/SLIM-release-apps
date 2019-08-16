n = [100 100];
o = [0 0];
d = [10 10];
nb = [20 20];
Nt = prod(n);

k = randn(n) / 100;
k = k(:);
x = randn(prod(n),1);


[H, A, Bk, C, Dc] = Helm2D(k, o, d, n, nb);

dk = randn(prod(n),1)/1000;

dk = reshape(dk,n);
dk(1:20,:) = 0;
dk(:,1:20) = 0;
dk(end-19:end,:) = 0;
dk(:,end-19:end) = 0;
dk = dk(:);

h  = 10.^([-6:1:0]);


for i = 1: length(h)
    ki = k + h(i) * dk;
    Hi = Helm2D(ki, o, d, n, nb);
    r1(:,i) = Hi*x - H*x;
    r2(:,i) = A * (Bk * h(i) * dk) + (Dc.*spdiags([x,x,x],[-1,0,1],Nt,Nt)) * h(i) * dk;
    nr1(i)  = norm(r1(:,i));
    nr2(i)  = norm(r1(:,i) - r2(:,i));

end

figure;loglog(h,h);hold on; loglog(h,h.^2); loglog(h,nr1); loglog(h,nr2);
legend('h','h^2','e0','e1')

