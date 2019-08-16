model.n    = [101 101];
model.d    = [10 10];
model.o    = [0 0];
model.nb   = [20 20];
model.f0   = 10;
model.t0   = .5;
model.freq = 0:0.2:29.8;
model.xsrc = 500;
model.zsrc = 100;
model.xrec = 0:10:1000;
model.zrec = 100;

params.computeLU = 0;
params.nthreads  = 4;

Q = speye(length(model.xsrc));

v = 2000 * ones(model.n);
m = 1e6./v(:).^2;

D = F(m,Q,model,params);
D = gather(D);
D = reshape(D,length(model.xrec),length(model.freq));
D = D.';
for i = 1:size(D,2)
    TD(:,i) = real(ifft((D(:,i))));
end
figure;imagesc(TD);colormap gray

