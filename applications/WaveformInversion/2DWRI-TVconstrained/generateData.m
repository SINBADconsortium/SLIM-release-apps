function d = generateData(model,q,Ps,Pr)
%returns simulated data and operator for definine M(u) in the Helmoltz eq.

%notation
nf = model.nf;
ns = model.ns;
nr = model.nr;
freq = model.freq;
mtrue = model.mtrue;
dsig = model.dsig;

%generate data matrices (one for each frequency)
Q = speye(ns); %source weight
d = cell(nf,1);
for v = 1:nf
    Uv = (model.L + model.M(freq(v),mtrue(:)))\(q(v)*(Ps'*Q));
    d{v} = Pr*Uv;
    %generate noise proportional to norm of clean data
    sig = dsig*norm(d{v},'fro')/sqrt(ns*nr);
    rng(model.rngnoise);
    data_noise = sig*randn(size(d{v}));
    d{v} = d{v} + data_noise;
end