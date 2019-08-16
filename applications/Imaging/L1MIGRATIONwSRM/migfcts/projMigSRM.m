function f = projMigSRM(x, g, params)

% declare global and persistent variable
global S_f S_f_full fidx
persistent pri_mod mul_mod x_prev fidx_prev sem_prev

% check whether x is updated, if not, it may not be necessary to
% remodel pri_mod and mul_mod
update_x_flag = not(isequal(x_prev,x));
x_prev = x;

% load para
opLeft = params.opLeft;
opRight = params.opRight;
opBornI = params.opBornI;
opBornP = params.opBornP;
I_sub = params.I_sub;
P_sub = params.P_sub;
modelI = params.modelI;
modelP = params.modelP;
m0 = params.m0;
b = params.b;
nf = params.nf;
fidx = params.fidx;
update_f_flag = not(isequal(fidx_prev,fidx));
fidx_prev = fidx;
nsrc = params.nsrc;
nrec = params.nrec;
sem = params.source_encoding_mat;
update_s_flag = not(isequal(sem_prev,sem));
sem_prev = sem;

update_flag = update_x_flag || (update_f_flag || update_s_flag);

% form operator
A_I = opLeft*opBornI*opRight;
A_P = opLeft*opBornP*opRight;
if isempty(pri_mod) || isempty(mul_mod)
    pri_mod = A_I*x;
    mul_mod = -A_P*x;
else
    if update_flag
        pri_mod = A_I*x;
        mul_mod = -A_P*x;
    end
end

% compute frequency spectrum of q
dm = opRight*x;
if not(all(not(dm)))
    pri_mod = invvec(pri_mod,[nrec*nsrc,nf]);
    res = invvec(b-mul_mod,[nrec*nsrc,nf]);
    for i=1:nf
        S_f(i) = gather(pri_mod(:,i)'*res(:,i)/(pri_mod(:,i)'*pri_mod(:,i)));
    end
    S_f_full(fidx) = S_f;
    pri_mod = pri_mod(:);
end

opConvS_I = oppKron2Lo(opDiag(S_f),opDirac(size(I_sub,1)*size(I_sub,2)));
S_sub = opConvS_I*I_sub(:);
S_sub = invvec(S_sub,size(I_sub));
if isequal(modelI.xsrc(:), modelP.xsrc(:))
    Q_sub = S_sub-P_sub;
    modelQ = modelI;
    opBornQ = oppDF_old(m0(:),Q_sub,modelQ);
    A = opLeft*opBornQ*opRight;
else
    opBornS = oppDF_old(m0(:),S_sub,modelI);
    A = opLeft*(opBornS-opBornP)*opRight;
end
opConvS = oppKron2Lo(opDiag(S_f),opDirac(nrec*nsrc));
if isempty(g)
    f = opConvS*pri_mod+mul_mod;
else
    f = A'*g;
end