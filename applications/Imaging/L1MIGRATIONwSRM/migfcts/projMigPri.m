function f = projMigPri(x, g, params)

% declare global and persistent variable
global S_f S_f_full fidx
persistent pri_mod x_prev fidx_prev sem_prev

update_x_flag = not(isequal(x_prev,x));
x_prev = x;

% load para
opLeft = params.opLeft;
opRight = params.opRight;
opBornI = params.opBornI;
nf = params.nf;
nrec = params.nrec;
nsrc = params.nsrc;
b = params.b;
fidx = params.fidx;
update_f_flag = not(isequal(fidx_prev,fidx));
fidx_prev = fidx;
sem = params.source_encoding_mat;
update_s_flag = not(isequal(sem_prev,sem));
sem_prev = sem;

update_flag = update_x_flag || (update_f_flag || update_s_flag);

% form operator
A_I = opLeft*opBornI*opRight;
if isempty(pri_mod)
    pri_mod = A_I*x;
else
    if update_flag
        pri_mod = A_I*x;
    end
end

% compute frequency spectrum of q
dm = opRight*x;
if not(all(not(dm)))
    pri_mod = invvec(pri_mod,[nrec*nsrc,nf]);
    b = invvec(b,[nrec*nsrc,nf]);
    for i=1:nf
        S_f(i) = gather(pri_mod(:,i)'*b(:,i)/(pri_mod(:,i)'*pri_mod(:,i)));
    end
    pri_mod = pri_mod(:);
    b = b(:);
end

opConvS = oppKron2Lo(opDiag(S_f),opDirac(nrec*nsrc));
A = opConvS*A_I;
if isempty(g)
    f = opConvS*pri_mod;
else
    f = A'*g;
end

% update global S_f_full
S_f_full(fidx) = S_f;