function [x, info] = inv_RM(m0, Q, b, model, inv_para)

if not(isfield(inv_para,'opts'))
    inv_para.opts = [];
end
if not(isfield(inv_para,'params'))
    inv_para.params = [];
end

% global parameters
iter = inv_para.iter;
nf_full = length(model.freq);
nf = inv_para.nf_batchsize;
ns_grid = length(model.xsrc);
ns = inv_para.ns_batchsize;
Q = repmat(Q, [1,1,nf_full]);
nrec = length(model.xrec);
spmd
    Q = codistributed(Q, codistributor1d(3));
end
Q = vec(Q);

innerloop = 0;
outerloop = 0;
success = 0;
infi_loop = 0;
innerloop_list=zeros(iter,1);
residual_list=zeros(2*iter,1);
tau_list=zeros(2*iter,1);
lambda_list=zeros(2*iter,1);

% make subsampling operator
fidx_full = randperm(nf_full);
fidx = sort(fidx_full(1:nf),'ascend');
opFreqSub = opRestriction(nf_full,fidx);

source_encoding_mat = make_encoding_mat(ns_grid,ns,inv_para.encoding_type);
opSourceSub = opRightMulMat([1,ns_grid], source_encoding_mat);

% inversion
while (innerloop < iter) && (not(success)) && (not(infi_loop))
    % re-make subsampling operator if using renrewal
    if inv_para.f_renewal
        fidx_full = randperm(nf_full);
        fidx = sort(fidx_full(1:nf),'ascend');
        opFreqSub = opRestriction(nf_full,fidx);
    end
    if inv_para.s_renewal
        source_encoding_mat = make_encoding_mat(ns_grid,ns,inv_para.encoding_type);
        opSourceSub = opRightMulMat([1,ns_grid], source_encoding_mat);
    end
    
    opSub_Q = oppKron2Lo(opFreqSub, opKron(opSourceSub,opDirac(ns_grid)));
    opSub_b = oppKron2Lo(opFreqSub, opKron(opSourceSub,opDirac(nrec)));

    % subsample Born operator
    Q_sub = opSub_Q*Q;
    Q_sub = invvec(Q_sub, [ns_grid, ns, nf]);
    model_sub = model;
    model_sub.freq = transpose(opFreqSub*vec(model.freq));
    opBorn_sub = oppDF_old(m0(:), Q_sub, model_sub);
    A_sub = inv_para.scale_op*opBorn_sub*inv_para.opGetReal*inv_para.C;

    % subsample data
    b_sub = opSub_b*b;
    
    % sparse inversion
    [x, ~, ~, info] = spgl1(A_sub, b_sub, inv_para.tau, inv_para.sigma, inv_para.x0, ...
                            inv_para.opts, inv_para.params);
    
    outerloop = outerloop + 1;
    innerloop = innerloop + info.iter;
    innerloop_list(outerloop) = info.iter;
    residual_list(innerloop-info.iter+1:innerloop) = info.rNorm2;
    tau_list(innerloop-info.iter+1:innerloop) = info.xNorm1;
    lambda_list(innerloop-info.iter+1:innerloop) = info.lambda;
    inv_para.x0 = x;
    inv_para.tau = norm(x, 1);
    success = ((info.stat >= 1)&&(info.stat <= 4))||(info.stat == 7);
    infi_loop = isequal(info.iter,0);
end

% summerize renewal info
info.innerloop_list=innerloop_list(1:outerloop);
info.residual_list=residual_list(1:innerloop);
info.tau_list=tau_list(1:innerloop);
info.lambda_list=lambda_list(1:innerloop);

if innerloop >= iter
    info.quitStat = 'full_iter';
elseif success == 1
    info.quitStat = 'BP or subOpBP solution found';
elseif infi_loop ==1
    info.quitStat = 'Trapped in infinite loop.';
end