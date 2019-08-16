function [x, info] = Mig_SRM_IG(m0, S, P, model, inv_para)

if not(isstruct(P))
    Ps = P;
    Pr = Ps;
else
    Ps = P.Ps;
    Pr = P.Pr;
end
clear P

% global parameters
iter = inv_para.iter;
src_polarity = inv_para.src_polarity;
nf_full = length(model.freq);
nf = inv_para.nf_batchsize;
ns_grid = length(model.xsrc);
ns = inv_para.ns_batchsize;
nrec = length(model.xrec);
if isfield(inv_para,'snapshot_file')
    opSave = opSaveSnapshot(model.n(1)*model.n(2),inv_para.snapshot_file);
end

if isfield(inv_para,'usefracdiff')
	usefracdiff = inv_para.usefracdiff;
else
	usefracdiff = 0;
end

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

if isequal(inv_para.encoding_type,'dirac')
    srmask = inv_para.srmask*source_encoding_mat;
	if usefracdiff
		freq_sub = opFreqSub*vec(model.freq);
		opSRMask = oppKron2Lo(opDiag((2*pi*freq_sub(:)).^(-0.5)), opDiag(srmask(:)));
	else
		opSRMask = oppKron2Lo(opDirac(nf), opDiag(srmask(:)));
	end
else
	if usefracdiff
		freq_sub = opFreqSub*vec(model.freq);
		opSRMask = oppKron2Lo(opDiag((2*pi*freq_sub(:)).^(-0.5)), opDirac(nrec*ns));
	else
		opSRMask = 1;
	end
end

% inversion
while (innerloop < iter) && (not(success)) && (not(infi_loop))
    % re-make subsampling operator if using renrewal
    if inv_para.f_renewal && (outerloop > 0)
        fidx_full = randperm(nf_full);
        fidx = sort(fidx_full(1:nf),'ascend');
        opFreqSub = opRestriction(nf_full,fidx);
    end
    if inv_para.s_renewal && (outerloop > 0)
        source_encoding_mat = make_encoding_mat(ns_grid,ns,inv_para.encoding_type);
        opSourceSub = opRightMulMat([1,ns_grid], source_encoding_mat);
		if isequal(inv_para.encoding_type,'dirac')
		    srmask = inv_para.srmask*source_encoding_mat;
			if usefracdiff
				freq_sub = opFreqSub*vec(model.freq);
				opSRMask = oppKron2Lo(opDiag((2*pi*freq_sub(:)).^(-0.5)), opDiag(srmask(:)));
			else
				opSRMask = oppKron2Lo(opDirac(nf), opDiag(srmask(:)));
			end
		else
			if usefracdiff
				freq_sub = opFreqSub*vec(model.freq);
				opSRMask = oppKron2Lo(opDiag((2*pi*freq_sub(:)).^(-0.5)), opDirac(nrec*ns));
			else
				opSRMask = 1;
			end
		end
    end
    inv_para.opLeft = inv_para.opLeft*opSRMask;
    
    opSub_S = oppKron2Lo(opFreqSub, opKron(opSourceSub,opDirac(ns_grid)));
    opSub_P = oppKron2Lo(opFreqSub, opKron(opSourceSub,opDirac(nrec)));

    % subsample Born operator
    S_sub = opSub_S*S;
    S_sub = sourcegrid(ns_grid, ns, src_polarity, S_sub);
    P_sub = opSub_P*Ps;
    P_sub = sourcegrid(nrec, ns, src_polarity, P_sub);
    
    model_subS = model;
    model_subS.freq = transpose(opFreqSub*vec(model.freq));
    model_subP = model_subS;
    model_subP.xsrc = model_subP.xrec;

    opBornS = oppDF_old(m0(:), S_sub, model_subS);
    opBornP = oppDF_old(m0(:), P_sub, model_subP);

    % subsample data
    b_sub = opSRMask*opSub_P*Pr;
    
    % params
    if isfield(inv_para,'hub')
        params.hub = inv_para.hub;
    end
    
    % sparse inversion
    if not(inv_para.estS)
        % S should contain the true source signature
        if isequal(model_subS.xsrc(:), model_subP.xsrc(:))
            Q_sub = S_sub-P_sub;
            model_subQ = model_subS;
            opBornQ = oppDF_old(m0(:), Q_sub, model_subQ);
            A_sub = inv_para.opLeft*opBornQ*inv_para.opRight;
        else
            A_sub = inv_para.opLeft*(opBornS-opBornP)* ...
                    inv_para.opRight;
        end
        [x,~,~,info] = spgl1(A_sub, b_sub, inv_para.tau, inv_para.sigma, ...
                             inv_para.x0, inv_para.opts, params);
    else
        % S is the impulsive wavefield I
        params.opLeft = inv_para.opLeft;
        params.opRight = inv_para.opRight;
        params.opBornI = opBornS;
        params.opBornP = opBornP;
        params.I_sub = S_sub;
        params.P_sub = P_sub;
        params.modelI = model_subS;
        params.modelP = model_subP;
        params.m0 = m0;
        params.b = b_sub;
        params.nf = nf;
        params.fidx = fidx;
        params.nrec = nrec;
        params.nsrc = ns;
        params.source_encoding_mat = source_encoding_mat;
        [x,~,~,info] = spgl1(@projMigSRM_IG, b_sub, inv_para.tau, ...
                             inv_para.sigma, inv_para.x0, inv_para.opts, ...
                             params);
    end
    
    outerloop = outerloop + 1;
    innerloop = innerloop + info.iter;
    innerloop_list(outerloop) = info.iter;
    residual_list(innerloop-info.iter+1:innerloop) = info.rNorm2;
    tau_list(innerloop-info.iter+1:innerloop) = info.xNorm1;
    lambda_list(innerloop-info.iter+1:innerloop) = info.lambda;
    inv_para.x0 = x;
    if not(isequal(inv_para.tau,inf))
        inv_para.tau = norm(x,1);
    end
    success = ((info.stat >= 1)&&(info.stat <= 4))||(info.stat == 7);
    infi_loop = isequal(info.iter,0);
    if isfield(inv_para,'snapshot_file')
        opSave*inv_para.opRight*x;
    end
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