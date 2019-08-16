function d = TraceNorm_dual_macq(x,weights,params)

% dual of trace norm is operator norm i.e maximum singular value
if params.mode == 1
   E = reshape(x,params.numr,params.numc);
   d = svds(gather(E),1);

elseif params.mode == 2 % svd on full matrix
   fps1 = reshape(x(1:params.mhnumr*params.mhnumc),params.mhnumr,params.mhnumc);
   fps2 = reshape(x(params.mhnumr*params.mhnumc+1:end),params.mhnumr,params.mhnumc);
   E = [fps1;fps2];
   d = svds(gather(E),1);
    
elseif params.mode == 3 % random SVD
   E = reshape(x,params.numr,params.numc);
   row_oversample = ceil(2*params.nr*log(params.nr));
   row_order = randperm(params.numr);
   ind = row_order(1:row_oversample);
   Ered = E(ind,:);
   d = svds(gather(Ered),1);

else % random SVD with irbl svd library
   E = reshape(x,params.numr,params.numc);
   row_oversample = ceil(2*params.nr*log(params.nr));
   row_order = randperm(params.numr);
   ind = row_order(1:row_oversample);
   Ered = E(ind,:);
   opts.K = params.svdnum;
   opts.tol = params.svdtol;
   d = irblsvds(Ered,opts);
end

end % function end

