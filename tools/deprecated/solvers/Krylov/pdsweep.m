function x = pdsweep(R,idx,x,b,w,n,n_threads)

   N  = prod(n(1:end-1));

   % setup
   dist = getCodistributor(x);
   Rloc = getLocalPart(R);
   xloc = getLocalPart(x);
   bloc = getLocalPart(b);
   idxloc = idx + (labindex~=1)*N;

   % forward sweep
   xloc = sweepR_mex(Rloc,idxloc,xloc,bloc,w,1,n_threads);  
      assert(not(any(isnan(xloc))),'dsweep: an element of xloc, the local part of x, is NaN after the forward sweep.');

   % get halos
   if labindex<numlabs; labTo = labindex + 1; else labTo = []; end;
   if labindex>1; labFrom = labindex - 1; else labFrom = []; end;
   halol = labSendReceive(labTo, labFrom, xloc(end-2*N+1:end));

   if labindex>1; labTo = labindex - 1; else labTo = []; end;
   if labindex<numlabs; labFrom = labindex + 1; else labFrom = []; end;
   halor = labSendReceive(labTo, labFrom, xloc(1:2*N));

   % update halos
   if labindex > 1; xloc(1:2*N)=(xloc(1:2*N)+halol)/2; end;
   if labindex < numlabs; xloc(end-2*N+1:end)=(xloc(end-2*N+1:end)+halor)/2; end;

   % backward sweep
   xloc = sweepR_mex(Rloc,idxloc,xloc,bloc,w,-1,n_threads);
      assert(not(any(isnan(xloc))),'dsweep: an element of xloc, the local part of x, is NaN after the backward sweep.');

   % get halos
   if labindex<numlabs; labTo = labindex + 1; else labTo = []; end;
   if labindex>1; labFrom = labindex - 1; else labFrom = []; end;
   halol = labSendReceive(labTo, labFrom, xloc(end-2*N+1:end));

   if labindex>1; labTo = labindex - 1; else labTo = []; end;
   if labindex<numlabs; labFrom = labindex + 1; else labFrom = []; end;
   halor = labSendReceive(labTo, labFrom, xloc(1:2*N));

   % update halos
   if labindex > 1; xloc(1:2*N)=(xloc(1:2*N)+halol)/2; end;
   if labindex < numlabs; xloc(end-2*N+1:end)=(xloc(end-2*N+1:end)+halor)/2; end;

   % wrap up
   x = codistributed.build(xloc,dist,'noCommunication');

end
