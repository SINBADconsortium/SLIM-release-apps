function x = dsweep_old(R,idx,w,x,b,ns,n_threads)
   for k=1:ns
      % forward sweep
      x = sweepR_mex(R,idx,x,b,w,1,n_threads);
      % backward sweep
      x = sweepR_mex(R,idx,x,b,w,-1,n_threads);
   end
end