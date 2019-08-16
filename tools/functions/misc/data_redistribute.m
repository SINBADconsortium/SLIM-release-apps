function D = data_redistribute(D,nrec,nsrc,nfreq,mode)
    
    switch mode
      case 'srcfreq'
        D = reshape(D,nrec,nsrc*nfreq);
      case 'freq'
        D = reshape(D,nrec*nsrc,nfreq);
      otherwise 
        error('Unrecognized mode');
    end    
end
