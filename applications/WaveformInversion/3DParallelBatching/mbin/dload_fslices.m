function D = dload_fslices(fslice_name_func,dsize)
% DLOAD_FSLICES - Loads (in a distributed fashion) frequency slices
% + assembles the entire data volume in parallel
%
% Curt Da Silva, 2015
%
%   D = dload_fslices(fslice_name_func,dsize);
%
% Input:
%   fslice_name_func    - function handle that returns the full path to the frequency slice for each specified frequency index
%   dsize               - 3 length vector of nrec x nsrc x nfreqs
% 
% Output:
%   D                   - codistributed data volume
%
    nrec = dsize(1); nsrc = dsize(2); nfreq = dsize(3);    
    nlabs = parpool_size();        
    assert(nlabs > 0, 'need open matlabpool');            
    
    spmd
        if labindex <= nfreq
            filepath = fslice_name_func(labindex);
            [~,~,ext] = fileparts(filepath);
            switch ext
              case '.mat'
                Dloc = load(filepath);
                fields = fieldnames(Dloc);
                if length(fields) > 1
                    error('Too many variables in file');
                end
                Dloc = Dloc.(fields{1});            
              case '.rsf'
                Dloc = rsf_read_all(filepath);
                Dloc = reshape(Dloc,nrec,nsrc);
              otherwise 
                error('Unsupported file extension');
            end
             assert(size(Dloc,1)==nrec && size(Dloc,2)==nsrc,['Unexpected freq slice size, expected ' num2str(nrec) ' x ' num2str(nsrc) ' but got ' num2str(size(Dloc,1)) ' x ' num2str(size(Dloc,2))]);
        else
            Dloc = zeros(nrec,0);
        end
        part = zeros(1,nlabs);
        part(1:nfreq) = nsrc;
        codistr = codistributor1d(2,part,[nrec,nsrc*nfreq]);
        D = codistributed.build(Dloc,codistr,'noCommunication');
        D = redistribute(D,codistributor1d(2,codistributor1d.unsetPartition,[nrec,nsrc*nfreq]));
    end
end
