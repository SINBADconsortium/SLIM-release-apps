function [output] = LowRank_2D(A,para,options,r,sigmafact,initfact)
% This is the main function for doing interpolation and / or denoising.
%
% use:
%   [output] = LowRank_2D(A,para,opts,rank,sigmafact,initfact)
%
% input:
% A         - observed data.
% para      - This will contain the information about data.
% opts      - parameters for SPGL1
% r         - rank of the system to be used in the processing. It can be a
%             single number or a vector equal to the number of frequency.
% sigmafact - defines the percentage of noise tolerance in the output
%             results
% initfact  - This factor control the amplitude of intial random guess. it
%             is used to make the amplitude of intial guess to be equal
%             to the amplitude of input data.
% output:
%  output   - Regularized/ interpolated/ denoised data

% Author: Rajiv Kumar
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: April, 2014
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

% Reading parameters for the experiment
K    = size(A);
nR   = K(1);
nC   = K(2);
nf   = K(3);
D    = reshape(A,nR,nC,nf);

if para.parallel==1 % run code in parallel mode
    freq = distributed(1:nf);
    D    = distributed(D);
    
    freqidx  = 1:1:nf;
    freqidx = distributed(freqidx);
    if para.rank==1
        Rank = distributed(r);
    end
    if para.penalty==1
        logname = [options.resultdir 'log_studentst_reg_freqindx'];
    else
        logname = [options.resultdir 'log_leastsquares_reg_freqindx'];
    end
    spmd
        if para.rank==1
            rankloc = getLocalPart(Rank);
        end
        codistr   = codistributor1d(2,[],[nR*nC,length(freq)]);
        freqloc   = getLocalPart(freq);
        freqidxloc = getLocalPart(freqidx);
        nfreqloc  = length(freqloc);
        outputloc = zeros(nR*nC,nfreqloc);
        Dloc      = getLocalPart(D);
        for k = 1:nfreqloc
            if para.rank==1
                rank = rankloc(k);
            else
                rank = r;
            end
            fid = fopen([logname num2str(freqidxloc(k)) '.dat'],'w');
            if para.algo==1
                Dlocal         = Interp_Denoise(Dloc(:,:,k),para,options,rank,sigmafact,initfact,fid);
            else
                Dlocal         = Reg_Interp(Dloc(:,:,k),para,options,rank,sigmafact,initfact,fid);
            end
            outputloc(:,k) = Dlocal(:);
        end
        output = codistributed.build(outputloc,codistr,'noCommunication');
    end
    output = gather(output);
    
else % run code in serial mode
    
    freqidx  = 1:1:nf;

    if para.penalty==1
        logname = [options.resultdir 'log_studentst_reg_freqindx'];
    else
        logname = [options.resultdir 'log_leastsquares_reg_freqindx'];
    end
    
    output = zeros(nR*nC,nf);
    for k = 1:nf
        if para.rank==1
            rank = r(k);
        else
            rank = r;
        end
        fid = fopen([logname num2str(freqidx(k)) '.dat'],'w');
        if para.algo==1
            Dlocal         = Interp_Denoise(D(:,:,k),para,options,rank,sigmafact,initfact,fid);
        else
            Dlocal         = Reg_Interp(D(:,:,k),para,options,rank,sigmafact,initfact,fid);
        end
        output(:,k) = Dlocal(:);
    end
       
    output = gather(output);
    
end
end
