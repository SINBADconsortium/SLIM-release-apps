function [R,idx] = pRtransp(R,idx)
% parallel version of transpose band-storage matrix
%
% use:
%   [R,idx] = pRtransp(R,idx)
%
% Author: Zhilong Fang
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


n = size(R,1);
%keyboard

spmd
    codist = getCodistributor(R);
    Part   = codist.Partition;
    dim    = codist.Dimension;
    RLocal = getLocalPart(R);
    id     = labindex;
    nlab   = numlabs;
    for k = 1:length(idx);
        RLocal(:,k) = circshift(RLocal(:,k),idx(k));
        if idx(k) > 0
            if id < nlab
                labSend(RLocal(1:idx(k),k),id+1);
            end
            if id > 1
                Rrec = labReceive(id-1);
            end

            if id == nlab
                labSend(RLocal(1:idx(k),k),1);
            end
            
            if id == 1
                Rrec = labReceive(nlab);
            end
            RLocal(1:idx(k),k) = Rrec;
        else if idx(k) < 0
                if id > 1
                    labSend(RLocal(size(RLocal,1)+idx(k)+1:size(RLocal,1),k),id-1);
                end
                if id < nlab
                    Rrec = labReceive(id+1);
                end
 
                if id == 1
                    labSend(RLocal(size(RLocal,1)+idx(k)+1:size(RLocal,1),k),nlab);
                end
                if id == nlab
                    Rrec = labReceive(1);
                end
                RLocal(end+idx(k)+1:end,k) = Rrec;
            end
        end
        
        idx(k) = -idx(k);
    end
    R = codistributed.build(RLocal,codist,'noCommunication');
    
end



