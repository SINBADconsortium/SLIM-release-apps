function [U] = output_shot_record(U,fwdpara,ntrace_loc,mode)
% how to output shot record, into a long distribute array or save to disk shot by shot

%% output U and rec_idx in composite
% Author : Xiangli, edited by Mathias Louboutin
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
% Date : July, 2015

if mode==1
	if strcmp(fwdpara.type,'marine')
		nr=size(fwdpara.xrec,1);
	else
    	nr = prod(size(fwdpara.xrec));
	end
    nt=length(fwdpara.NyqT);
    nsrc= size(fwdpara.xsrc,2);
    if isfield(fwdpara,'tfire')
        nsrc=1;
    end
elseif mode==2
    nr =size(fwdpara.xsrc,1);
    nt=length(0:fwdpara.dt:fwdpara.T);
else
    nr =prod(fwdpara.n);
    nt=length(fwdpara.nsub);
end

spmd

	if iscell(U)
		U=cell2mat(U);
	end
    codistr	= codistributor1d(1,[ntrace_loc.*nt],[nr*nt*nsrc,1]);
    U		= codistributed.build(U(:),codistr,'noCommunication');
end


end