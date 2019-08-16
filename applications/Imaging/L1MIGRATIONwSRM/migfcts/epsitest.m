function epsitest(epsiresult)
% Syntax
% epsitest(epsiresult)
%
% Description
% Test EPSI result. Will display the relative residule and a warn if it is larger
% than 0.2. This program will also write variable "multiple" to EPSI result.
%
% Input list:
% epsiresult: the full path of EPSI solmat.
%
% Output list: None.
% 
% Author: Ning Tu
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: Feb/14/2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

%load data
load(epsiresult,'q_est');
q=q_est(:,end);
load(epsiresult,'initial_data');
load(epsiresult,'primaryIR');
load(epsiresult,'options');
load(epsiresult,'dt')

if isempty(options.topmuteT)
    topmute=0;
else
    topmute=floor(options.topmuteT/dt);
end

use_oblique=options.useOblique;
%% test time domain EPSI operators

disp('Time domain EPSI operator test...')

[nt,nrec,nshot] = size(initial_data);

if options.parallel
	initial_data = distSort_shot2timeSlice(distributed(initial_data));
	primaryIR = distSort_shot2timeSlice(distributed(primaryIR));

	[OP_EPSI OP_P_TERM] = opEPSI_dist(initial_data,q,topmute,dt,use_oblique);

	disp('Testing OP_EPSI...')
	p_recover = OP_EPSI*primaryIR(:);
	ratio_p = norm(p_recover(:)-initial_data(:))/norm(initial_data(:));
	fprintf('Relative data misfit : %f%%\n',100*ratio_p)

	if ratio_p > 0.2
    	disp('Relative residual above 0.2. You may want to double check your EPSI result.')
    end

    multiple = -OP_P_TERM*primaryIR(:);
    multiple = reshape(undist(multiple), nrec, nshot, nt);
    multiple = shiftdim(multiple,2);
    save(epsiresult,'multiple','-append')
else
	[OP_EPSI OP_P_TERM] = opEPSI_serial(initial_data,q,topmute,dt,use_oblique);

	disp('Testing OP_EPSI...')
	p_recover = OP_EPSI*primaryIR(:);
	ratio_p = norm(p_recover(:)-initial_data(:))/norm(initial_data(:));
	fprintf('Relative data misfit : %f%%\n',100*ratio_p)

	if ratio_p > 0.2
    	disp('Relative residual above 0.2. You may want to double check your EPSI result.')
    end

    multiple = -OP_P_TERM*primaryIR(:);
    multiple = reshape(multiple, nt, nrec, nshot);
    save(epsiresult,'multiple','-append')
end