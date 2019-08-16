function [status,result] = fwdwrapper(fname_buoy,fname_bulk,model,params,options)
%
% Wrapper for forward modeling, input= model; output=data (.su format)
%             Usage:fwdwrapper(buoy,bulk,model,params)
% buoy= buoyancy; bulk=bulkmod,model=model;
% params=changing the parameters by the user as you would do in a parameter file(.par).


% setting the path
% result directory of the forward modelling test

if isfield(model,'randseq')
    isrand = model.randseq;
else
    isrand = 0;
end

if isfield(model,'simultsrc')
    ismulti = model.simultsrc;
else
    ismulti = 0;
end

if isfield(options,'label')
    label = [options.label, '_fwd'];
else
    label = 'fwd_results_temp';
end



if ismulti == 0 & isrand == 0
    if ~exist(label,'dir')
        mkdir(label);
    end
    curdir = pwd;
    cd(label);
    flagReg = 1;
else
    flagReg = 0;
end

%  matlab scripts to generate input files for iWAVE


% write script to generate the header file
if ismulti == 0
	iwave_wrdatahdr(model,label);
else
	iwave_wrdatahdr_sim(model,label);
end
% Generate source file if it is needed
if isfield(model,'t0')
	if ismulti == 0
    		iwave_wrsrc(model,label);
	else
		iwave_wrsrc_sim(model,label);
	end
end

fwd_wrparam(model,label,params,fname_buoy,fname_bulk);


nlab    = options.nlab;

% make the header file for output file
mkh = ['. ./' label '_mkhdr.sh'];
[status,result]=system(mkh);  % execute the shell script to generate the header file


% Call the iwave function

%mpirun_Local = options.MPIRUN_Command;
if isfield(options,'MPICommand')
        mpirun_Local = options.MPICommand;
else
        mpirun_Local = iwave_MPIcommand(options.nlab,flagReg);
end

if nlab == 1
        pr = ['asgfwd.x par=' label '_par.par'];
else
        pr  = [mpirun_Local ' asgfwd.x par=' label '_par.par'];
end

% if nlab == 1 || isrand ~= 0 || ismulti ~= 0
%     pr= ['asgfwd.x par=' label '_par.par'];    % for single processor
%     [status,result]=system(pr)
% else
%     if options.torque < 1
%        pr = ['mpirun -np ' num2str(nlab) ' asgfwd.x par=' label '_par.par'];   %for mpi
%     else
%     %getenv('PBS_NODEFILE')
%         pr = ['mpiexec.hydra -np ' num2str(nlab) ' -hostfile ' getenv('PBS_NODEFILE') ' asgfwd.x par=' label '_par.par'];   %for mpi
%     end
%     [status,result]=system(pr)
% end
[status,result]=system(pr);

if status~=0
        error(result);
end


end
