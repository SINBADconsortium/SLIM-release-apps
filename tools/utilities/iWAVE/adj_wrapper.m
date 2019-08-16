function [status,result] = adj_wrapper(fname_buoy,fname_bulk,model,datafile,params,options)

% setting the path
if isfield(options,'label')
    label = [options.label,'_adj'];
else
    label = 'adj_results_temp';
end

if ~exist(label,'dir')
    mkdir(label);
end

% write source file if possible
if isfield(model,'romb')
   	isromb = model.romb;
else
    	isromb = 0;
end

if isfield(model,'simultsrc')
	ismulti = model.simultsrc;
else
	ismulti = 0;
end

if isfield(model,'encod')
	isphase = model.encod;
else
	isphase = 0;
end

if isfield(model,'randseq')
      isrand = model.randseq;
else
      isrand = 0;
end

if ismulti == 0 && isrand == 0 && isphase == 0 && isromb == 0
   curdir = pwd;
   cd(label);
   flagReg = 1;
else
   flagReg = 0;
end



if isfield(model,'t0')
	if ismulti == 0 && isphase == 0 && isromb ==0
    		iwave_wrsrc(model,label);
	else
		iwave_wrsrc_sim(model,label);
	end
end

% write parameter file
adj_wrparam(model,label,params,fname_buoy,fname_bulk,datafile);

% detect number of labs
nlab    = options.nlab;



% make header file
mkh = ['. ./' label '_mkhdr.sh'];
[status,result]=system(mkh);  % execute the shell script to generate the header file

% call iwave function
%mpirun_Local = options.MPIRUN_Command;
if isfield(options,'MPICommand')
        mpirun_Local = options.MPICommand;
else
        mpirun_Local = iwave_MPIcommand(options.nlab,flagReg);
end

if nlab == 1
        pr = ['asgadj.x par=' label '_par.par'];
else
        pr  = [mpirun_Local ' asgadj.x par=' label '_par.par'];
end
[status,result]=system(pr);
if status~=0
        error(result);
end


% if nlab == 1 || isrand ~= 0 || ismulti ~= 0 || isphase ~= 0 || isromb ~= 0
%     pr= ['asgadj.x par=' label '_par.par'];    % for single processor
%     [status,result]=system(pr)
% else
%     pr = ['mpirun -np ' num2str(nlab) ' asgadj.x par=' label '_par.par'];   %for mpi
%     [status,result]=system(pr)
% end


end
