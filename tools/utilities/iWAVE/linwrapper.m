function [status,result] = linwrapper(fname_buoy,fname_bulk,fname_dbulk,fname_dbuoy,model,params,options)
%
% Usage: iwavewrapper(v,rho,model)
% v refers to the velocity model;rho refers to the density model
% iwave_wrmodel,iwave_wrdatahdr, iwave_wrsrc, iwave_wrparam

% setting the path
if isfield(options,'label')
    label = [options.label,'_lin'];
else
    label = 'lin_results_temp';
end

if ~exist(label,'dir')
    mkdir(label);
end

%  matlab scripts to generate input files for iWAVE
if isfield(model,'simultsrc')
    ismulti = model.simultsrc;
else
    ismulti = 0;
end

if isfield(model,'romb')
    isromb = model.romb;
else
    isromb = 0;
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

if ismulti == 0 && isphase == 0 && isromb == 0 && isrand == 0
    curdir = pwd;
    cd(label);
    flagReg = 1;
else
    flagReg = 0;
end

% Generate the script to generate the header file
if ismulti == 0 && isphase == 0 && isromb == 0
    iwave_wrdatahdr(model,label);
else
    iwave_wrdatahdr_sim(model,label);
end

% Generate the parameter file
lin_wrparam(model,label,params,fname_buoy,fname_bulk,fname_dbulk,fname_dbuoy);


% Generate the source file if it is needed
if isfield(model,'t0')
    if ismulti == 0 && isphase == 0 && isromb == 0
        iwave_wrsrc(model,label);
    else
        iwave_wrsrc_sim(model,label);
    end
end

% Detect the number of labs
nlab    = options.nlab;


% Make the header file
mkh = ['. ./' label '_mkhdr.sh'];
[status,result]=system(mkh);  % execute the shell script to generate the header file


% call the iwave function
%mpirun_Local = options.MPIRUN_Command;
if isfield(options,'MPICommand')
        mpirun_Local = options.MPICommand;
else
        mpirun_Local = iwave_MPIcommand(options.nlab,flagReg);
end

if nlab == 1
    pr = ['asglin.x par=' label '_par.par'];
else
    pr  = [mpirun_Local ' asglin.x par=' label '_par.par'];
end
[status,result]=system(pr);
if status~=0
    error(result);
end
% if nlab == 1
%     pr= ['asglin.x par=' label '_par.par'];    % for single processor
%     [status,result]=system(pr)
% else
%     pr = ['mpirun -np ' num2str(nlab) ' asglin.x par=' label '_par.par'];   %for mpi
%     [status,result]=system(pr)
% end

end
