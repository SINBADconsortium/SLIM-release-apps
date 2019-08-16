function [status,result] = inv_wrapper(init_buoyancy,init_bulkmod,model,datafile,params)
%
% Usage: iwavewrapper(v,rho,model)
% v refers to the velocity model;rho refers to the density model
% iwave_wrmodel,iwave_wrdatahdr, iwave_wrsrc, iwave_wrparam

% setting the path
label = 'inv_results';     %Results directory
mkdir(label);
curdir = pwd;
cd(label);

% read velocity [m/s], select subset
[buoy,n,d,o]=rsf_read_all(init_buoyancy);
[bulk,n,d,o]=rsf_read_all(init_bulkmod);
n = size(buoy);

% grid params
model.o = o; model.d = d; model.n = n;

% % matlab scripts to generate input files for iWAVE
iwave_wrmodel(buoy,bulk,model,label);
iwave_wrdatahdr(model,label);
iwave_wrsrc(model,label);
inv_wrparam(model,label,params,init_buoyancy,init_bulkmod,datafile);

spmd
    nlabs = numlabs;                                % function numlabs-Total number of workers operating in parallel on current job
end
nlab=nlabs{1};

%%%%% For job submission %%%%%
%scriptname = [label '_job.sh'];              % job script produced from the iwave scripts
%fprintf('scriptname is: %s',scriptname);
%[status,result]=unix(sprintf('qsub %s',scriptname));
%%%%%%%%%%%%

mkh = ['. ./' label '_mkhdr.sh'];
[status,result]=system(mkh);  % execute the shell script to generate the header file


if nlab == 1
    pr= ['asginv.x par=' label '_par.par'];    % for single processor
    [status,result]=system(pr)
else
    pr = ['mpirun -np ' num2str(nlab) ' asginv.x par=' label '_par.par'];   %for mpi
    [status,result]=system(pr)
end


end
