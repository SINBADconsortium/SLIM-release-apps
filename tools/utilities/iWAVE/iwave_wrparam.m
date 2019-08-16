function [] = iwave_wrparam(model,label,arg)
% iwave_wrparam(model,label,arg)
% scheme_phys = 24         %scheme (22, 24, 44 - 1D only)
% cfl = 0.4
% cmin = 1.0
% cmax = 4.5
% fpeak = 0.015      %central frequency
% %Model info:
% velocity = vp.rsf
% density = rho.rsf
% %Source info:
% srctype = point
% sampord = 1          %sampling order
% phase = zerophase  %wavelet phase
% refdist = 100.0     %reference distance
% refamp = 1.0        %reference pressure
% %PML
% nl1 = 0.0         %z - neg
% nr1 = 0.5         %z - pos
% nl2 = 0.5         %x - neg
% nr2 = 0.5         %x - pos
% nl3 = 0.5         %y - neg
% nr3 = 0.5         %y - pos
% scheme_npml = 24          %scheme for PML (22, 24)
% %Trace info:
% hdrfile  = [label '_hdr.su'];
% datafile = [label '_data.su'];
% %MPI info:
% mpi_np1 = 1           %n_doms along axis 1
% mpi_np2 = 1           %n_doms along axis 2
% mpi_np3 = 1           %n_doms along axis 3
% partask = 0;          %shot parallelism
% %Output info:
% dump_pi = 1           %dump parallel/dom. decomp info
% dump_lda = 1          % dump grid data for allocated arrays
% dump_ldc = 1           %dump grid data for computational arrays
% dump_term = 1          % dump terminator data

nd = length(model.n);

%% default parameters
scheme_phys = 24 ;        %scheme (22, 24, 44 - 1D only)
cfl = 0.4;
cmin = 1.0;
cmax = 6;
dmin = 0.5;
dmax = 3;
fpeak = model.f0*1e-3 ;     %central frequency
%Model info:
buoyancy= [label '_buoy.rsf'];
bulkmod = [label '_bulk.rsf'];
%Source info:
srctype = 'point';
sampord = 1    ;      %sampling order
phase = 'zerophase' ; %wavelet phase
refdist = 100.0;     %reference distance
refamp = 1.0  ;      %reference pressure
%PML
nl1 = .5 ;        %z - neg
nr1 = .5 ;        %z - pos
nl2 = .5  ;       %x - neg
nr2 = .5  ;       %x - pos
if nd>2
    nl3 = .5 ;       %y - neg
    nr3 = .5  ;       %y - pos
end
scheme_npml = 24 ;         %scheme for PML (22, 24)
%Trace info:
hdrfile  = [label '_hdr.su'];
datafile = [label '_data.su'];
%MPI info:
mpi_np1 = 1 ;          %n_doms along axis 1
mpi_np2 = 1 ;          %n_doms along axis 2
mpi_np3 = 1  ;         %n_doms along axis 3
partask = 1;
%Output info:
dump_pi = 1   ;        %dump parallel/dom. decomp info
dump_lda = 1  ;        % dump grid data for allocated arrays
dump_ldc = 1  ;         %dump grid data for computational arrays
dump_term = 1  ;        % dump terminator data

%%aux. params
if nargin > 2
    eval(arg);
end
clear arg model nd;
S = whos;

%% write
fid = fopen([label '_par.par'],'w');
for k = 1:length(S)
    fprintf(fid,[S(k).name ' = %s\n'],num2str(eval(S(k).name)));
end
fclose(fid);

%% job script
fid = fopen([label '_job.sh'],'w');
np  = mpi_np1*mpi_np2*mpi_np3*max(1,partask);
if np > 1
    nnodes = ceil(np/2);
    ncores = 2;
    if nnodes==1
        fprintf(fid,'#PBS -q smp\n');
    else
        fprintf(fid,'#PBS -q mpi\n');
    end
    fprintf(fid,'#PBS -l nodes=%d:ppn=%d\n',nnodes,ncores);
    fprintf(fid,'#PBS -l walltime=20:00:00\n');
    fprintf(fid,'#PBS -j oe\n');
    fprintf(fid,'#PBS -N iWAVE_%s\n',label);
    fprintf(fid,'module load iWAVE/RSF/10133\n');
    fprintf(fid,'cd %s\n',[pwd]);
    % fprintf(fid,'. ./%s_mkhdr.sh\n',label);
    fprintf(fid,'mpirun -np %d asg.x par=%s\n',np,[label '_par.par']);
else
    fprintf(fid,'#PBS -q sngl\n');
    fprintf(fid,'#PBS -l nodes=1:ppn=1\n');
    fprintf(fid,'#PBS -l walltime=05:00:00\n');
    fprintf(fid,'#PBS -j oe\n');
    fprintf(fid,'#PBS -N iWAVE_%s\n',label);
    %fprintf(fid,'module load iWAVE/srl/253\n');
    fprintf(fid,'module load iWAVE/RSF/10133\n');
    fprintf(fid,'cd %s\n',[pwd]);
    %fprintf(fid,'. ./%s_mkhdr.sh\n',label);
    fprintf(fid,'asg.x par=%s\n',[label '_par.par']);
end
system(['chmod a+rwx ' label '_job.sh']);

fclose(fid);
