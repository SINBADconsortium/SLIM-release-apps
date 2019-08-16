function [] = adj_wrparam(model,label,arg,fname_buoy,fname_bulk,datafile)
nd = length(model.n);

%% default parameters
order       = 2;
scheme_phys = 24 ;        %scheme (22, 24, 44 - 1D only)
cfl         = 0.4;
cmin        = 1.0;
cmax        = 4.5;
dmin        = 0.5;
dmax        = 6;
fpeak       = model.f0/1000;;     %central frequency
max_step    = 0;


%Model info:
bulkonly      = 1;
bulkmod       = fname_bulk;
buoyancy      = fname_buoy;
dbulkmod      = fname_bulk;
dbuoyancy     = fname_buoy;
mbulkmod      = [label '_cammbulk.rsf'];
mbuoyancy     = [label '_cammbuoy.rsf'];
window_grid   = fname_bulk;
window_width  = 0;

%Source info:
srctype = 'point';
if isfield(model,'t0')
    source = [label '_src.su'];
end
sampord = 1    ;      %sampling order
refdist = 1000.0;     %reference distance
refamp = 1.0  ;      %reference pressure


%PML
nl1 = .5 ;        %z - neg
nr1 = .5 ;        %z - pos
nl2 = .5 ;        %x - neg
nr2 = .5 ;        %x - pos
if nd>2
    nl3 = .5;        %y - neg
    nr3 = .5;        %y - pos
end
scheme_npml = 24 ;         %scheme for PML (22, 24)
pmlampl = 100.0;

%Trace info:
datafile = datafile;
hdrfile  = datafile;
%MPI info:
mpi_np1 = 1 ;          %n_doms along axis 1
mpi_np2 = 1 ;          %n_doms along axis 2
mpi_np3 = 1 ;          %n_doms along axis 3
partask = 1;

%Checkpoint info
nsnaps = 20;          % number of live checkpoints
maxbuffernum = 20;    %     max num of incore buffers

%Output info:
dump_pi    = 0   ;        %dump parallel/dom. decomp info
dump_lda   = 0   ;        % dump grid data for allocated arrays
dump_ldc   = 0   ;        %dump grid data for computational arrays
dump_term  = 0   ;        % dump terminator data
dump_steps = 0   ;
printact   = 0   ;
dump_pars  = 0   ;


if nargin > 2
    for i = 1:length(arg)
        eval(arg{i});
    end
end
clear arg model nd;
S = whos;

%% write
fid = fopen([label '_par.par'],'w');
for k = 1:length(S)
    fprintf(fid,[S(k).name ' = %s\n'],num2str(eval(S(k).name)));
end
fclose(fid);


