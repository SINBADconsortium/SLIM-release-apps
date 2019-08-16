function [] = fwd_wrparam(model,label,arg,fname_buoy,fname_bulk)
%% Default setting
order = 2;            %  accurate order for time
cfl   = 0.4;          %  cfl condition number
cmin  = 1.0;          %  min velocity - checked
cmax  = 4.5;          %  max velocity - checked
dmin  = 0.5;          %  min density - checked
dmax  = 6.0 ;         %  max density - checked
max_step = 0;         %  1 = set adaptively, 0 = use standard cfl from cmax
fpeak = model.f0/1000 ;       %nominal central frequency

% Model information
bulkmod  = fname_bulk;
buoyancy = fname_buoy;
mpi_np1  = 1 ;     %n_doms along axis 1
mpi_np2  = 1 ;     %n_doms along axis 2
mpi_np3  = 1 ;     %n_doms along axis 3
partask  = 1 ;     %task parallelization
sampord  = 1 ;     % sampling order

% source information
refdist  = 1000.0; %calibration distance (3D)
refamp  = 1.0;     %target amplitude at calibration distance (3D)
srctype = 'point'; %type - point or array (array requires wavelet)
if isfield(model,'t0')
    source = [label '_src.su'];
end
nd  = length(model.n);
if nd > 2
    nl3 = 0.5;
    nr3 = 0.5;
end
nl1 = 0.5;        % z - neg
nr1 = 0.5 ;       % z - pos
nl2 = 0.5 ;       % x - neg
nr2 = 0.5 ;       % x - pos
pmlampl    = 100.0;
%datafile   = 'data.su';    %output data file
datafile   = [label '_data.su'];
hdrfile    = [label '_hdr.su'];
printact   = 0;
dump_pi    = 0;
dump_lda   = 0;
dump_ldc   = 0 ;
dump_term  = 0 ;
dump_pars  = 0  ;
dump_steps = 0   ;

%%aux. params
%if nargin > 2
%    eval(arg);
%end
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

