function [] = lin_wrparam(model,label,arg,fname_buoy,fname_bulk,fname_dbulk,fname_dbuoy)

%% Default setting

order  = 2;                % Accurate order in time
cfl    = 0.4;              % CFL stable condition number
cmax   = 4.5;              % max velocity
cmin   = 1;                % min velocity
dmax   = 6;                % max density
dmin   = 0.5;              % min density
fpeak  = model.f0/1000;    % center frequncy
max_step = 0;


% Model information
bulkonly     = 1;             % 1 means only use bulk modules perturbation
bulkmod      = fname_bulk;    % bulk modules
buoyancy     = fname_buoy;    % buoyancy
dbulkmod     = fname_dbulk;   % perturbation of bulk modules
dbuoyancy    = fname_dbuoy;   % perturbation of buoyancy
window_grid  = fname_bulk;    % window grid
window_width = 0;             % window width

datafile     = [label '_data.su'];     % output data file
hdrfile = [label '_hdr.su'];  % header file

dump_lda  = 0;
dump_ldc  = 0;
dump_pi   = 0;
dump_term = 0;

% Parallel parameters
mpi_np1 = 1;    % domain decomposition for 1st direction
mpi_np2 = 1;    % domain decomposition for 2nd direction
mpi_np3 = 1;    % domain decomposition for 3rd direction
partask = 1;    % number of shots simulate parallelly

% pml setting
nd      = length(model.n);
if nd > 2
    nl3 = 0.5;   % y - neg
    nr3 = 0.5;   % y - pos
end
nl1 = 0.5;      % z - neg
nl2 = 0.5;      % z - pos
nr1 = 0.5;      % x - neg
nr2 = 0.5;      % x - pos

% source setting
refamp  = 1;       % reference amplitude
refdist = 1000;    % reference distance
sampord = 1;
pmlampl = 100.0;
scheme_npml = 24;  % accurate order in pml
scheme_phys = 24;  % accurate order in phys
srctype = 'point'; % source type
if isfield(model,'t0')
    source = [label '_src.su'];
end

% read input setting
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


