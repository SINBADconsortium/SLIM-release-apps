%% 3D constant-density acoustic frequency-domain modeling: WaveEquationFD function
%
% This script aims at demonstrating the simplest possible use of the forward
% modeling kernel using WaveEquationFD function. As usual, this script starts 
% cleaning the enviroment and opening the log file
% 
% The new version of the code is rather verbose if the optional parameter |flog|
% is passed to WaveEquationFD; it prints information both on the screen and on 
% the log file. If multiple Matlab workers are available, only the master worker
% prints the information while the other workers ignore it.
% 
% Also the new code is very simple and requires virtually no set up other than
% reading the data itself. The new routine fixes several issues present in the 
% legacy version such as 
% 
% * automatic and abstract use of PML
% * properly rounded (up) grid spacing for the computational grid to satisfy 
%   stability conditions.
% * dynamic computational grid chosen depending on the frequency - the user
%   does not need to manually shrink/stretch the model for lower/higher 
%   frequency.
% * Computed wavefields match analitical wavefields
% * Fixed the sign of the wavefield - legacy code produces -conj(.) of the 
%   expected data
%
vfile='../data/m_true.rsf';
dfile='../data/Dobs.rsf';
flog = fopen('../results/demo_MLCRMN_log','w');

% Read data
plog(flog,'* Reading data... \n');
[model.v,model.nv,model.dv,model.ov] = rsf_read_all(vfile);
[acq.data,acq.nd,acq.dd,acq.od] = rsf_read_all(dfile);
model.unit = 'm/s';

% Create sources for synthetic data
plog(flog,'* ',prod(acq.nd(4:6)),' shots read; drawing receivers grid...\n');
[ acq.xrec,acq.yrec,acq.zrec,...
  acq.xsrc,acq.ysrc,acq.zsrc, acq.freq] = odn2grid(acq.od,acq.dd,acq.nd);
acq.sources  = speye( prod(acq.nd(4:6)) );

% Setting up matrix, preconditioner and solver
plog(flog,'* Setting up preconditioner...\n');
par_solver.maxcy = 1;        % Just because this is an example; for real cases, 
                             % usually you shouldn't set this.
par_solver.name = 'mlcrmn';  % Just to ensure that MLCRMN will be used. In most
                             % cases, the script will automagically choose a 
                             % solver for you
H = helmholtz_solver(model,4,[],par_solver,flog);
% ...and done

% Setting up right hand side; putting a point source in the center
Q = zeros(H.nt);
Q(ceil(H.nt(1)/2),ceil(H.nt(2)/2),ceil(H.pml.x+2)) = 1;
Q = Q(:);

% Solve...
plog(flog,'\n\n* Starting solver...\n');
u = H.solve(Q,0*Q);
% ... solved... (in this example we only run 1 cycle, which means 5 iterations).
% Unset the "maxcy" parameter in par_solver to allow it to run as many cycles as
% necessary.

fclose(flog);
