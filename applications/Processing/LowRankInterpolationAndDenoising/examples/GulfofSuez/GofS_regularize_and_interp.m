% This script performs the regularization and/or interpolation on Gulf of Suez data set.
% Define the input data path and results directory in the script setpath.m
%
% It is not advisable to run this in serial. In parallel mode, 
% when use 4 nodes and 3 workers (for e.g., torque4x3 configuration), this particular 
% experiment will take approximately 30 hours.
%
% In this example, we assume data to be at irregular grid along sources and
% at regular grid along receivers and time.

%% set the input and output directory
setpath
%% Set environmnet and input parameters
options.expdir     = datareg; % define the experiment directory of input data
options.resultdir  = resultsdir; % define the output directory
options.result     = 'GulfofSuez_reg_interp'; % output file name                           
options.nrow       = 356;   % number of sources
options.ncol       = 356;   % number of receivers
options.ntime      = 1024;  % number of time samples
options.d          = [0.004 12.5 12.5];   % sampling interval for time-receiver-source axis
options.flag       = 1;     % input data has odd or even number of source/receiver
options.display    = 1;     % display the results
options.parallel   = sign(parpool_size);  % 0 = run low-rank algorithm in serial mode
                                          % 1 = run low-rank algorithm in parallel mode
options.interp     = 1;     % testing interpolation only , no regularization
options.subsamp    = 0;     % 1 = testing with missing data.
                            % 0 = testing with fully sampled data.
options.algo       = 2;     % 1 = when performing interpolation and denoising
                            % 2 = when performing regularization also
% display the Time,Source,Receiver and frequency index to display the data
if options.display ==1 
   options.idx     = 250;   % define the index of time or frequency axis
   options.ridx    = 150;   % define the index of receiver axis
   options.sidx    = 150;   % define the index of source axis
end
options.rank       = floor(linspace(20,50,floor(options.ntime/2)+1)); % rank can be a single number or vector.
                            % if rank is a vector then it should be equal
                            % to the length of fequency axis along which
                            % regularization and interplation has to be
                            % perfomed                           
options.GoS        = 1;     % 1 = Gulf of Suez data set 
                            % 0 = User defined input data set
if options.GoS == 1
options.factor     = 2;     % define what percent of data to be removed for synthetic test case (2 = 50%)
options.noise      = 0.1;   % define what percent of synthetic data to be replaced by noise (0.1 = 10%)
options.noiselevel = 1e6;   % define the amplitude of noisy traces replacing the GoS data. If amplitude value of noise is very high then use
                            % options.penalty = 1 else use options.penalty = 0.
end
options.penalty    = 0;     % 0 = use least-squares penalty (when amplitude of noise is in the range of data samplitude)
                            % 1 = use students't penalty (when amplitude of noise is
                            % very high compare to data amplitude)
options.sigmafact  = 0.001; % sigma = sigmafact*norm(input_data,2); (error throedshold to solve || Ax-B || <= sigma)
options.initfact   = 1;     % This factor make the amplitude of initial guess in the range of the amplitude
                            % present in the data 
if options.penalty ==1
options.nu         = 8e4;   % The degrees of freedom parameter that controls the sensitivity of students't penalty to outlier noise.
end

% define irregular input grid along sources and regular output grid.Define
% these parameters when testing on your own data set. The imput and output gris position 
% should be a vector where input irregular grid vector should be a size
% equal to number of sources and output grid vector should be a size equal
% to the desired finer source grid. 
if options.interp == 1
    if options.GoS==1
        options.regpos   = load([options.expdir 'regpos.mat']); % periodic sampled grid
        options.irregpos = load([options.expdir 'irregpos.mat']); % unstructured grid 
    else
        options.regpos   = [];
        options.irregpos = [];
    end
end
%% parameters for SPGL1                            
options.iteration  = 150;
options.optTol     = 1e-5;
options.bpTol      = 1e-5;
options.decTol     = 1e-4;

%% Define data file (Gulf of Suez or User-defined)
% The input data file should be in seismic unix format. Time axis should be
% the first coordinate.

if options.GoS == 1   %  Read the Gulf of Suez data set
    options.datafile = [options.expdir 'SuezShots125_355shots_Irreg.mat']; % Run Sconstruct in your experiment directory to get fully sampled Gulf of Suez data set at irregular grid
else
    options.datafile = '../..';  % define here the real data set file (Seismic Unix-Format)
end

% run interpolation and/ or denoising
runinterp(options);

