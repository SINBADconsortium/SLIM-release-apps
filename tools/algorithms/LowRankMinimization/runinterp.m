function [] = runinterp(options)

if ~exist(options.expdir,'dir')
    mkdir(options.expdir);
end
cd(options.expdir);

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
rank           = options.rank;
para.penalty   = options.penalty;
sigmafact      = options.sigmafact;
initfact       = options.initfact;

if para.penalty ==1
    para.nu        = options.nu;
end

para.odd       = options.flag;
para.GoS       = options.GoS;
para.parallel  = options.parallel;

if para.GoS == 1
    para.factor    = options.factor;
    para.noise     = options.noise;
    para.noiselevel= options.noiselevel;
end

nrow           = options.nrow;
ncol           = options.ncol;
ntime          = options.ntime;
para.subsamp   = options.subsamp;
para.interp    = options.interp;

if para.interp==1
    if options.GoS==1
        para.regpos    = options.regpos.pos;
        para.irregpos  = options.irregpos.newpos;
    else
        para.regpos    = options.regpos;
        para.irregpos  = options.irregpos;
    end
end

para.algo     = options.algo;
%% Read Input data
if para.interp==1
    D              = load(options.datafile);
    D              = reshape(D.Dirreg,ntime,nrow,ncol);
else
    D              = ReadSuFast(options.datafile);
    D              = reshape(D,ntime,nrow,ncol);
end
%% Make the test data set 
if para.GoS == 1 
    if para.subsamp==1
        [D] = syndata(D,para);
    end
end
%% Tranform data into frequency domain
[D,K]     = timetofreq(D,ntime,nrow,ncol);
para.nf   = K(1);
para.nrow = K(2);
para.ncol = K(3);

if length(rank)==1
    para.rank = 0;
else
    para.rank = 1;
end
%% generate the unstrcutred operator at regular (output) and irregular (input) grid
if para.interp==1
    if para.GoS==1
        para.regout = sort(para.regpos);
        para.regout = para.regout - min(para.regout);
        para.regout = para.regout ./ max(para.regout);
        para.regout = para.regout - 0.5;
        
        para.regin = sort(para.irregpos);
        para.regin = para.regin - min(para.regin);
        para.regin = para.regin ./ max(para.regin);
        para.regin = para.regin - 0.5;
        
        para.Freg = opKron(opNFFT(para.ncol,para.regout),opNFFT(para.nrow,para.regout));
        para.Fireg = opKron(opNFFT(para.ncol,para.regin),opNFFT(para.nrow,para.regout));
    else
        para.regout = sort(para.regpos);
        para.regout = para.regout - min(para.regout);
        para.regout = para.regout ./ max(para.regout);
        para.regout = para.regout - 0.5;
        
        para.regin = sort(para.irregpos);
        para.regin = para.regin - min(para.regin);
        para.regin = para.regin ./ max(para.regin);
        para.regin = para.regin - 0.5;
        
        para.Freg = opKron(opNFFT(para.ncol,para.regout),opNFFT(para.nrow,para.regout));
        para.Fireg = opKron(opNFFT(para.ncol,para.regin),opNFFT(para.nrow,para.regout));
    end
end
%% Perform Interpolation and/ or denoising
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[data] = LowRank_2D(D,para,options,rank,sigmafact,initfact);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Transform Recovered data into the time-domain
data       = reshape(data,para.nrow,para.ncol,para.nf);
[output,K] = freqtotime(data,para,ntime);

%% write results
if ~exist(options.resultdir,'dir')
    mkdir(options.resultdir);
end

cd(options.resultdir)
save([options.result '.mat'],'output');

end
