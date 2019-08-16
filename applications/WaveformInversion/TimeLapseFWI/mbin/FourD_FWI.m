function  [results] = FourD_FWI(vel_initial,Dbase,Dmon,model_base,model_monitor,inv_params,expdir)

% 	This funtion performs the full-waveform inversion of an observed baseline and monitor data set.
%     It uses the modified Gauss-Newton inversion approach where each gauss-newton subproblem is solved with 
%     the l1 solver (SPGL1) developed by Ewout van den Berg and Michael P. Friedlander. 
% 	
% 	
%     usage:
%     [results] = TimeLapse_FWI(vel_initial,Dbase,Dmon,model_base,model_monitor,inv_params,options)
% 
%     
%     Input:
%     vel_initial: same starting model for FWI on baseline and monitor, 2D matrix;
%                  vertical axis should be depth in meters. Unit is m/s
%     Dbase: observed baseline data, 3D data cube, sorted as nrec*nsrc*nfreq
%            i.e. number of receivers, sources, frequencies
%     Dmon: observed monitor data, 3D data cube, sorted as nrec*nsrc*nfreq
%            i.e. number of receivers, sources, frequencies
%     model_base: data struct with the modeling parameters for baseline
%     model_monitor: data struct with the modeling parameters for monitor 
%     inv_params: jdata struct with joint inversion parameters used in the main inversion operator
%               beta: the weight on the l1 norm of the common part, default is 1
%               nproc: number of frequencies to use in each subproblem
%               nsim: number of simultaneous shots to use
%               nexp: number of GN subproblems to be solved
%               maxiter: maximum number of iterations for each subproblem
%     expdir: location of data, must be defined in the script
%     
%     
%     Author: Felix Oghenekohwo
%          Seismic Laboratory for Imaging and Modeling
%          Department of Earth, Ocean and Atmospheric Sciences
%          The University of British Columbia
% 	Date: February 2015.
	
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.	
	
tic
%%===================initializing=======================
% setup a random seed
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1));


%% load velocity, observed data and model parameters
vel_initial = 1e6./(vel_initial(:)).^2;
m0 = vel_initial(:);
m1f = m0;m2f = m0;
mbf = m0;mmf = m0;
model = model_base;

nrec = length(model.xrec);
nsrc = length(model.xsrc);
nfreq = length(model.freq);
Q = speye(length(model_base.xsrc));

Dbase = reshape(Dbase,nrec,nsrc,nfreq);
Dmon =  reshape(Dmon,nrec,nsrc,nfreq);

%% Define the sparsifying transform and the 2-by-3 block
nz = size(vel_initial,1);nx = size(vel_initial,2);
C = opCurvelet(nz,nx,max(1,ceil(log2(min(nz,nx)) - 3)),16,1,'ME',0);
N = size(C,1);

beta = inv_params.beta;
Ic = [opStack(beta*C',beta*C') opBlockDiag(C',C')];

%% load inversion parameters
iter = inv_params.maxiter;
opts = spgSetParms('optTol',1e-3, ...
                   'bpTol' , 1e-3 , ...
                   'decTol',5e-2,...
                   'iterations',iter,...
                    'minPareto',ceil(iter/2),...
                    'quitPareto',1);

nsim = inv_params.nsim;
nexp = inv_params.nexp;


%% create the freq. batches for the inversion : same for base/monitor

a    = model.freq(1);
nproc = inv_params.nproc;
deltaf = model.freq(2) - model.freq(1);


b = a+(nproc-1)*deltaf;
i = 1;
while b < model.freq(end)
    batch(:,i) = a:deltaf:b;
    i = i+1;
    a = model.freq(1)+(i-1)*ceil((nproc/2))*deltaf;  
    b = a+(nproc-1)*deltaf;
end

if batch(end) < model.freq(end)
    batch(:,i) = sort(model.freq(end):-deltaf:model.freq(end)+(nproc-1)*-(deltaf));
end   
nbatch = size(batch,2);   
clear i

IBase = zeros(nz*nx,nexp,nbatch);
IMon  = zeros(nz*nx,nexp,nbatch);
JBase = zeros(nz*nx,nexp,nbatch);
JMon  = zeros(nz*nx,nexp,nbatch);

%% Inversion loop starts here


for j = 1:nbatch
    
    batch_freq = batch(:,j);
    [~,index_freq,~] = intersect(model.freq,batch_freq);
    opSubdata = opKron(opRestriction(nfreq,index_freq),opKron(opDirac(nsrc),opDirac(nrec)));
    
    % redraw over frequency baseline
    model_bsub = model_base;
    model_bsub.freq = batch_freq;
    Dbsub = opSubdata*Dbase(:);
    
    
    % redraw over frequenncy monitor
    model_msub = model_monitor;
    model_msub.freq = batch_freq;
    Dmsub = opSubdata*Dmon(:);
    
    for i = 1:nexp
        
        
        disp(['=== Now processing frequency band ',num2str(j), 'for the run', num2str(i), '==='])
        disp(' ')
        
        % restriction
        G = opGaussian(nsim,nsrc);
        RM = opKron(opDirac(length(batch_freq)), opKron(G,opDirac(nrec)));
        
        %% Parallel Inversion
        % baseline imaging operator
        Jb = oppDF_4D(mbf,Q*G',model_bsub,1);
        tic;dDb = F_4D(mbf,Q*G',model_bsub,1);toc;
        % monitor imaging operadDmtor
        Jm = oppDF_4D(mmf,Q*G',model_msub,1);
        tic;dDm = F_4D(mmf,Q*G',model_msub,1);toc;
        % right hand side ( b in Ax = b)
        dDbsub = (RM*Dbsub) - dDb;
        dDmsub = (RM*Dmsub) - dDm;
        tic;
        xestj = spgl1(Jb*C', dDbsub, 0, 0, [], opts);
        toc;
        z1 = real(C'*xestj(1:N)); mbf = mbf+z1;
        tic;
        xestj = spgl1(Jm*C', dDmsub, 0, 0, [], opts);
        toc;
        z2 = real(C'*xestj(1:N)); mmf = mmf+z2;
        
        IBase(:,i,j) = mbf(:);
        IMon (:,i,j) = mmf(:);
        
        
        %save Result_Updates nbatch IBase IMon
        %% Joint Inversion
        % baseline imaging operator
        Jb = oppDF_4D(m1f,Q*G',model_bsub,1);
        tic;dDb = F_4D(m1f,Q*G',model_bsub,1);toc;
        
        % monitor imaging operadDmtor
        Jm = oppDF_4D(m2f,Q*G',model_msub,1);
        tic;dDm = F_4D(m2f,Q*G',model_msub,1);toc;
        
        % right hand side ( b in Ax = b)
        dDbsub = (RM*Dbsub) - dDb;
        dDmsub = (RM*Dmsub) - dDm;
        
        A = opBlockDiag(Jb,Jm)*Ic;
        b = [dDbsub ; dDmsub];
        
        tic;
        xestj = spgl1(A, b, 0, 0, [], opts);
        toc;
        
        z0 = real(C'*beta*xestj(1:N));
        z1 = real(C'*xestj(N+1:2*N));
        z2 = real(C'*xestj(2*N+1:end));

        m1f  = m1f + z0 + z1;
        m2f  = m2f + z0 + z2;
        
        JBase(:,i,j) = m1f(:);
        JMon (:,i,j) = m2f(:);
        
        
        %save Result_Updates nbatch JBase JMon
       
        
    end
    
    save([expdir '/IRS_UPdates.mat'], 'j', 'nbatch', 'IBase', 'IMon')
    save([expdir '/JRM_UPdates.mat'], 'j', 'nbatch', 'JBase', 'JMon')
    
    
end
toc;
results.ind_base = mbf;
results.ind_mon  = mmf;
results.joint_base = m1f;
results.joint_mon  = m2f;


end

