%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  L1 constrained least-squares Migration with BG compass model(2D) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a basic LSM migration example. In this example we do not estimate
% Source wavelet.
% 
% For this example, you have to use 'nf' workers "parpool(nf)", otherwise it
% will not work; 'nf' stands for number of frequencies used for least-squares
% imaging. 10 in this case;
%
% 
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1));
close all;

%  Set parameters
tic
% create smooth intial model.
load BG_compass_model409x1401

% true pertubation
dm0         = 1./vel - 1./vels;

% wavelet
cf          = 30;   % central frequency
dt          = .006; % time interval
tt          = 3;    % total sampling time
[wav,taxis] = wvlet(cf,dt,tt);

% frequency range for imaging
fmin = 20;fmax = 50;
% compute source wavelet in frequency domain
[fwav faxis]= fftreal(wav,taxis);
idex        = find(faxis>fmin & faxis<fmax);
fwav        = fwav(idex);
faxis       = faxis(idex);

% imaging parameters
Grid_length = 10;
[nz,nx]     = size(vel);                        % model size
xd          = Grid_length*(nx-1);               % Real distance between horizon grid points
zd          = Grid_length*(nz-1);               % Real distance between depth grid points
nx_border   = floor(0.1*nx);                    % boundary grid points in x
nz_border   = floor(0.1*nz);                    % boundary grid points in z
nx_tot      = nx+2*nx_border;
nz_tot      = nz+2*nz_border;

nshot       = round(nx./2);                     % # of shots
sp          = round(linspace(1,nx,nshot));      % shot locations
sourcet     = 2;                                % source type, 1 sim source; 2 sequential shot.
nsim        = 3;                                % # of sim-shots for each update
nr          = round(nx./2);                     % # of receivers
nr_dis      = 2 * Grid_length;                  % receiver interval              

offset      = 1;                                % marrine data constriction; 1 yes, 2 no.
maxoffst    = 3000;                             % meters
minoffst    = 100;                              % meters

nf          = 10;                               % # of frequencies used in imaging.

rec_dep     = 6;                                % receiver depth
rdp         = 0;                                % right depth precondition, 1 yes; 0 no.
ldp         = 0;                                % left depth precondition, 1 yes; 0 no.
snapshot    = 0;                                % save snapshot, 1 yes; 0 no.
ivc         = 1;                                % inversion crime, 1 yes; 0 no.
mlinm       = 50;                               % max linearization for the imaging problem
l1iter      = 20;                               % max iterations for the l1 solver
l1mode      = 6;                                % 
name    = ['LSM_BGcompass_model2d_L1_WR'];      % name of result

% start a loop here

% results and updates container
results     = zeros(nz,nx,mlinm);
updates     = zeros(nz,nx,mlinm);
dmc         = [];
res         = [];

%---------- Imaging loop start here -------------
dmw_c = [];  % initial guess of update

% number of frequencies
numworkers = parpool_size;
assert(numworkers==nf,...
    'Wrong number of workers in parallel pool.\n It should be %d but is %d.',nf,numworkers)
idex    = randperm(length(faxis));
freq    = faxis(idex(1:nf));
famp    = fwav(idex(1:nf));


for linm = 1:mlinm
	disp(['--- calculating the ',num2str(linm),'update ---'])
	disp(' ')
	
	
	rate    = 1 - nsim./nshot; % subsampling rate
	% compute source term
	if sourcet == 1
		% simultaneous shot source Terms
		[b,nsr] = brm_Border(xd,zd,nx,nz,nx_border,nz_border,nshot,rate,0); 
		b       = sparse(b);
		if  linm == 1,name = [name,'Sim'];end
	elseif  sourcet == 2
		%  single Shot Source Terms
		nsr     = round(nshot * (1 - rate));
		index   = randperm(length(sp));
		sx      = sp(sort(index(1:nsr)));
		sz      = 1;  
		b       = SIG_Source(xd,zd,nx,nz,2,sx,sz,nx_border,nz_border);
		if linm == 1,name = [name,'Seq'];end
		% offset restiction operator, only if sequential shots are used
		if offset
			ofstind = offsetmask(nr,nr_dis,sx,minoffst,maxoffst);
			ofstR   = opColumnRestrict(1,nr*nsr,ofstind);
			ofstR   = opKron(opDirac(nf),ofstR);
			clear ofstind
			if  linm == 1,name = [name,'Offset'];end
		end 
	end
	
	
	% receiver restriction operator 
	du_nl   = zeros(nx,nshot,nf);
	idex    = (nx_border*nz_tot+nz_border +4):2 * nz_tot:(nx_border*nz_tot+nx*nz_tot);
	nr      = length(idex);
	R       = opRestriction(nx_tot*nz_tot,idex);

	spmd
		ifreq = labindex;
		% forward modeling with hard model
	    f            = freq(ifreq);
	    A            = opHelm2D9pt(vel,xd,zd,nx_border,nz_border,f);
	    uh   = famp(ifreq).*full(A\b);
	    uh   = R * uh;
		% forward modeling with smooth model
	    A            = opHelm2D9pt(vels,xd,zd,nx_border,nz_border,f);
	    u     = famp(ifreq).*full(A\b);
	    du_nl = uh  - R* u;
		u     = codistributed.build(u,codistributor1d(3,[ones(1,numlabs)],[size(u) numlabs]),'noCommunication');
		du_nl = codistributed.build(du_nl,codistributor1d(3,[ones(1,numlabs)],[size(du_nl) numlabs]),'noCommunication');	
	end		

	clear b
	clear uh
	clear temp
	if  linm == 1,name = [name,num2str(nsr),'Shots',num2str(nf),'Fs'];end
		

	du_nl  = gather(du_nl);
    % wavefield difference
	du_nl  = du_nl(:);
	
	% define the operator (Jacobian)
	Mig         = opMigSpmd(u,R,freq,nx,nz,nx_border,nz_border,xd,zd,vels);
	
	
	if offset
		du_nl  = ofstR * du_nl;
		Mig = ofstR * Mig;
	end
	
	
	if ivc
		du_nl  = Mig * dm0(:);
		if  linm == 1,name = [name,'Ivc'];end
	end
       % save the result
	if snapshot
		name        = ['../Results/AdMigsnapshotwithrenewal30HZ',num2str(linm)];
		S           = opSnapshot(nz*nx,name);
		Mig         = Mig*S;
	end
	
	% depth precondition
	if rdp
		rdep= sqrt(0:Grid_length:(nz-1)*Grid_length);
		D   = opKron(opEye(nx),opMatrix(spdiags(rdep(:),0,nz,nz)));
		Mig =  Mig * D;
		if linm == 1,name = [name,'RDepPct'];end
	end
	if ldp
		Fdiag = opDiag(freq);
		ldep  = opKron(Fdiag,opEye(nr * nshot));
		Mig   = ldep * Mig;
		du_nl = ldep * du_nl;
		if linm == 1,name = [name,'LDepPct'];end
	end

	% nonlinearized l1 recovery
	[dmw_c,iter_num,C,info] = l1recovery(Mig,du_nl,nx,nz,l1mode,l1iter,dmw_c);
	
	dmw_nl      = C*dmw_c;
	dmw_nl      = reshape(real(dmw_nl),nz,nx);
	
	
	
	
	
	% save results and updates
	results(:,:,linm) = vels;
	updates(:,:,linm) = dmw_nl;
	imagesc(dmw_nl);pause(.1);
	% dmc(:,end + 1)     = dmw_c;
	dmc                 = dmw_c;
	res{linm}       = info;
	
	% save the result
	eval(['save ../Results/',name,' results updates res dmc']);
end


time        = toc;
eval(['save ../Results/',name,' results updates time res dmc']);
