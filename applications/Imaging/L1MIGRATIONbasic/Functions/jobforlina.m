addpath(genpath('/scratch/slic/xiang/matlab/lixiang/FWI_Spot/'))
addpath(genpath('/scratch/slic/xiang/matlab/lixiang/Lina_code/'))
addpath(genpath('/users/slic/xiang/matlab/CREWESmatlabrelease'))
addpath(genpath('/users/installs/slic/mat_toolbox/sparco-1.2/tools'))
% least square migration of BG model with renewals
RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
clear;  close all;
%  Set parameters
tic
load ./result/FwiBGmodel7Shots29HzSmartBdr
vel1 = results(:,:,end);
clearvars -except vel1

load bg2dmodel.mat
vel         = 1000*vp;
vel1        = lpf(vel,0.1);
clear vp v0 marmhard marmhard256128
vmin  = min(vel(:));
vmax  = max(vel(:));


% interp vel
[X,Y]   = meshgrid(1:size(vel,2),1:size(vel,1));

fg      = .5; 
[X1,Y1] = meshgrid(1:fg:size(vel,2),1:fg:size(vel,1));
vel1    = interp2(X,Y,vel1,X1,Y1,'cubic');
vel     = interp2(X,Y,vel,X1,Y1,'cubic');
vels    = vel1;
clear X Y X1 Y1

% project the initial
index       = find(vel1 < vmin);
vel1(index) = vmin;
index       = find(vel1 > vmax);
vel1(index) = vmax;
dm0         = 1./vel - 1./vel1;

% wavelet
cf          = 30; % hz
dt          = .006; % time interval
tt          = 3;  % total sampling time
[wav,taxis] = wvlet(cf,dt,tt);

% frequency range
fmin = 20;fmax = 50;
% frequency parameter setup
[fwav faxis]= fftrl(wav,taxis);
idex        = find(faxis>fmin & faxis<fmax);
fwav        = fwav(idex);
faxis       = faxis(idex);

Grid_length = 10;
% rate        = 0.95; 
nsim        = 3;  % No of sim-shots for each update
[nz,nx]     = size(vel);
xd          = Grid_length*(nx-1);               %Real distance between horizon grid points
zd          = Grid_length*(nz-1);               %Real distance between depth grid points
nx_border   = floor(0.1*nx);
nz_border   = floor(0.1*nz);
nx_tot      = nx+2*nx_border;
nz_tot      = nz+2*nz_border;

nshot       = round(nx./2);
sp          = round(linspace(1,nx,nshot));
sourcet     = 1; % source type, 1 sim source; 2 sequential shot.

nr          = round(nx./2);
nr_dis      = 2 * Grid_length;

offset      = 0; % marrine data constriction; 1 yes, 2 no.
maxoffst    = 3000; % meters
minoffst    = 100;  % meters



rec_dep     = 6;
rdp         = 0;  % right depth precondition, 1 yes; 0 no.
ldp         = 0;  % left depth precondition, 1 yes; 0 no.
snapshot    = 1;  % save snapshot, 1 yes; 0 no.
ivc         = 1;  % inversion crime, 1 yes; 0 no.
miter       = 1;  % max iterations of the FWI Loop
mlinm       = 10; % max linearization for one FWI loop
l1iter      = 300; % max iterations for the l1 solver
l1mode      = 13;  % solve l1 in which domain
name    = ['Jobsforlina1']; % name of result
%  single Shot Source Terms
sx = 128;   sz = 5;    
% [b,ixyshot] = rhs2d(xd,zd,nx,nz,nshot);
% [b] = SIG_Source(xd,zd,nx,nz,2,sx,sz,nx_border,nz_border);      
% nshot = 1;

% start a loop here

% % get rid of the low freqs
% lowcut          = 4;
% faxis(1:lowcut) = [];
% fwav(1:lowcut)  = [];


% results and updates container
results     = vel1;
updates     = [];
res         = [];
dmc         = [];

%---------- FWI loop start here -------------
for iter = 1:miter
	dmw_c = [];  % initial guess of update
	disp(['=== Now Caculating the frequency band ',num2str(iter),'==='])
	disp(' ')
	
	% number of frequencies
	nf      = 10;
	idex    = randperm(length(faxis));
	freq    = faxis(idex(1:10));
	famp    = fwav(idex(1:10));
	% freq    = faxis((nf./2)*(iter-1)+1:(nf./2)*(iter+1));
	% famp    = fwav((nf./2)*(iter-1)+1:(nf./2)*(iter+1));
	rate    = 1 - nsim./nshot;
	if sourcet == 1
		% simultaneous shot source Terms
		[b,nsr] = brm_Border(xd,zd,nx,nz,nx_border,nz_border,nshot,rate,0); 
		b       = sparse(b);
		if iter == 1, name = [name,'Sim'];end
	elseif  sourcet == 2
		%  single Shot Source Terms
		nsr     = round(nshot * (1 - rate));
		index   = randperm(length(sp));
		sx      = sp(sort(index(1:nsr)));
		sz      = 1;  
		b       = SIG_Source(xd,zd,nx,nz,2,sx,sz,nx_border,nz_border);
		if iter == 1 && linm == 1,name = [name,'Seq'];end
		if offset
			ofstind = offsetmask(nr,nr_dis,sx,minoffst,maxoffst);
			ofstR   = opColumnRestrict(1,nr*nsr,ofstind);
			ofstR   = opKron(opDirac(nf),ofstR);
			clear ofstind
			if iter == 1 && linm == 1,name = [name,'Offset'];end
		end 
	end
	
	
	
	for linm = 1:mlinm
		disp(['---calculating the ',num2str(linm),'update for frequency band ',num2str(iter),'---'])
		disp(' ')
		% 
		% % % simultaneous shot source Terms
		% nshot   = round(nx./4);                               % Shot position numbers
		% [b,nsr] = brm_Border(xd,zd,nx,nz,nx_border,nz_border,nshot,rate,0); 
		% b       = sparse(b);
		% nshot   = nsr;
		% 	
	
		% forward modeling with hard model
 		du_nl   = zeros(nx,nshot,nf);
		% u       = zeros(nx_tot*nz_tot,nshot,nf);
		idex    = (nx_border*nz_tot+nz_border +4):2 * nz_tot:(nx_border*nz_tot+nx*nz_tot);
		nr      = length(idex);
		R       = opRestriction(nx_tot*nz_tot,idex);

		spmd
			ifreq = labindex
		    f            = freq(ifreq);
		    A            = opHelm2D9pt(vel,xd,zd,nx_border,nz_border,f);
		    uh   = famp(ifreq).*full(A\b);
		    uh   = R * uh;
			% forward modeling with smooth model
		    A            = opHelm2D9pt(vel1,xd,zd,nx_border,nz_border,f);
		    u     = famp(ifreq).*full(A\b);
		    du_nl = uh  - R* u;
			u     = codistributed.build(u,codistributor1d(3,[ones(1,numlabs)],[size(u) numlabs]),'noCommunication');
			du_nl = codistributed.build(du_nl,codistributor1d(3,[ones(1,numlabs)],[size(du_nl) numlabs]),'noCommunication');	
		end		


		clear uh
		clear temp
		if iter == 1 && linm == 1,name = [name,num2str(nsr),'Shots',num2str(nf),'Fs'];end
			
	
		du_nl  = gather(du_nl);
	    % wavefield difference
		du_nl  = du_nl(:);
		
		% define the operator
		Mig         = opMigSpmd(u,R,freq,nx,nz,nx_border,nz_border,xd,zd,vel1);
		
		
		if offset
			du_nl  = ofstR * du_nl;
			Mig = ofstR * Mig;
		end
		
		
		if ivc
			du_nl  = Mig * dm0(:);
			if iter == 1 && linm == 1,name = [name,'Ivc'];end
		end
        % save the result
		if snapshot
			sname        = ['./result/jobforlina',num2str(linm)];
			S           = opSnapshot(nz*nx,sname);
			Mig         = Mig*S;
		end
		
		% depth precondition
		if rdp
			rdep= sqrt(0:Grid_length:(nz-1)*Grid_length);
			D   = opKron(opEye(nx),opMatrix(spdiags(rdep(:),0,nz,nz)));
			Mig =  Mig * D;
			if iter == 1 && linm == 1,name = [name,'RDepPct'];end
		end
		if ldp
			Fdiag = opDiag(freq);
			ldep  = opKron(Fdiag,opEye(nr * nshot));
			Mig   = ldep * Mig;
			du_nl = ldep * du_nl;
			if iter == 1 && linm == 1,name = [name,'LDepPct'];end
		end

		% nonlinearized l1 recovery
		[dmw_c,iter_num,C,info] = l1recovery(Mig,du_nl,nx,nz,l1mode,l1iter,dmw_c);
		
		dmw_nl      = C*dmw_c;
		dmw_nl      = reshape(real(dmw_nl),nz,nx);
				
		% save results and updates
		results(:,:,end+1) = vel1;
		updates(:,:,end+1) = dmw_nl;
		% dmc(:,end + 1)   = dmw_c;
		dmc                = dmw_c;
		res{end + 1}       = info;
		
		% save the result
		eval(['save ./result/',name,' results updates res dmc']);
	end
end

time        = toc;
eval(['save ./result/',name,' results updates time res dmc']);