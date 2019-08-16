
%% load data. 
% please download Chevron 2014 SEG workshop benchmark dataset 
% <https://s3.amazonaws.com/open.source.geoscience/open_data/seg_workshop_fwi_2014/seg_workshop_fwi_2014.html>
% please put the SEG file in the 'data' folder named as below
filename1	= '../data/SEG14.Pisoelastic.segy'; % Chevron released SEGY dataset.


%% hardware requirement
% To run the script, you need at least 40 cores. each core needs a least 2GB memory. 


%% open matlab pool (160 sugguested, you can also use 80,40 according to the number of cores).


%% reading some infomation from the headfile
nfile		= 1;
[h]			= ReadSegyHeader(filename1);
ntrace		= h.ntraces;
nr			= h.DataTracePerEnsemble;
taxis		= h.time;
dt			= taxis(2) - taxis(1);
taxis		= 0:dt:8;
ns			= ntrace./nr;

load ../data/rec_src_pos

% choose subset of the shot
x1 = 0;
x2 = 46000;

% load model
load ../data/model_fwi_mask
idx  = find(mmx>x1 & mmx<x2);
% mmx      = mmx(idx);
% vel_freq = vel_freq(:,idx);
[nz,nx] = size(vel_init);
nbz = 50;
nbx = 50;
Pab = opExtend2(nx-2*nbx,nz-nbz,nbx,nbz,1,1);


model.vel = vel_init(:);
model.den = 1;
% find sea bottom
% sea_bottom = find(vel_freq <1515);
% den 			= .31 .*vel_freq.^.25;
% den(sea_bottom,:) 		= 1;

shot_start = 1;
shot_all   = 1:1599;
shot_idx   = shot_start + shot_all;
NS		   = length(shot_idx);
ridx	   = nr*(shot_idx(1)-1)+1:nr*shot_idx(end);
ridx	   = reshape(ridx,nr,NS);

%
ns_sub   = 160 ; % 50 shots for each subproblem

fwdpara.fmode = 1;

% turning wave modeling
% nfb = [5,8,12,15,20,25];% 3,5,8,9,10
nfb = [5,7,9,11,13]; 
ewc  = [100,80,60,50,40,30,20,10];
ewr  = ones(nr,1);
ewr(1:10) = 0;ewr(end-10:end) = 0;
nsub = 3;
l2iter= 6;

% load ../data/mask
% clear Tmask
% Tmask = Rmask';
% nt 	  = size(Tmask,1);
% --- make time mask ----

for mmm = 1:length(nfb)
	tidx = (100:2700);
	% P1   = opKron(opSmoothmod(nr,round(nfb(mmm)./50)),opSmoothmod(length(tidx),nfb(mmm)));
	% P_lp = oppKron2Lo(opDiag(ones(ns_sub,1)),P1);
	

	% edge effect operator
	ew = ones(length(tidx),1);
	ew(1:ewc(mmm)) =0;
	ew(end-ewc(mmm):end) = 0;
	
	Pe =  opKron(opDiag(ewr),opDiag(ew));
	Pe = oppKron2Lo(opDiag(ones(ns_sub,1)),Pe);
	
	
	for nnn = 1:nsub

		ns_sub_idx = randperm(NS);
		ns_sub_idx = sort(ns_sub_idx(1:ns_sub));


		fwdpara.slx = slx_all(ns_sub_idx+shot_start);
		fwdpara.slz = slz_all(ns_sub_idx+shot_start);
		fwdpara.rtype = 'marine';
		fwdpara.rlx = rlx_all(:,ns_sub_idx+shot_start);

		fwdpara.rlz = rlz_all(:,ns_sub_idx+shot_start);

		ridx1 = vec(ridx(:,ns_sub_idx));

		D = XiangLi_ReadSegy(filename1,'traces',ridx1);

		outb_idx 			  = find(fwdpara.rlx>x2-1000);
		fwdpara.rlx(outb_idx) = x2 - 60;
		D(:,outb_idx) 		 = 0;

		% setup for modeling
		fwdpara.mmy = 1;


		fwdpara.stype = 'seq_same';
		fwdpara.sly   = ones(size(fwdpara.slx));
		fwdpara.rly   = ones(size(fwdpara.rlx));

		fwdpara.saved 	= 1;
		fwdpara.shot_id = 1;
		
		% prepare source wavelet
		load ../data/chevron_wavelet
		taxis 	= -1:dt:8;
		Pt 		= opLInterp1D(t,taxis);
		src 	= Pt * wt;

		D	= [zeros(250,nr*ns_sub);D];
		taxis = 0:dt:9;

		fwdpara.dt = .002;
		fwdpara.taxis = 0:fwdpara.dt:9;

		Pt = opLInterp1D(taxis,fwdpara.taxis);
		src = Pt * src;
		D   = Pt * D;
		fwdpara.fcent = 10;

		fwdpara.abx = 50;
		fwdpara.abz = 50;
		fwdpara.aby = 50;
		

		fwdpara.free = 1;
		fwdpara.cut_model = 1;
		
		fwdpara.ddcompx = 1;
		fwdpara.ddcompz = 1;
		fwdpara.ddcompy = 1;


		fwdpara.mmx = mmx;
		fwdpara.mmz = mmz;
		fwdpara.mmy = 1;

		fwdpara.chkp_space	= 100;
		fwdpara.chkp_save	= 2;

		fwdpara.disp_snapshot = 0;
		fwdpara.wave_equ = 1;
		fwdpara.v_up_type = 'slowness';
		fwdpara.space_order = 4;


		


		% do modeling
		
		fwdpara.taxis = fwdpara.taxis(tidx);
		nt = length(fwdpara.taxis);
		
		
		D			  = D(tidx,:);
		src			  = src(tidx);
		
		% P_lp = opBandPassFilter(length(fwdpara.taxis),1,fwdpara.taxis,[0 3 nfb(mmm) nfb(mmm)+1],'c');
		P_lp = opBandPassFilter(length(fwdpara.taxis),1,fwdpara.taxis,[0 3 nfb(mmm) 2.*nfb(mmm)-4],'c');
		D_lp = P_lp * D;
		src_lp  = P_lp * src(:);
		src_lp(1:100) = src_lp(1:100).*(linspace(0,1,100).^2)';
		
		[Ur,chkp] = XiangLi_Time(model,src_lp,fwdpara); 
		Ur = reshape(Ur,size(D_lp));
		Ur(:,outb_idx) = 0;
		Ur =  Ur(:);
		Ps			 	= op_XiangLi_Time(model,src_lp,fwdpara,chkp);

		alpa = -1.5e7;
		D_lp = alpa .* D_lp(:);
		
		dU   = D_lp - Ur;
		
		
		
		if mmm==1 &&nnn ==1
			Rmask = zeros(nt,nr);
			tcut1 = linspace(2,6.7,nr);
			tcut2 = linspace(5,5.2,nr);
			for tcut_index = 1:nr;
				tcut_idx = find(fwdpara.taxis>tcut1(tcut_index)&fwdpara.taxis<tcut2(tcut_index));
				Rmask(tcut_idx,tcut_index) = 1;
			end
			Ptcut_smooth = opSmooth(size(Rmask,1),200);
			Rmask = Ptcut_smooth * Rmask;
			
			P_select_data = oppKron2Lo(opDirac(ns_sub),opDiag(Rmask(:)));
			
			P_select_model = opDiag(mask1(:));
			
		end
		
		
		
		
		Plp = oppKron2Lo(opDirac(ns_sub),opKron(opSmooth(nr,10),P_lp));
		
		dU   = Plp * Pe * P_select_data * dU;

		P_all = Plp * Pe * P_select_data *  Ps * P_select_model;

		dm = lsqr_fwi(P_all,dU,l2iter,1); 

		% dm = l1curveletreg2(dm,nz,nx,mmz(end)./1e3,.8,1,'yes');
		
		
		
		model.vel = 1./(1./model.vel + dm);
		
		model.vel     = Pab*Pab'*	model.vel ;
		idx = find(	model.vel  < 1510);
		model.vel (idx) = 1510;
		idx = find(	model.vel  > 5500);
		model.vel (idx) = 5500;		
		% imagesc(reshape(vel_freq,467,1000));pause(1);
		
		
		% vel_freq = reshape(vel_freq,length(fwdpara.mmz),length(fwdpara.mmx));
		
		% water layer
			% model.vel(sea_bottom) = 1510;
		% vel_freq 		= vel_freq(:);
		% den 			= .31 .*vel_freq.^.25;
		% den(sea_bottom) 		= 1;
		
		eval(['save ../results/chevron2_reflec2_',num2str(mmm),'_nsp_',num2str(nnn),' model mmm nnn '])

	end
end
















