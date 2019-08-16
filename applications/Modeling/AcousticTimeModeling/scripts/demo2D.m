
tic
mmx = linspace(0,6000,601);
mmz = linspace(0,3000,301);
mmy = 1;

nz = length(mmz);
nx = length(mmx);
ny = length(mmy);

vel = 1500 * ones(nz,nx,ny);
vel(100:end,:) = 1800; 
vel(150:end,:) = 2000; 
% initial model
Ps = opSmooth(nz,50);

vel_init = Ps * vel;
den = ones(size(vel));
% den(300:end,:) = 2;

vel 	 = vel(:);
vel_init = vel_init(:);
dv		 = 1./vel - 1./vel_init;
den		 = den(:);
% den = .23 * vel.^.25;


model.vel = vel;
model.den = den;


model_init.vel = vel_init;
model_init.den = den;


fwdpara.stype = 'seq_same';

fwdpara.slz   = [ mmz(2) ];
fwdpara.slx   = [ mmx(300) ];
fwdpara.sly   = [1];






fwdpara.saved = 1;
fwdpara.shot_id = 1;
dr  = 5;
rlx = mmx(10:5:end-10);
% rlx = fwdpara.slx:10:fwdpara.slx+1000;
rly = 1;
rlz = mmz(2);

[rlx,rlz,rly] = meshgrid(rlx,rlz,rly);

fwdpara.rtype = 'marine';
fwdpara.rlx = [rlx(:)];
fwdpara.rlz = [rlz(:)];
fwdpara.rly = [rly(:)];

fwdpara.fcent = 15;
fwdpara.dt    = .003;
t     = 4;
t0    = 1;
[src,fwdpara.taxis] = wvlet2(fwdpara.fcent,fwdpara.dt,t,t0);

% [src1,fwdpara.taxis] = wvlet(3,fwdpara.dt,t,t0);
 % src = zeros(size(src1));
  % src(1:2:end) = src1(1:2:end);
% src = src1;
fwdpara.abx = 40;
fwdpara.abz = 40;
fwdpara.aby = 40;

fwdpara.free = 1;

fwdpara.ddcompx = 1;
fwdpara.ddcompz = 1;
fwdpara.ddcompy = 1;

fwdpara.mmx = mmx;
fwdpara.mmz = mmz;
fwdpara.mmy = mmy;

fwdpara.chkp_space	= 50;
fwdpara.chkp_save	= 2; 


fwdpara.disp_snapshot = 1;
fwdpara.wave_equ = 1;
fwdpara.v_up_type = 'slowness';
fwdpara.space_order = 4;


fwdpara.fmode = 1;
[Ur,chk] = XiangLi_Time(model,src,fwdpara);
%keyboard
[Ur1,chk] = XiangLi_Time(model_init,src,fwdpara);
Ps = op_XiangLi_Time(model_init,src,fwdpara,chk);


dU = Ur - Ur1;


g = Ps' * dU;

%keyboard
g_lsrtm = lsqr_fwi(Ps,dU(:),10,1,20);

fwdpara.fmode = 2;
[Ur1] = XiangLi_Time(model,src,fwdpara);

toc
figure(1);subplot(1,2,1);imagesc(reshape(Ur,4001,117));subplot(1,2,2);imagesc(reshape(Ur1,4001,117))
figure(2);imagesc(reshape(Ur-Ur1,4001,117))

[Ur1,chk] = XiangLi_Time(model_init,src,fwdpara);

dUr1 = Ur - Ur1;
%

Ps = op_XiangLi_Time(model_init,src,fwdpara,chk);

dUr = Ps * dv;

g = Ps' * dUr1;



a = dUr1(:).' * dUr(:);
b = g' * dv;

disp(['dottest, difference is ',num2str(gather(a)-b),'; ratio is ', num2str(gather(a)/b)]);



% % [Ur1] = XiangLi_Time_old(vel,den,src,fwdpara);
% % norm(Ur - Ur1)
%
% % keyboard
% %
% % fwdpara.fmode = -1;
% % [src_back] = XiangLi_Time(vel,den,Ur,fwdpara);
% % [src_back1] = XiangLi_Time_old(vel,den,Ur1,fwdpara);
% % norm(src_back - src_back1)
%
% fwdpara.fmode = 3;
% tic
% [dUr] = XiangLi_Time(vel_init,den,src,fwdpara,0,dv);
% toc
% % tic
% %  [dUr1] = XiangLi_Time_old(vel_init,den,src,fwdpara,0,dv);
% % toc
% % % norm(dUr-dUr1)
% %
% fwdpara.fmode = -3;
% [g] = XiangLi_Time(vel_init,den,src,fwdpara,0,dUr);
% % % [g1] = XiangLi_Time_old(vel_init,den,src,fwdpara,0,dUr1);
% %
% % norm(g-g1)
%
%
%
% a = dUr1(:).' * dUr(:)
% b = g' * dv;
%
% disp(['dottest, difference is ',num2str(gather(a)-b),'; ratio is ', num2str(gather(a)/b)]);

% fwdpara.fmode = 2;
% [Ur2] = XiangLi_Time(vel,den,src,fwdpara,chk);
% keyboard;
%
% [Uri,chkp] = XiangLi_Time(vel_init,den,src,fwdpara);
%
%
%
% fwdpara.fmode = 3;
% [dUr] = XiangLi_Time(vel_init,den,src,fwdpara,0,dv);
%
% dUr1 = Ur - Uri;
%
% %
% Ps = op_XiangLi_Time(vel_init,den,src,fwdpara,chkp);
% %
% % % fwdpara.fmode = -3;
% % % [g,chkp] = XiangLi_Time(vel_init,den,src,fwdpara,0,dUr1);
% %
% g = Ps' * dUr1;
% %
%
% dottest for jacobian

%













