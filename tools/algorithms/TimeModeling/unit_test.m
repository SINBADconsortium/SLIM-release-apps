clear


%% test op_Src_Rcv_intep
x = 1:100;
y = 100:200;
z = 50:200;


x0 = 30:20:80;
y0 = 120:21:190;
z0 = 60:12:180;


Pr1 = opKron(opLInterp1D(y,y0),opLInterp1D(x,x0),opLInterp1D(z,z0));

[x1,z1,y1] = meshgrid(x0,z0,y0);


Pr = op_Src_Rcv_intep(z,x,y,z1(:),x1(:),y1(:));

dottest(Pr,1,1);

a = randn(size(Pr,2),1);
norm(Pr * a - Pr1*a);

a = randn(size(Pr,1),1);
norm(Pr' * a - Pr1'*a);


%% adjoint test for 2D forward modeling kernel 

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


fwdpara.disp_snapshot = 0;
fwdpara.wave_equ = 1;
fwdpara.v_up_type = 'slowness';
fwdpara.space_order = 4;


fwdpara.fmode = 1;
[Ur,chk]  = XiangLi_Time(model,src,fwdpara);
[Ur1,chk] = XiangLi_Time(model_init,src,fwdpara);
Ps 		  = op_XiangLi_Time(model_init,src,fwdpara,chk);


dUr1 = Ur - Ur1;

g = Ps' * dUr1;
dUr = Ps * dv;

a = dUr1(:).' * dUr(:);
b = g' * dv;

disp(['dottest, difference is ',num2str(gather(a)-b),'; ratio is ', num2str(gather(a)/b)]);



%% adjoint test for 3D forward modeling kernel

mmx = linspace(0,1000,101);
mmy = linspace(0,1000,101);
mmz = linspace(0,1000,101);

nz = length(mmz);
nx = length(mmx);
ny = length(mmy);

vel = 1500 * ones(nz,nx,ny);
vel(round(nz./2):end,:,:) = 2000;

vel = vel(:);

den = .23 * vel.^.25;

model.vel = vel;
model.den = den;

Ps = opKron(opDirac(101),opDirac(101),opSmooth(101,20));
model_init.vel = Ps * vel;
model_init.den = den;


fwdpara.stype = 'seq_same';
% fwdpara.slz   = [mmz(15)];
% fwdpara.slx   = [mmx(50)];
% fwdpara.sly   = [mmy(56)];
fwdpara.slz   = [mmz(20)];
fwdpara.slx   = [mmx(round(nx./2))];
fwdpara.sly   = [mmy(round(ny./2))];
fwdpara.saved = 1;
fwdpara.shot_id = 1;

dr  = 50;
rlx = mmx(15):dr:mmx(end-15);
rly = mmy(15):dr:mmy(end-15);
rlz = mmz(10);

[rlx,rlz,rly] = meshgrid(rlx,rlz,rly);

fwdpara.rtype = 'full';
fwdpara.rlx = rlx(:);
fwdpara.rlz = rlz(:);
fwdpara.rly = rly(:);

fwdpara.fcent = 14;
fwdpara.dt    = .002;
t     		  = 1;
t0    		  = .2;
[src,fwdpara.taxis] = wvlet2(fwdpara.fcent,fwdpara.dt,t,t0);


fwdpara.abx = 20;
fwdpara.abz = 20;
fwdpara.aby = 20;

fwdpara.free = 0;

fwdpara.ddcompx = 1;
fwdpara.ddcompz = 1;
fwdpara.ddcompy = 1;

fwdpara.mmx = mmx;
fwdpara.mmz = mmz;
fwdpara.mmy = mmy;

fwdpara.v_up_type = 'slowness';
fwdpara.chkp_space	= 20;
fwdpara.chkp_save	= 3; 

fwdpara.disp_snapshot = 0;
fwdpara.wave_equ	= 1;

fwdpara.fmode = 1;

[Ur,chkp] = XiangLi_Time(model,src,fwdpara);

% fwdpara.fmode = 2;
% [Ur1] = XiangLi_Time_old(vel,den,src,fwdpara,chkp);
[Uri,chkp] = XiangLi_Time(model_init,src,fwdpara);

dUr1 = Ur - Uri;


Ps = op_XiangLi_Time(model_init,src,fwdpara,chkp);

dv = model.vel - model_init.vel;

dUr = Ps * dv;

g = Ps' * dUr1;

% dottest for jacobian
a = dUr1(:).' * dUr(:);
b = g' * dv;
% dm = lsqr_fwi(Ps,dUr1,5,1);
disp(['dottest, difference is ',num2str(gather(a)-b),'; ratio is ', num2str(gather(a)/b)]);