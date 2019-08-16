

mmx = linspace(0,3000,301);
mmz = linspace(0,2000,201);
mmy = 1;

nz = length(mmz);
nx = length(mmx);
ny = length(mmy);

model.vel = 2000 * ones(nz,nx,ny);
model.vel(1:100,:,:) = 1500;
model.vel(150:end,:,:) = 2200;


load ../data/Water_Mar
model.vel= vel_true( 1:201,300:600);
mmx = linspace(0,3000,301);
mmz = linspace(0,2000,201);

% initial model
Ps = opSmooth(nz,100);

model_init.vel = Ps * model.vel;
den = 1;
model.vel 	 = model.vel(:);vel = model.vel;
model_init.vel = model_init.vel(:);

model.den		 = den(:);
% den = .23 * vel.^.25;
model_init.den		 = den(:);

fwdpara.stype = 'seq_same';
% fwdpara.slz   = [mmz(15)];
% fwdpara.slx   = [mmx(50)];
% fwdpara.sly   = [mmy(56)];
fwdpara.slz   = [mmz(20)]+5;
fwdpara.slx   = [mmx(100)+5];
fwdpara.sly   = [1];



fwdpara.saved = 1;
fwdpara.shot_id = 1;
dr  = 5;
rlx = mmx(20:2:end-20);
rly = 1;
rlz = mmz(15);

[rlx,rlz,rly] = meshgrid(rlx,rlz,rly);

fwdpara.rtype = 'full';
fwdpara.rlx = rlx(:);
fwdpara.rlz = rlz(:);
fwdpara.rly = rly(:);

fwdpara.fcent = 14;
fwdpara.dt    = .0013;
t     = 2;
t0    = .1;
[src,fwdpara.taxis] = wvlet2(fwdpara.fcent,fwdpara.dt,t,t0);


fwdpara.abx = 10;
fwdpara.abz = 10;
fwdpara.aby = 10;

fwdpara.free = 1;

fwdpara.ddcompx = 1;
fwdpara.ddcompz = 1;
fwdpara.ddcompy = 1;

fwdpara.mmx = mmx;
fwdpara.mmz = mmz;
fwdpara.mmy = mmy;

fwdpara.chkp_space	= 80;
fwdpara.chkp_save	= 2; 

fwdpara.disp_snapshot = 0;
fwdpara.wave_equ = 1;
fwdpara.v_up_type = 'slowness';
fwdpara.space_order = 4;


fwdpara.fmode = 1;




dv = model.vel - model_init.vel;


h = 10.^[0:-1:-6];
err = zeros(length(h),2);


% Jacobian test

for mm = 1:length(h)

	model.vel  = model_init.vel + h(mm) * dv;

	dm = 1./model.vel - 1./model_init.vel;

	[Ur1] 		= XiangLi_Time(model,src,fwdpara);
	[Ur,chk] 	= XiangLi_Time(model_init,src,fwdpara);
	J   		= op_XiangLi_Time(model_init,src,fwdpara,chk);

	

	dU1 = J * dm;


	err(mm,1) = norm(Ur1 - Ur);

	err(mm,2) = norm(Ur1 - Ur  - dU1);




	% dot test
	dU 	= Ur - Ur1;
	g 	= J' * dU;



	% dottest for jacobian
	a = dU(:).' * dU1(:)
	b = g' * dm;

	disp(['dottest, difference is ',num2str(gather(a)-b),'; ratio is ', num2str(gather(a)/b)]);


end
figure(1);set(gca,'fontsize',18);
loglog(h, [err(:,1) err(:,2) h(:) h(:).^2],'*-','lineWidth',2);
title('Jacobian test')
legend('U(m+h\delta dm)-U(m)','U(m+h\delta dm)-U(m)-hJ(m)\delta m','h','h^2')

	
% gradient test

model.vel = vel;
[Urt]      = XiangLi_Time(model,src,fwdpara);
[Ur,chk] = XiangLi_Time(model_init,src,fwdpara);
J    = op_XiangLi_Time(model_init,src,fwdpara,chk);

dU = Ur - Urt;

% dU1 = J * dm;
g = J' * dU;


dm = 1./model.vel - 1./model_init.vel;
h1 = -10.^[0:-1:-6];

err_fwiobj = zeros(length(h1),2);
for mm = 1:length(h1)
	
	model.vel   = 1./(1./model_init.vel + h1(mm)*dm);
	
%	dm    = 1./vel - 1./vel_init;
	
	[Ur1] = XiangLi_Time(model,src,fwdpara);

	
	err_fwiobj(mm,1) = 0.5*norm(Ur1 - Urt).^2 - 0.5*norm(Ur-Urt).^2;
	
	err_fwiobj(mm,2) = 0.5*norm(Ur1 - Urt).^2 - .5*norm(Ur-Urt).^2 - h1(mm) *g.' * dm;
	
	
	
	
	% g = J' * dU;




	% dottest for jacobian
	% a = dU(:).' * dU1(:)
% 	b = g' * dm;
%
% 	disp(['dottest, difference is ',num2str(gather(a)-b),'; ratio is ', num2str(gather(a)/b)]);
%
	
end
figure(2);set(gca,'fontsize',18);
loglog(-h1, [err_fwiobj abs(h1') (h1.^2)'],'*-','lineWidth',2);
title('gradient test')
legend('\Phi(m+h\delta dm)-\Phi(m)','\Phi(m+h\delta dm)-\Phi(m)- h\nabla\Phi(m)\delta m','h','h^2')














