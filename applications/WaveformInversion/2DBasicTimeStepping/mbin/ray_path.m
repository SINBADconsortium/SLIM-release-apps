

x1 = 20000;
x2 = 32000;

% load model
load ./data/mode_for_fwi
idx  = find(mmx>x1 & mmx<x2);
mmx      = mmx(idx);
vel_init = vel_init(:,idx);


x = mmx;
z = mmz;




nx=length(x);
nz=length(z);
dg=x(2)-x(1);

v = vel_init;

 rayvelmod(vel_init,dg);

%estimate tmax,dt,tstep
% vlow=min(min(vel_init));
% tmax=max(z)/vlow;dt=.004;tstep=0:dt:tmax;

%specify a fan of rays
angles=[0.5:2:80]*pi/180;
x0=mmx(40)-mmx(1);z0=0;
indx=near(x,x0);indz=near(z,z0);
v0=v(indz,indx);

%trace the rays
tstep = 0:0.004:9;

imagesc(x-x(1),z,((v)));
for k=1:length(angles)
	r0=[x0 z0 sin(angles(k))/v0 cos(angles(k))/v0];
	[t,r]=shootrayvxz(tstep,r0);
	line(r(:,1),r(:,2),ones(size(t)),'color','w','lineWidth',1);
end

