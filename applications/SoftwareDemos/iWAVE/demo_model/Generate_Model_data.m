%% test spot operator for iwave
function Generate_Model_data()

curdir   = pwd;
%mbinpath = [curdir '/../../../../tools/utilities/iWAVE'];
%addpath(mbinpath);
addpath(curdir);
% Important file names and their paths
fname_buoy    = 'cambuoy.rsf';
fname_bulkt   = 'cambulkt.rsf';
fname_bulkh   = 'cambulkh.rsf'; 
fname_dbulk   = 'camdbulk.rsf';
fname_dbuoy   = 'camdbuoy.rsf';
fname_ddata1  = 'camd1_data.su';
fname_data    = 'cam_data.su';

fname_bulklayert = 'layerbulkt.rsf';
fname_buoylayer  = 'layerbuoy.rsf';
fname_bulklayers = 'layerbulks.rsf';
fname_dbulklayer = 'layerdbulk.rsf';
fname_dbuoylayer = 'layerdbuoy.rsf';
fname_ddatalayer = 'layerd1_data.su';
fname_datalayer  = 'layer_data.su';

% Generate the cam model
n             = [201, 201];
d             = [5,5];
o             = [0,0];
v             = [2,2.5]; % low and high velocity

r             = 200;

cambuoy        = ones(n);
cambulkt       = (v(1)^2)*ones(n);
cambulkh       = cambulkt;

for i = 1:n(1)
	for j = 1:n(2)
		z = (i-1)*d(1);
		x = (j-1)*d(2);
		opoint = (n-1).*d/2;
		if norm(opoint - [z,x]) < r
			cambulkt(i,j) = v(2)^2;
		end
	end
end



rsf_write_all(fname_buoy,{'out=stdout'},cambuoy,d,o);
rsf_write_all(fname_bulkt,{'out=stdout'},cambulkt,d,o);
rsf_write_all(fname_bulkh,{'out=stdout'},cambulkh,d,o);
rsf_write_all(fname_dbulk,{'out=stdout'},cambulkh-cambulkt,d,o);
rsf_write_all(fname_dbuoy,{'out=stdout'},cambuoy-cambuoy,d,o);

m1    = [cambulkt(:);cambuoy(:)];
m2    = [cambulkh(:);cambuoy(:)];


model.t = 0:.004:2;
model.zsrc = 40;
model.xsrc = 500;
model.zrec = 40;
model.xrec = 0:10:1000;
model.f0   = 30;
model.t0   = .2 ;
model.n    = n;
model.o    = o;
model.d    = d;

options.fwdpara = 'fwdparams1.txt';
options.adjpara = 'adjparams.txt';
options.linpara = 'linparams.txt';
options.delete  = 1;
options.nlab     = 1;



D1     = Fd2DT(m1,model,options);
D2     = Fd2DT(m2,model,options);

label  = 'camd1';
iwave_wrtrac(model,label,D1-D2);

model.xsrc = 0:50:1000; 
model.zrec = 980;
options.fwdpara = 'fwdparamsm.txt';

clear D1

options.nlab     = 4;
D1    = Fd2DT(m1,model,options);
label = 'cam';
iwave_wrtrac(model,label,D1);




%% --------------------------------------------------------------------------------------
%% --------------------------------------------------------------------------------------
clear D1
clear D2


% Generate the layer model
layerbulkt             = ones(n)*v(1)^2;
layerbulkt(40:70,:)    = (v(1)*1.3)^2;
layerbulkt(71:120,:)   = (v(1)*1.6)^2;
layerbulkt(121:150,:)  = (v(1)*2)^2;
layerbulkt(151:end,:)  = (v(1)*2.5)^2;
opS                          = opSmooth(size(layerbulkt,1),16);
layerbulks             = opS * layerbulkt;
layerbuoy              = cambuoy;

rsf_write_all(fname_buoylayer,{'out=stdout'},cambuoy,d,o);
rsf_write_all(fname_bulklayert,{'out=stdout'},layerbulkt,d,o);
rsf_write_all(fname_bulklayers,{'out=stdout'},layerbulks,d,o);
rsf_write_all(fname_dbulklayer,{'out=stdout'},layerbulks-layerbulkt,d,o);
rsf_write_all(fname_dbuoylayer,{'out=stdout'},cambuoy-cambuoy,d,o);

m1    = [layerbulkt(:);layerbuoy(:)];
m2    = [layerbulks(:);layerbuoy(:)];

model.t = 0:.004:2;
model.zsrc = 40;
model.xsrc = 500;
model.zrec = 40;
model.xrec = 0:10:1000;
model.f0   = 30;
model.t0   = .2 ;
model.n    = n;
model.o    = o;
model.d    = d;

options.fwdpara = 'fwdparams1.txt';
options.adjpara = 'adjparams.txt';
options.linpara = 'linparams.txt';
options.delete  = 1;
options.nlab     = 1;




D1     = Fd2DT(m1,model,options);
D2     = Fd2DT(m2,model,options);

label  = 'layerd1';
iwave_wrtrac(model,label,D1-D2);

model.xsrc = 0:50:1000; 
options.fwdpara = 'fwdparamsm.txt';
options.nlab     = 4;
options.delete  = 0;


clear D1



D1    = Fd2DT(m1,model,options);
label = 'layer';
iwave_wrtrac(model,label,D1);



