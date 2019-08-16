function test_suite = test_TimeStepping_parallel
	% unit tests for 2D Time-domain modelling
	% THis is a Serial code but calling spmd inside the function therfor it will open a parpool
	% Mathias Louboutin,2015
	% mloubout@eos.ubc.ca
	%
	test_suite=buildFunctionHandleTestSuite(localfunctions);
end

function model = setup
	v=1.5*ones(201,201);
	v(100:end,:)=2.5;
	model.m=1./v.^2;
	model.m0=1/(1.5)^2*ones(size(model.m));
	%% Model
	model.o=[0 0 0]; %Origins of the axes [m]
	model.n=[201 201 1]; %Number of grid points  for each dimension (excluding boundaries) [z x y]
	model.ddcompx=1; % Domain decomposition x direction
	model.ddcompy=1; % Domain decomposition y direction
	model.ddcompz=1; % Domain decomposition z direction
	model.f0=.008; % Source peak frequency in Hz
	model.d=[15 15 15]; % Grid size in each direction (same in every direction unless really necessary)
	model.T=2000; %Acquisition duration [s]
	model.freesurface=0; % Freesurface ( 0 : no freesurface, 1 : freesurface)
	model.space_order=2; % Space discretization order (2 or 4 only for now)
	model.type='full';
	model.gppwl = 6; 
	N=prod(model.n);
	ani.delta = ones(N,1)*0.0;
	ani.epsilon = ones(N,1)*0.0;
	ani.theta = ones(N,1)*0;
	disp('setup');
	model.ani=ani;

	model.xsrc=(model.n(2)-1)*model.d(2)/2;
	model.ysrc=(model.n(3)-1)*model.d(3)/2;
	model.zsrc=10*model.d(1);

	model.xrec=0:model.d(2):(model.n(2)-1)*model.d(2);
	model.yrec=(model.n(3)-1)*model.d(3)/2;
	model.zrec=10*model.d(1);
end


%% Adjoint test for forward modelling

function testAdjointF(model)
	disp('F^T');
	[m,model1,m0,ani]=Setup_CFL(model.m,model,model.m0,model.ani);
	model1.NyqT=0:model1.dt:model1.T;
	nt=length(0:model1.dt:model1.T);
	%% Adjoint test F
	source=sp_RickerWavelet(model1.f0,1/model1.f0,model1.dt,model1.T);
	
	data=Gen_data(m,model1,source,[],ani);
	data2=randn(size(data));
	source2=Bak_propag(m,model1,data2,[],ani);

	dd=data2'*data;
	ss=source2*source;
	assert(abs(dd/ss-1)<1e-3,'Forward modeling does not pass the adjoint test');
end


% Adjoint test Jacobian
function testAdjointJ(model)
	disp('J^T');
	[m,model1,m0]=Setup_CFL(model.m,model,model.m0);
	model1.save='RAM';
	dm=m-m0;
	model1.NyqT=0:model1.dt:model1.T;
	nt=length(0:model1.dt:model1.T);
	%% Source
	source=sp_RickerWavelet(model1.f0,1/model1.f0,model1.dt,model1.T);

	%% Adjoint test J
	data2=Born(m,model1,source,dm,1);
	Im=Born(m,model1,source,data2,-1);

	dd=data2'*data2;
	ss=dm'*Im;
	assert(abs(dd/ss-1)<1e-1,'Jacobian does not pass the adjoint test');
end

% Gradient test

function testGrad(model)
	disp('G');
	[m,model1,m0]=Setup_CFL(model.m,model,model.m0);
	model1.NyqT=0:model1.dt:model1.T;
	model1.save='RAM';
	nt=length(0:model1.dt:model1.T);
	q=10*   sp_RickerWavelet(model1.f0,1/model1.f0,model1.dt,model1.T);
	%% Gradient test
	h   = 10.^(-1:-1:-6);
	dm=(m-m0);

	% Forward
	D1=Gen_data(m0,model1,q);
	du=Born(m0,model1,q,dm,1);
	fprintf('Gradient value - |J(mo)*dm|                 = %12.8f \n',norm(du,1));
	for i = 1:length(h)
	    mloc=m0+h(i)*dm;
	    d=Gen_data(mloc,model1,q);
	    error1(i)      = sum(abs(d - D1));
	    fprintf('Perturbation - |h|                                      = %12.8f \n', h(i));
	    fprintf('First order error - |f(m0+h*dm) - f(m0) |               = %12.8f \n', error1(i));
	    error2(i)  = sum(abs(d - D1 - h(i)*du));
	    fprintf('Second order correction - |h*J(m0)*dm|                  = %12.8f \n',sum(abs(h(i)*du)) );
	    fprintf('Second order error - |f(m0+h*dm) - f(m0) - h*J(m0)*dm| = %12.8f \n', error2(i));
	    fprintf('=========================================================================\n\n');
	end
	
    error1=log10(error1);
	P1=polyfit(log10(h),error1,1);
	error2=log10(error2);
	P2=polyfit(log10(h),error2,1);
    assert( abs( P1(1) - 1) < 1e-1 ,'First order error is not O(h)');
    assert( abs( P2(1) - 2) < 1e-1 ,'Second order error is not O(h^2)');
end

