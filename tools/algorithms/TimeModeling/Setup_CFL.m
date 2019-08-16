function [v,model,den,ani]=Setup_CFL(v,model,den,ani)
    %% CFL setup function 
    % This function will compute the grid size and time step according to
    % the source peak frequency, the minimum and maximum velocities and the
    % discretization order in space
    %
    % INPUT :
    %       v : square slowness(in s^2/km^2)
    %       model :  structure containing model parameters
    %       den : density if needed
    %       ani : anisotropy (sructure with epsilon,delta,theta)
    %
    % OUTPUT : output are the input parameters on thee new grid with the
    % new grid size and time step inside the model structure
    %       v : square slowness(in s^2/km^2)
    %       model :  structure containing model parameters
    %       den : density if needed
    %       ani : anisotropy (sructure with epsilon,delta,theta)
	% Author : Mathias Louboutin, Philipp Witte
	%         Seismic Laboratory for Imaging and Modeling
	%         Department of Earth, Ocean, and Atmosperic Sciences
	%         The University of British Columbia
    % Date : July, 2015
	f=model.f0;
	if model.space_order~=2 && model.space_order~=4 && model.space_order~=6
		error(['Only 2nd, 4th and 6th order are suported']);
	end
	
	if model.freesurface
	    nbz1=model.space_order/2;
	else
	    nbz1=40;
	end


	minm=min(v(:));
	maxm=max(v(:));
	%% minimum and maximum velocities
	minv=1e3*min(sqrt(1./v(:)));
	maxv=1e3*max(sqrt(1./v(:)));
	%% Discretization constants
	if model.space_order == 2
		a		=sqrt(1/2);
		if model.n(3)>1
			a	= sqrt(1/3); 
		end
	elseif model.space_order == 4 
		a		=sqrt(4/11);
		if model.n(3)>1
			 a	=.5; 
		end
	elseif model.space_order == 6 
		a		=0.5752;
		if model.n(3)>1
			 a	=0.4697; 
		end
	end

	%% Grid size
    if nargin==4
        % anisotropic modeling with spectral derivative in space
        G = model.gppwl;
    else
        G = 10;
    end
    d  = minv/(G*1e3*f);
    d=floor(d);
    
	%% Time stepping rate
    if nargin==4
        dt = 2* a * d /maxv/pi;
    else
        dt = a * d /maxv;
    end
	dt=1e3*(dt-dt/10);
	d=[d d d];
	disp(['CFL conditions gives dt = ' num2str(dt) 'ms and d = ' num2str(d) ' m ']);

	%% Interpolate to new grid

	if model.n(3)>1
		[dz,dx,dy]=odn2grid(model.o,model.d,model.n);

		nxnew=floor((dx(end)-model.o(2))/d(2))+1;
		nznew=floor((dz(end)-model.o(1))/d(1))+1;
		nynew=floor((dy(end)-model.o(3))/d(3))+1;

		[dzn,dxn,dyn]=odn2grid(model.o,d,[nznew nxnew nynew]);

		P=opKron(opLInterp1D(dy,dyn),opLInterp1D(dx,dxn),opLInterp1D(dz,dzn));

	else
		[dz,dx,dy]=odn2grid(model.o,model.d,model.n);

		nxnew=floor((dx(end)-model.o(2))/d(2))+1;
		nznew=floor((dz(end)-model.o(1))/d(1))+1;
		nynew=1;
        
        % for anisotropy with spectral derivative, force even # of grid
        % points
        if nargin==4
			if model.ddcompx>1
	            while mod((nxnew+80)/model.ddcompx,2)~=0
	                nxnew = nxnew+1;
	            end
			else
	            if mod(nxnew,2)~=0
	                nxnew = nxnew+1;
	            end
			end
			if model.ddcompz>1
	            while mod((nznew+40+nbz1)/model.ddcompz,2)~=0
	                nznew = nznew+1;
	            end
			else
	            if mod(nznew,2)~=0
	                nznew = nznew+1;
	            end
			end
        end
        
		[dzn,dxn,dyn]=odn2grid(model.o,d,[nznew nxnew nynew]);

		P=opKron(opLInterp1D(dx,dxn),opLInterp1D(dz,dzn));

	end
	v=P*v(:);

	v(v<minm)=minm;
	v(v>maxm)=maxm;

	if nargin==3
		mind=min(den(:));
		maxd=max(den(:));
		den=P*den(:);

		den(den<mind)=mind;
		den(den>maxd)=maxd;
        
    elseif nargin==4
        mind=min(ani.epsilon(:));
		maxd=max(ani.epsilon(:));
		ani.epsilon=P*ani.epsilon(:);
        
		ani.epsilon(ani.epsilon<mind)=mind;
		ani.epsilon(ani.epsilon>maxd)=maxd;
        
        mind=min(ani.delta(:));
		maxd=max(ani.delta(:));
		ani.delta=P*ani.delta(:);

		ani.delta(ani.delta<mind)=mind;
		ani.delta(ani.delta>maxd)=maxd;
        
        mind=min(ani.theta(:));
		maxd=max(ani.theta(:));
		ani.theta=P*ani.theta(:);

		ani.theta(ani.theta<mind)=mind;
		ani.theta(ani.theta>maxd)=maxd;
        
        if ~isempty(den)
            mind=min(den(:));
            maxd=max(den(:));
            den=P*den(:);

            den(den<mind)=mind;
            den(den>maxd)=maxd;
        end
	end

	model.d=d;
	model.dt=dt;
	model.n=[nznew nxnew nynew];

	disp('Velocity interpolated on new grid');
