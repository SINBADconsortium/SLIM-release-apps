classdef opMig_DeMig < opSpot
% opMig_DeMig(U,freq,nx,nz,nx_border,nz_border,xd,zd,vel)
%
%   OPMIGRATION Apply migration (demigration) matrix on data (image). 
%   this operator creates an N*NS*NF by N operator based on the 
%   forward modeling wavefields U(N,NS,NF) and FREQ(NF) that maps
%   back-propagated wavefields vec(V(N,NS,NF)) into an image
%   vector of size N. MODE = 1 acts as the forward (de-migration) 
%   operator. MODE = 2 (migration) brings the image vector back 
%   into the wavefield space.
%
%   Input 
%         U(N,NS,NF): forward modeling wavefieds
%                      with N  : number of gridpoints nx*nz
%                           NS : number of shots
%                           NF : number of frequencies
%
%           FREQ(NF)   frequencies used to generate forward/
%                      backpropagated wavefields
%           
%           vel         model, according to which we have U
%
%%   Author: Xiang LI
%
%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
       funHandle = []; % Multiplication function
    end % Properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opMig_DeMig(U,R,freq,nx,nz,nx_border,nz_border,xd,zd,vel)

		if      ndims(U) == 2
			if length(freq) ~= 1
				[n nf]       = size(U);
				ns           = 1;
				disp('Only one shot is used in the experment')
			else
				[n ns]       = size(U);
				nf           = 1;
				disp('Only one frequency is used in the experment')
			end
		elseif  ndims(U) == 3
		    [n ns nf]    = size(U);
		    disp([num2str(ns),' shots and ',num2str(nf),' frequencies are used in the experment'])
		end
		if nf ~= max(size(freq))
		  error('Number of frequencies in U and Freq does not match');
		end

		nx_tot        = nx + 2*nx_border;
		nz_tot        = nz + 2*nz_border;
		% picky         = zeros(n,1);
		% idex          = (nx_border*nz_tot+nz_border +4):nz_tot:(nx_border*nz_tot+nx*nz_tot);
		% picky(idex,1) = 1;
		% picky         = repmat(picky,1,ns);
		sigma         = 1./vel(:);
		sigma         = fixborder(sigma,nx,nz,nx_border,nz_border,1,2);
		sigma         = repmat(sigma,1,ns);

		fun = @(x,mode) opMig_DeMig_intrnl(U,R,freq,n,ns,nf,nx,nz,nx_border,nz_border,xd,zd,sigma,vel,x,mode);

          % Construct operator
          op = op@opSpot('Mig&DeMig',n*ns*nf,nx*nz);
          op.cflag     = 1;
          op.funHandle = fun;
       end % Constructor

    end % Methods
       
 
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = multiply(op,x,mode)
          y = op.funHandle(x,mode);
       end % Multiply          

    end % Methods
   
end % Classdef


%=======================================================================


function y = opMig_DeMig_intrnl(U,R,freq,n,ns,nf,nx,nz,nx_border,nz_border,xd,zd,sigma,vel,x,mode) 


if (mode == 1)  % demigration mode 
    % y              = zeros(n,ns,nf);
    x              = fixborder(x,nx,nz,nx_border,nz_border,1,2);
    x              = repmat(x,1,ns);
    parfor ifreq = 1:nf
        f          = freq(ifreq); 
        A          = opHelm2D9pt(vel,xd,zd,nx_border,nz_border,f);
        omg        = -2*(2*pi*f).^2;
        y1         = A\(omg.*U(:,:,ifreq).*sigma.*x);
        % y1         = y1.*picky; 
        y(:,:,ifreq) = R * y1;
    end
   y = y(:);
elseif (mode == 2)   % migration mode
    x              = reshape(x,n,ns,nf);
    y              = zeros(n,nf);
    parfor ifreq = 1:nf
        xidex      = x(:,:,ifreq).*picky;
        f          = freq(ifreq);
        A          = opHelm2D9pt(vel,xd,zd,nx_border,nz_border,f);       
        omg        = -2*(2*pi*f).^2;
        y1         = omg.*conj(U(:,:,ifreq)).*sigma.*(A'\xidex);
        y(:,ifreq) = sum(y1,2);
    end
    y              = sum(y,2);
    y              = fixborder(y,nx,nz,nx_border,nz_border,1,1);
end
end