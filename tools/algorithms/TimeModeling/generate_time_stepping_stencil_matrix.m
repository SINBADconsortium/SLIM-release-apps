function [A1,A1_inv,A2,A3,A3_inv,dA1,dA2,dA3] = generate_time_stepping_stencil_matrix(v,model,den,thomsen)
%% This function can give you a stencil matrix for time stepping
% A3 Ut-1 + A2 Ut + A1 Ut+1 = Qt
% code up difference mode.
% mode = 1, 4 th order finite difference with acoustic wave equation of constant density
% INPUTS
%
%   'model' : stucture containing the model parameters
%       'v'  :  square slowness [s^2/m^2]
%     'den' : density [g/cm^3]
%     'ani' : anisotropy parameters (epsilon,delta,theta in a
% Author : Mathias Louboutin, Philipp Witte
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
% Date : July, 2015

% Size of the model
nx = model.n(2);
ny = model.n(3);
nz = model.n(1);
v=reshape(v,nz,nx,ny);
% define PML parameter
% create pml boundary tensor (based on paper :
%Efficient PML for the wave equation
%Marcus J. Grote1 , and Imbo Sim2?)
dampf       = log(1/0.001)/(model.nb(2)*model.d(2)); % damping factor
yeta_x      = [fliplr(linspace(0,1,model.nb(2))),zeros(1,nx-2*model.nb(2)),linspace(0,1,model.nb(2))];
yeta_x=dampf.*(abs(yeta_x)-sin(2*pi*abs(yeta_x))/(2*pi));
if model.freesurface == 0
    yeta_z	= [fliplr(linspace(0,1,model.nb(1))),zeros(1,nz-2*model.nb(1)), linspace(0,1,model.nb(1))];
else
    yeta_z	= [zeros(1,nz-model.nb(1)),linspace(0,1,model.nb(1))];
end
yeta_z=dampf.*(abs(yeta_z)-sin(2*pi*abs(yeta_z))/(2*pi));

yeta_y      = 0;
if ny ~= 1
    yeta_y	= [fliplr(linspace(0,1,model.nb(3))),zeros(1,ny-2*model.nb(3)),linspace(0,1,model.nb(3))];
    yeta_y=dampf.*(abs(yeta_y)-sin(2*pi*abs(yeta_y))/(2*pi));
end
yeta_x = yeta_x(model.x_id_loc);
yeta_y = yeta_y(model.y_id_loc);
yeta_z = yeta_z(model.z_id_loc);
yeta=zeros(nz,nx,ny);
yeta2=zeros(nz,nx,ny);
for i=1:length(yeta_x)
    for j=1:length(yeta_y)
        for k=1:length(yeta_z)
            yeta(k,i,j)		=  1/sqrt(v(k,i,j))*(yeta_x(i) + yeta_z(k) + yeta_y(j));
            yeta2(k,i,j)		= 1/v(k,i,j)*(yeta_x(i) * yeta_z(k)+yeta_y(j) * yeta_z(k)+yeta_x(i) * yeta_y(j));
        end
    end
end

% size of the local model and wavefield
nx_loc_m 	= length(model.u_idx);
ny_loc_m 	= length(model.u_idy);
nz_loc_m 	= length(model.u_idz);

nx_loc_u 	= length(model.mmx_loc);
ny_loc_u 	= length(model.mmy_loc);
nz_loc_u 	= length(model.mmz_loc);
n_loc_u  	= nx_loc_u * ny_loc_u * nz_loc_u;

[zz,xx,yy] = ndgrid(model.u_idz,model.u_idx,model.u_idy);
U_index    = sub2ind([nz_loc_u,nx_loc_u,ny_loc_u],zz(:),xx(:),yy(:));
clear xx yy zz

if isempty(den) && isempty(thomsen)
    % if labindex==1
        % disp('Acoustic isotropic');
    % end
    %% ===============================================
    %    acoustic wave equation with only velocity
    %% ===============================================
    %              1   dU
    %     L U   + --- ---- = q
    %             V^2  dt
    %	discretization
    %
    %	A1 U_t-1 + A2 U_t + A3 U_t+1 = Q_t+1
    %
    % Including PML boundary
    %	A1 = m^2(1+yeta*dt)/dt^2;
    %	A2 = m^2(-2/dt^2  + yeta^2) - L;
    %	A3 = m^2(1-yeta*dt)/dt^2;
    %
    % -----------------------------------------------
    
    
    A1_inv		= (model.dt^2./v(:))./(1+yeta(:).*model.dt);
    A1_inv		= spdiags(A1_inv,0,n_loc_u,n_loc_u);

    A3 			= zeros(n_loc_u,1);
    A3(U_index)	= (1-vec(yeta(U_index)).*model.dt).*vec(v(U_index))./(model.dt^2);
    A3			= spdiags(A3,0,n_loc_u,n_loc_u);

    A1 			= zeros(n_loc_u,1);
    A1(U_index)	= (1+vec(yeta(U_index)).*model.dt).*(vec(v(U_index)))./model.dt^2;
    A1			= spdiags(A1,0,n_loc_u,n_loc_u);

    A3_inv		=  (model.dt^2./v(:))./(1-yeta(:).*model.dt);
    A3_inv		= spdiags(A3_inv,0,n_loc_u,n_loc_u);
    
    
    
    % ------------------- finite difference stencil coefficients ----------------------------
    switch model.space_order
        case 2
            stencil_ceo = [0 1 -2 1 0]; % 2nd order
            G_np_extra  = 2;
        case 4
            stencil_ceo = [-1, 16,-30,16,-1]./12; % 4 th order
            % stencil_ceo = [-0.1,1.4,-2.6,1.4,-0.1]; % Lluis's coefficient, contact him for more info
            %  "Guasch, Lluis" <l.guasch08@imperial.ac.uk>
            G_np_extra  = 2;
        case 6
            stencil_ceo = [			1/90	-3/20	3/2	-49/18	3/2	-3/20	1/90]; % 6th order
            G_np_extra  = 3;
        case 8
            stencil_ceo = [-1/560	8/315	-1/5	8/5	-205/72	8/5	-1/5	8/315	-1/560];
            G_np_extra  = 4;
    end
    %
    % ---------------------------------------------------------------------------------------
    dx 									= model.d(2);
    if ny~=1
        dy							    = model.d(3);
    end
    dz 									= model.d(1);
    Lx									= zeros(nx_loc_u,length(stencil_ceo));
    if ny~=1
        Ly							    = zeros(ny_loc_u,length(stencil_ceo));
    end
    Lz									= zeros(nz_loc_u,length(stencil_ceo));
    Lx(model.u_idx,:)				    = ones(nx_loc_m,1)*stencil_ceo./(dx^2);
    if ny~=1
        Ly(model.u_idy,:)	            = ones(ny_loc_m,1)*stencil_ceo./(dy^2);
    end
    Lz(model.u_idz,:)				    = ones(nz_loc_m,1)*stencil_ceo./(dz^2);
    
    diag_U_ind_x						= zeros(nx_loc_u,1);
    diag_U_ind_x(model.u_idx)		    = 1;
    diag_U_ind_y						= zeros(ny_loc_u,1);
    diag_U_ind_y(model.u_idy)		    = 1;
    diag_U_ind_z						= zeros(nz_loc_u,1);
    diag_U_ind_z(model.u_idz)		    = 1;
    
    Lx			    = spdiags(Lx,-G_np_extra:G_np_extra,nx_loc_u,nx_loc_u)';
    Lx 		 		= kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),kron(Lx,spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
    
    Lz			= spdiags(Lz,-G_np_extra:G_np_extra,nz_loc_u,nz_loc_u)';
    Lz 			= kron(kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u)),Lz);
    
    if ny ~= 1
        Ly		= spdiags(Ly,-G_np_extra:G_np_extra,ny_loc_u,ny_loc_u)';
        Ly			= kron(Ly,kron(spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u),spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
        L			= Lx + Lz + Ly;clear Lx Ly Lz
    elseif ny == 1
        L			= Lx + Lz;clear Lx Ly Lz
    end
    
    M				=  zeros(n_loc_u,1);
    M(U_index)		=  (-2/model.dt^2 + vec(yeta2(U_index))).*vec(v(U_index));
    M				=  spdiags(M,0,n_loc_u,n_loc_u);
    A2 			= M - L; clear L M
    
elseif ~isempty(den) && isempty(thomsen)
    % if labindex==1
        % disp('Elastic isotropic');
    % end
    %% ===============================================
    %    elastic wave equation
    %% ===============================================
    %        1           1       dU
    %     G --- G U   + ------- ---- = q
    %       rho         rho V^2  dt
    %	discretization
    %
    %	A1 U_t-1 + A2 U_t + A3 U_t+1 = Q_t+1
    %
    % Including PML boundary
    %	A1 = -(1+yeta*dt)/(v^2*dt^2);
    %	A2 = G diag(1/rho) G - (1/v^2)(-2/dt^2  + yeta^2);
    %	A3 = -(1-yeta*dt)/(v^2*dt^2);
    %
    % -----------------------------------------------
    % ------------------- finite difference stencil coefficients ----------------------------
    %  To include density, the real order is = fwdpara.space_order./2
    
    stencil_ceo = [-1 1]; % 2nd order
    G_np_extra    = 1;
    
    % create gradient operator 'G'
    dx 									= model.d(2);
    if ny~=1
        dy							    = model.d(3);
    end
    dz 									= model.d(1);
    Lx									= zeros(nx_loc_u,length(stencil_ceo));
    if ny~=1
        Ly							    = zeros(ny_loc_u,length(stencil_ceo));
    end
    Lz									= zeros(nz_loc_u,length(stencil_ceo));
    Lx(model.u_idx,:)				    = ones(nx_loc_m,1)*stencil_ceo./(dx);
    if ny~=1
        Ly(model.u_idy,:)	            = ones(ny_loc_m,1)*stencil_ceo./(dy);
    end
    Lz(model.u_idz,:)				    = ones(nz_loc_m,1)*stencil_ceo./(dz);
    
    diag_U_ind_x						= zeros(nx_loc_u,1);
    diag_U_ind_x(model.u_idx)		    = 1;
    diag_U_ind_y						= zeros(ny_loc_u,1);
    diag_U_ind_y(model.u_idy)		    = 1;
    diag_U_ind_z						= zeros(nz_loc_u,1);
    diag_U_ind_z(model.u_idz)		    = 1;
    
    G				= spdiags(vec(1./den(:)),0,n_loc_u,n_loc_u);
    % lx
    Lx1				= fliplr(Lx);
    Lx1				= spdiags(Lx1, -G_np_extra:0,nx_loc_u,nx_loc_u)';
    Lx1 		 	= kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),kron(Lx1,spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
    
    Gx				= .5 * ones(nx_loc_u,length(stencil_ceo));Gx(end,1) = 1;
    Gx				= spdiags(Gx, 0:G_np_extra,nx_loc_u,nx_loc_u);
    Gx 		 		= kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),kron(Gx,spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
    Gx				= spdiags(Gx*vec(1./den(:)),0,n_loc_u,n_loc_u);
    
    Lx1				= Gx * Lx1;
    
    Lx2				= fliplr(Lx);
    Lx2				= spdiags(Lx2, 0:G_np_extra,nx_loc_u,nx_loc_u)';
    Lx2 		 	= kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),kron(Lx2,spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
    
    Gx				= .5 * ones(nx_loc_u,length(stencil_ceo));Gx(1,2) = 1;
    Gx				= spdiags(Gx, -G_np_extra:0,nx_loc_u,nx_loc_u);
    Gx 		 		= kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),kron(Gx,spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
    Gx				= spdiags(Gx*vec(1./den(:)),0,n_loc_u,n_loc_u);
    Lx2				= Gx * Lx2;
    
    Lx				= (-Lx2  + Lx1)./dx;	clear Lx1 Lx2 Gx
    
    %lz
    
    Lz1				= fliplr(Lz);
    Lz1				= spdiags(Lz1, -G_np_extra:0,nz_loc_u,nz_loc_u)';
    Lz1 		 	= kron(kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u)),Lz1);
    
    Gz				= .5 * ones(nz_loc_u,length(stencil_ceo));Gz(end,1) = 1;
    Gz				= spdiags(Gz, 0:G_np_extra,nz_loc_u,nz_loc_u);
    Gz 		 		= kron(kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u)),Gz);
    Gz				= spdiags(Gz*vec(1./den(:)),0,n_loc_u,n_loc_u);
    
    Lz1				= Gz * Lz1;
    
    Lz2				= fliplr(Lz);
    Lz2				= spdiags(Lz2, 0:G_np_extra,nz_loc_u,nz_loc_u)';
    Lz2 		 	= kron(kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u)),Lz2);
    
    Gz				= .5 * ones(nz_loc_u,length(stencil_ceo));Gz(1,2) = 1;
    Gz				= spdiags(Gz, -G_np_extra:0,nz_loc_u,nz_loc_u);
    Gz 		 		= kron(kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u)),Gz);
    Gz				= spdiags(Gz*vec(1./den(:)),0,n_loc_u,n_loc_u);
    Lz2				= Gz * Lz2;
    
    Lz				= (- Lz2  + Lz1)./dz;	clear Lz1 Lz2 Gz
    
    
    if ny ~= 1
        
        Ly1				= fliplr(Ly);
        Ly1				= spdiags(Ly1, -G_np_extra:0,ny_loc_u,ny_loc_u)';
        Ly1 		 	= kron(Ly1,kron(spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u),spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
        
        Gy				= .5 * ones(ny_loc_u,length(stencil_ceo));Gy(end,1) = 1;
        Gy				= spdiags(Gy, 0:G_np_extra,ny_loc_u,ny_loc_u);
        Gy 		 		= kron(Gy,kron(spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u),spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
        Gy				= spdiags(Gy*vec(1./den(:)),0,n_loc_u,n_loc_u);
        
        Ly1				= Gy * Ly1;
        
        Ly2				= fliplr(Ly);
        Ly2				= spdiags(Ly2, 0:G_np_extra,ny_loc_u,ny_loc_u)';
        Ly2 		 	= kron(Ly2,kron(spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u),spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
        
        Gy				= .5 * ones(ny_loc_u,length(stencil_ceo));Gy(1,2) = 1;
        Gy				= spdiags(Gy, -G_np_extra:0,ny_loc_u,ny_loc_u);
        Gy 		 		= kron(Gy,kron(spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u),spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
        Gy				= spdiags(Gy*vec(1./den(:)),0,n_loc_u,n_loc_u);
        Ly2				= Gy * Ly2;
        
        Ly				= (- Ly2  + Ly1)./dy;	clear Lz1 Lz2 Gy
        L			= Lx + Ly + Lz;clear Lx Ly Lz G
    else
        
        L			= Lx + Lz;clear Lx Ly G
        
    end
    
    
    M			=  zeros(n_loc_u,1);
    M(U_index)	=  (-2/model.dt^2 + vec(yeta2(U_index)))./vec(den(U_index)).*vec(v(U_index));
    M			=  spdiags(M,0,n_loc_u,n_loc_u);
    
    A2 			= M - L ;clear L M
    
    A1 			= zeros(n_loc_u,1);
    A1(U_index)	= (1+vec(yeta(U_index)).*model.dt)./(model.dt^2*vec(den(U_index)./v(U_index)));
    A1			= spdiags(A1,0,n_loc_u,n_loc_u);

    A1_inv		= (model.dt^2*den(:)./v(:))./(1+yeta(:).*model.dt);
    A1_inv		= spdiags(A1_inv,0,n_loc_u,n_loc_u);
    
    A3 			= zeros(n_loc_u,1);
    A3(U_index)	= (1-vec(yeta(U_index)).*model.dt).*vec(v(U_index))./(model.dt^2*vec(den(U_index)));
    A3			= spdiags(A3,0,n_loc_u,n_loc_u);

    A3_inv		=  ((den(:)./v(:))*model.dt^2)./(1-vec(yeta(U_index)).*model.dt);
    A3_inv		= spdiags(A3_inv,0,n_loc_u,n_loc_u);
    
    
elseif isempty(den) && ~isempty(thomsen)
    % if labindex==1
        % disp('Acoustic anisotropic');
    % end
    
    epsilon = thomsen.epsilon(:);
    delta = thomsen.delta(:);
    theta = thomsen.theta(:);
        
    A1_inv		= (model.dt^2./v(:))./(1+yeta(:).*model.dt);
    A1_inv		= opDiag(A1_inv);

    A3 			= zeros(n_loc_u,1);
    A3(U_index)	= (1-vec(yeta(U_index)).*model.dt).*vec(v(U_index))./(model.dt^2);
    A3			= opDiag(A3);

    A1 			= zeros(n_loc_u,1);
    A1(U_index)	= (1+vec(yeta(U_index)).*model.dt).*(vec(v(U_index)))./model.dt^2;
    A1			= opDiag(A1);

    A3_inv		=  (model.dt^2./v(:))./(1-yeta(:).*model.dt);
    A3_inv		= opDiag(A3_inv);
    
    M				=  zeros(n_loc_u,1);
    M(U_index)		=  (-2/model.dt^2 + vec(yeta2(U_index))).*vec(v(U_index));
    M				=  opDiag(M);

    n_loc = [nz_loc_u, nx_loc_u,ny_loc_u];
    a = model.o;
    b = model.d.*(n_loc-1);
    eps = 1e-20;
    nz = n_loc(1);
    nx = n_loc(2);
    
    % 1D wavenumber vectors
    kz = 2*pi/(b(1) - a(1))*[eps:nz/2-1 eps -nz/2+1:-1];
    kx = 2*pi/(b(2) - a(2))*[eps:nx/2-1 eps -nx/2+1:-1];
    
    % 2d wavenumber fields
    Kz = kron(ones(1,nx),kz');
    Kx = kron(ones(nz,1),kx);
    Ksqr = Kx.^2+Kz.^2;
    K1 = -Kx.^2;
    K2 = -Kz.^2;
    K3 = -Kx.*Kz;
    K4 = Kx.^4./Ksqr;
    K5 = Kz.^4./Ksqr;
    K6 = Kx.^3.*Kz./Ksqr;
    K7 = Kx.^2.*Kz.^2./Ksqr;
    K8 = Kx.*Kz.^3./Ksqr;
    
    
    % Operators for anisotropic Laplacian
    cxx = ((1+2*epsilon).*cos(theta).^2 + (1+2*delta-2*epsilon).*sin(theta).^2);
    czz = ((1+2*epsilon).*sin(theta).^2 + (1+2*delta-2*epsilon).*cos(theta).^2);
    cxz = -2*(2*epsilon-delta).*sin(2*theta);
    cxxxx = 2*(delta-epsilon).*sin(theta).^4;
    czzzz = 2*(delta-epsilon).*cos(theta).^4;
    cxxxz = 8*(delta-epsilon).*sin(theta).^3.*cos(theta);
    cxxzz = 3*(delta-epsilon).*sin(2*theta).^2;
    cxzzz = 8*(delta-epsilon).*sin(theta).*cos(theta).^3;

    F2 = opDFT2(nz,nx);
    
    D1 = opDiag(cxx)*real(F2'*opDiag(K1)*F2);
    D2 = opDiag(czz)*real(F2'*opDiag(K2)*F2);
    D3 = opDiag(cxz)*real(F2'*opDiag(K3)*F2);
    D4 = opDiag(cxxxx)*real(F2'*opDiag(K4)*F2);
    D5 = opDiag(czzzz)*real(F2'*opDiag(K5)*F2);
    D6 = opDiag(cxxxz)*real(F2'*opDiag(K6)*F2);
    D7 = opDiag(cxxzz)*real(F2'*opDiag(K7)*F2);
    D8 = opDiag(cxzzz)*real(F2'*opDiag(K8)*F2);

    L = D1+D2+D3+D4+D5+D6+D7+D8;
    
    clear D1 D2 D3 D4 D5 D6 D7 D8
    clear cxx czz cxz cxxxx czzzz cxxxz cxxzz cxzzz
    
    % Operators for dA epsilon
    dcxx = 2*cos(theta).^2 - 2*sin(theta).^2;
    dczz = 2.*sin(theta).^2 - 2*cos(theta).^2;
    dcxz = -4*sin(2*theta);
    dcxxxx = -2*sin(theta).^4;
    dczzzz = -2*cos(theta).^4;
    dcxxxz = -8*sin(theta).^3.*cos(theta);
    dcxxzz = -3*sin(2*theta).^2;
    dcxzzz = -8*sin(theta).*cos(theta).^3;
    
    dD1 = opDiag(dcxx)*real(F2'*opDiag(K1)*F2);
    dD2 = opDiag(dczz)*real(F2'*opDiag(K2)*F2);
    dD3 = opDiag(dcxz)*real(F2'*opDiag(K3)*F2);
    dD4 = opDiag(dcxxxx)*real(F2'*opDiag(K4)*F2);
    dD5 = opDiag(dczzzz)*real(F2'*opDiag(K5)*F2);
    dD6 = opDiag(dcxxxz)*real(F2'*opDiag(K6)*F2);
    dD7 = opDiag(dcxxzz)*real(F2'*opDiag(K7)*F2);
    dD8 = opDiag(dcxzzz)*real(F2'*opDiag(K8)*F2);
    
    dL = dD1+dD2+dD3+dD4+dD5+dD6+dD7+dD8;
    
    clear dD1 dD2 dD3 dD4 dD5 dD6 dD7 dD8
    clear dcxx dczz dcxz dcxxxx dczzzz dcxxxz dcxxzz dcxzzz
    
    clear Kz Kx Kz2 Ksqr K1 K2 K3 K4 K5 K6 K7 K8
    
    
    % Laplacian
    dA2ani =-dL;
    A2 = -L + M; 
    clear L F2 M dL
    
else
    error('This mode is not currently supported');
end

% dA for jacobian. 
% Second time derivative, PML don't matter as it will be cut for the
% gradient
dA1 = 1/model.dt^2;
dA2 = -2/model.dt^2;
dA3=dA1;


end