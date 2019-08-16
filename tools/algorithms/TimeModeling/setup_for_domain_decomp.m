function [vel_dis, fwdpara_dis,dens_dis,ani_dis] = setup_for_domain_decomp(vel,model,dens,ani)
% this function will setup parameters for demain decomperzation.
% Author : XiabgLi edited by Mathias Louboutin
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
% Date : July, 2015
%1, split model
%2, split src and rec location.

%%
% global vel_dis den_dis dvel_dis fwdpara_dis


% --------------------------------------
% create distributed model
% index of labs for domain decomposition
%       ___ ___ ___
%     / 10/13/ 16 / |
%    --- --- ---   /|
%  /___/___/___/ |/ |
% | 1 | 4 | 7 | /|17|
%  --- --- ---|/ |/ |
% | 2 | 5 | 8 | /|18|
%  --- --- ---|/ |/
% | 3 | 6 | 9 | /
%  __________ /
% ---------------------------------------
Num_of_worker=parpool_size();
if Num_of_worker==0
	Num_of_worker=1;
end% number of labs
% disp(['numlabs = ' num2str(Num_of_worker)]);

% Domain decomposition size
if isfield(model,'ddcompx')
	Num_of_model_split = model.ddcompx * model.ddcompz * model.ddcompy;
else
	model.ddcompx=1;
	model.ddcompz =1;
	model.ddcompy=1;
	Num_of_model_split=1;
end


if  Num_of_model_split >  Num_of_worker
    error(['Number of works (',num2str( Num_of_worker),') can not be smaller ' ...
        'than number of model splits (',num2str( Num_of_model_split),'); Need to be N ' ...
        'times larger, then you do not waste works.']);
end

nx		= model.n(2);
ny		= model.n(3);
nz		= model.n(1);

ngx	= round(linspace(0,nx,model.ddcompx+1));
ngz	= round(linspace(0,nz,model.ddcompz+1));
ngy	= round(linspace(0,ny,model.ddcompy+1));


[model.mmz,model.mmx,model.mmy]=odn2grid(model.o,model.d,model.n);
% create some datas for difference works
fwdpara_dis = Composite();
vel_dis     = Composite();
% mmx mmy mmz for each model split
vel = reshape(vel,nz,nx,ny);
if ~isempty(ani)
    ani.epsilon=reshape(ani.epsilon,nz,nx,ny);
    ani.delta=reshape(ani.delta,nz,nx,ny);
    ani.theta=reshape(ani.theta,nz,nx,ny);
    ani_dis = Composite();
end
if ~isempty(dens)
    dens= reshape(dens,nz,nx,ny);
    dens_dis = Composite();
end

for labi = 1:Num_of_worker
    model_parition_index	= mod(labi-1,Num_of_model_split)+1;
    [ddcompz_id,ddcompx_id,ddcompy_id]	= ind2sub([model.ddcompz,model.ddcompx,model.ddcompy],model_parition_index);
    num_of_extra_point		= max([2,model.space_order/2]); % extra points needs for domain decompzation
    % find x,y,z grid coordinate of model that on this worker.
    
    if model.ddcompx == 1
        mmx_loc = model.mmx;
        u_idx   = 1:length(mmx_loc);
        x_id	= u_idx;
    else
        switch ddcompx_id
            case 1
                x_id	= ngx(ddcompx_id)+1:ngx(ddcompx_id+1)+num_of_extra_point;
                mmx_loc = model.mmx(x_id);
                u_idx   = 1:(length(mmx_loc)-num_of_extra_point);
            case model.ddcompx
                x_id	= ngx(ddcompx_id)-num_of_extra_point+1:ngx(ddcompx_id+1);
                mmx_loc = model.mmx(x_id);
                u_idx   = (num_of_extra_point+1):length(mmx_loc);
            otherwise
                x_id	= ngx(ddcompx_id)-num_of_extra_point+1:ngx(ddcompx_id+1)+num_of_extra_point;
                mmx_loc = model.mmx(x_id);
                u_idx   = num_of_extra_point+1:(length(mmx_loc)-num_of_extra_point);
        end
    end
    if model.ddcompz == 1
        mmz_loc = model.mmz;
        u_idz   = 1:length(mmz_loc);
        z_id	= u_idz;
    else
        switch ddcompz_id
            case 1
                z_id	= ngz(ddcompz_id)+1:ngz(ddcompz_id+1)+num_of_extra_point;
                mmz_loc = model.mmz(z_id);
                u_idz   = 1:(length(mmz_loc)-num_of_extra_point);
            case model.ddcompz
                z_id	= ngz(ddcompz_id)-num_of_extra_point+1:ngz(ddcompz_id+1);
                mmz_loc = model.mmz(z_id);
                u_idz   = (num_of_extra_point+1):length(mmz_loc);
            otherwise
                z_id	= ngz(ddcompz_id)-num_of_extra_point+1:ngz(ddcompz_id+1)+num_of_extra_point;
                mmz_loc = model.mmz(z_id);
                u_idz   = (num_of_extra_point+1):(length(mmz_loc)-num_of_extra_point);
        end
    end
    if model.ddcompy == 1
        mmy_loc  = model.mmy;
        u_idy   = 1:length(mmy_loc);
        y_id	= u_idy;
    else
        switch ddcompy_id
            case 1
                y_id	= ngy(ddcompy_id)+1:ngy(ddcompy_id+1)+num_of_extra_point;
                mmy_loc = model.mmy(y_id);
                u_idy   = 1:(length(mmy_loc)-num_of_extra_point);
            case model.ddcompy
                y_id	= ngy(ddcompy_id)-num_of_extra_point+1:ngy(ddcompy_id+1);
                mmy_loc = model.mmy(y_id);
                u_idy   = (num_of_extra_point+1):length(mmy_loc);
            otherwise
                y_id	= ngy(ddcompy_id)-num_of_extra_point+1:ngy(ddcompy_id+1)+num_of_extra_point;
                mmy_loc = model.mmy(y_id);
                u_idy   = (num_of_extra_point+1):(length(mmy_loc)-num_of_extra_point);
        end
    end
    
    slx_dis		= model.xsrc(:);
    sly_dis     = model.ysrc(:);
    slz_dis     = model.zsrc(:);
    
    % find rec location on this work
    fwdpara_dis{labi}	= struct('x_id_loc',x_id					...
        ,'y_id_loc',y_id					...
        ,'z_id_loc',z_id					...
        ,'u_idx',u_idx						...
        ,'u_idy',u_idy						...
        ,'u_idz',u_idz						...
        ,'np_extra',num_of_extra_point 		...
        ,'mmx_loc',mmx_loc					...
        ,'mmy_loc',mmy_loc					...
        ,'mmz_loc',mmz_loc                  ...
        ,'num_model',Num_of_model_split     ...
        ,'slx_band',slx_dis 				...
        ,'sly_band',sly_dis 				...
        ,'slz_band',slz_dis					...
        ,'ddcompx_id', ddcompx_id			...
        ,'ddcompy_id', ddcompy_id			...
        ,'ddcompz_id', ddcompz_id			...
        ,'o',model.o						...
        ,'d',model.d						...
        ,'n',model.n						...
        ,'nb',model.nb						...
        ,'space_order',model.space_order    ...
        ,'mmx',model.mmx					...
        ,'mmy',model.mmy					...
        ,'mmz',model.mmz    				...
        ,'ddcompx',model.ddcompx			...
        ,'ddcompy',model.ddcompy			...
        ,'ddcompz',model.ddcompz			...
        ,'freesurface',model.freesurface	...
        ,'dt',model.dt					    ...
        ,'T',model.T);
    
		% nloc=[length(u_idz) length(u_idx)]
    
        vel_dis{labi} 		= vel(z_id,x_id,y_id);
        
        if ~isempty(ani)
            eps_dis=ani.epsilon(z_id,x_id,y_id);
            delta_dis=ani.delta(z_id,x_id,y_id);
            theta_dis=ani.theta(z_id,x_id,y_id);
            aniloc=struct('epsilon',eps_dis,'delta',delta_dis,'theta',theta_dis);
            ani_dis{labi}       = aniloc;
        end
        
        if ~isempty(dens)
             dens_dis{labi} 		= dens(z_id,x_id,y_id);
        end
    
end

if isempty(ani)
    ani_dis=[];
end

if isempty(dens)
    dens_dis=[];
end
    
end