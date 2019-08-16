function P = setup_constraints_freq_FWI_3D(constraint,model,m,freq_batch_nr)
% Assumes input is in m/s, output is also m/s.
% Script designed for 3D frequency domain FWI
%
% This script it work in progress, as not all constraints that are
% available in 2D are ready to be used in 3D.
%
% Bas Peters, March 2016.

if isfield(model,'n')
    n=model.n;
elseif isfield(model,'nt')
    n=model.nt;
else
    error('no information about the grid geometry is available, specify model.n or model.nt')
end
if isfield(model,'d')
    d=model.d;
elseif isfield(model,'dt')
    d=model.dt;
else
    error('no information about the grid geometry is available, specify model.n or model.nt')
end

nr_constraints = constraint.use_bounds +  constraint.use_TV ;
counter=1;

if nr_constraints==0
    error('no constraints selected')
else
    for i=1:nr_constraints
        
        if constraint.use_bounds == 1
             %% bound constraints
        LB_glob = ones(n).*(constraint.v_min);
        UB_glob = ones(n).*(constraint.v_max);
        
        %UB_rel  = 1e6 * ( 1./(1e3./sqrt(constraint.m_start)-constraint.v_rel_min).^2) ;
        %LB_rel  = 1e6 * ( 1./(1e3./sqrt(constraint.m_start)+constraint.v_rel_plus).^2) ;
        
        if isfield(constraint,'water_depth')==1
            water_bottom_index = floor(constraint.water_depth/d(1));
            if water_bottom_index==0;
                water_bottom_index=1;
            end
            LB_water = ones(n).*(constraint.water_min);
            UB_water = ones(n).*(constraint.water_max);
            
            LB_water(water_bottom_index+1:end,:,:)=LB_glob(water_bottom_index+1:end,:,:);
            UB_water(water_bottom_index+1:end,:,:)=UB_glob(water_bottom_index+1:end,:,:);
            
            %tightest combination
            LB=max(LB_glob(:),LB_water(:));
            UB=min(UB_glob(:),UB_water(:));
            
        else
        LB = ones(n).*(constraint.v_min);
        UB = ones(n).*(constraint.v_max);
        end
        
            P{counter}            = @(input) boundProject(input,LB(:),UB(:));
            counter               = counter+1;
            constraint.use_bounds = 0;
        end
        
        %% Total-variation constraints
        if constraint.use_TV == 1
            [Dx_3D,Dy_3D,Dz_3D] = getDiscreteGrad3D(n,d);
            D=[Dx_3D;Dy_3D;Dz_3D];
            %h=comp_grid.d(1);
            %[D,~] = getDiscreteGrad(comp_grid.n(1),comp_grid.n(2),h,ones(size(m,1),1));
            %D =@(input) vec(gradient_op3d((reshape(input,n)));
            
            
            if freq_batch_nr==1; TV = constraint.TV_initial.*sum(abs(D*m(:))); end;
            if freq_batch_nr>1;  TV = constraint.TV_increment.*sum(abs(D*m(:))); end;

        
        if constraint.ADMM_options.adjust_rho==0
            R=chol(speye(length(m))+constraint.ADMM_options.rho*(D'*D));
        else
            R=[];
        end
        funProj = @(input) NormL1_project(input,TV);
        P{counter} = @(input)  Compute_projection_ADMM(input,D,funProj,constraint.ADMM_options,R);
        
        counter                = counter+1;
        constraint.use_TV = 0;
    end
    
end
end
