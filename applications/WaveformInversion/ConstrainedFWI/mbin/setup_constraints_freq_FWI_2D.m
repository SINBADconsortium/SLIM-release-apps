function [P,LB,UB] = setup_constraints_freq_FWI_2D(constraint,model,comp_grid,params,m,freq_batch_nr)
m=m(:);
nr_constraints =constraint.use_cardinality_2+constraint.use_bounds+constraint.use_min_smooth+constraint.use_nuclear+constraint.use_rank+constraint.use_TV+constraint.use_cardinality+constraint.use_TD_bounds+constraint.use_cvxbody+constraint.use_subspace;
counter=1;

if nr_constraints==0
    %P{1}=@(input) 1*input;
    error('no constraints selected')   
else 
    for i=1:nr_constraints
        
        %% bound constraints
        
        if isfield(constraint,'water_depth')==1
            LB_glob = ones(comp_grid.n)*(1e6 * ( 1./(constraint.v_max.^2) ));
            UB_glob = ones(comp_grid.n)*(1e6 * ( 1./(constraint.v_min.^2) ));
            %UB_rel  = 1e6 * ( 1./(1e3./sqrt(constraint.m_start)-constraint.v_rel_min).^2) ;
            %LB_rel  = 1e6 * ( 1./(1e3./sqrt(constraint.m_start)+constraint.v_rel_plus).^2) ;
            water_bottom_index = floor(constraint.water_depth/comp_grid.d(1));
            if water_bottom_index==0;
                water_bottom_index=1;
            end
            LB_water = ones(comp_grid.n).*(1e6 * ( 1./(constraint.water_max.^2) ));
            UB_water = ones(comp_grid.n).*(1e6 * ( 1./(constraint.water_min.^2) ));
            LB_water(water_bottom_index+1:end,:)=LB_glob(water_bottom_index+1:end,:);
            UB_water(water_bottom_index+1:end,:)=UB_glob(water_bottom_index+1:end,:);
            
            %tightest combination
            %LB=max(LB_glob(:),max(LB_rel(:),LB_water(:)));
            LB=max(LB_glob(:),LB_water(:));
            %UB=min(UB_glob(:),max(UB_rel(:),UB_water(:)));
            UB=min(UB_glob(:),UB_water(:));
        else
            LB = ones(length(m),1)*(1e6 * ( 1./(constraint.v_max.^2) ));
            UB = ones(length(m),1)*(1e6 * ( 1./(constraint.v_min.^2) ));
        end
        %figure;imagesc(reshape(LB,comp_grid.n));
        %figure;imagesc(reshape(UB,comp_grid.n));
        %pause(1)
        if constraint.use_bounds == 1
            P{counter}            = @(input) boundProject(input,LB,UB);
            counter               = counter+1;
            constraint.use_bounds = 0;
        end
        
        %% Nuclear norm constraints
        if constraint.use_nuclear == 1
            error('nuclear norm temporarily not supported')
%             if freq_batch_nr==1; PN = proj_nuclear(constraint.NN_initial*norm(svd(reshape(m,comp_grid.n)),1)); end
%             if freq_batch_nr>1;  PN = proj_nuclear(constraint.NN_increment*norm(svd(reshape(m,comp_grid.n)),1)); end
%             P{counter}             = @(input) vec(PN(reshape(input,comp_grid.n),1));
%             counter                = counter+1;
%             constraint.use_nuclear = 0;
        end
        
        %% Rank constraints
        if constraint.use_rank == 1
            if freq_batch_nr==1; rank = constraint.r_initial; end;
            if freq_batch_nr>1;  rank = floor(constraint.r_initial+freq_batch_nr*constraint.r_increment); end;
            P{counter} = @(input) vec(proj_Rank(reshape(input,comp_grid.n),rank));
            counter                = counter+1;
            constraint.use_rank = 0;
        end
        
        %% Total-variation constraints
        if constraint.use_TV == 1
            h1=comp_grid.d(1);
            h2=comp_grid.d(2);
            [D,~] = getDiscreteGrad(comp_grid.n(1),comp_grid.n(2),h1,h2,ones(size(m,1),1));
            
            if strcmp(constraint.TV_units,'m/s')==1
                if freq_batch_nr==1; TV = constraint.TV_initial.*sum(abs(D*(1e3./sqrt(P{1}(m(:)))))); end;
                if freq_batch_nr>1;  TV = constraint.TV_increment.*sum(abs(D*(1e3./sqrt(P{1}(m(:)))))); end;
            elseif strcmp(constraint.TV_units,'s^2/m^2')==1
                if freq_batch_nr==1; TV = constraint.TV_initial.*sum(abs(D*m(:))); end;
                if freq_batch_nr>1;  TV = constraint.TV_increment.*sum(abs(D*m(:))); end;
            end
            
            if constraint.ADMM_options.adjust_rho==0
                R=chol(speye(length(m))+constraint.ADMM_options.rho*(D'*D));
            else
                R=[];
            end
            funProj = @(input) NormL1_project(input,TV);
            if strcmp(constraint.TV_units,'m/s')==1
                P{counter} = @(input) (1e6./((Compute_projection_ADMM((1e3./sqrt(P{1}(input))),D,funProj,constraint.ADMM_options,R)).^2));
            elseif strcmp(constraint.TV_units,'s^2/m^2')==1
                P{counter} = @(input)  Compute_projection_ADMM(input,D,funProj,constraint.ADMM_options,R);
            end
            
            counter                = counter+1;
            constraint.use_TV = 0;
        end
        
        %% Cardinality constraints
        if constraint.use_cardinality==1
            if freq_batch_nr==1;
                [S]=get_TD_operator(m,comp_grid,constraint);
                %card = constraint.card_initial*nnz(S*m(:));
                card = constraint.card_initial*comp_grid.n(2);
            end;
            if freq_batch_nr>1;
                [S]=get_TD_operator(m,comp_grid,constraint);
                card = constraint.card_increment*constraint.card_initial*comp_grid.n(2);
            end
            if constraint.ADMM_options.adjust_rho==0
                R=chol(speye(length(m))+constraint.ADMM_options.rho*(S'*S));
            else
                R=[];
            end
            card=round(card);
            funProj = @(input) proj_card(input,card);
            if strcmp(constraint.card_units,'m/s')==1
                P{counter} = @(input) 1e6./(Compute_projection_ADMM(real(1e3./sqrt(P{1}(input))),S,funProj,constraint.ADMM_options,R).^2);
            elseif strcmp(constraints.card_units,'s^2/m^2')==1
                P{counter} = @(input) Compute_projection_ADMM(input,S,funProj,constraint.ADMM_options,R);
            end
            counter                = counter+1;
            constraint.use_cardinality=0;
        end
        
          %% Cardinality constraints (another type)
        if constraint.use_cardinality_2==1
            if freq_batch_nr==1;
                [S]=get_TD_operator(m,comp_grid,constraint.card_operator_2);
                %card = constraint.card_initial*nnz(S*m(:));
                card =2*comp_grid.n(1);
            end;
            if freq_batch_nr>1;
                [S]=get_TD_operator(m,comp_grid,constraint);
                card = constraint.card_increment*constraint.card_initial*comp_grid.n(2);
            end
            if constraint.ADMM_options.adjust_rho==0
                R=chol(speye(length(m))+constraint.ADMM_options.rho*(S'*S));
            else
                R=[];
            end
            card=round(card);
            funProj = @(input) proj_card(input,card);
            if strcmp(constraint.card_units,'m/s')==1
                P{counter} = @(input) 1e6./(Compute_projection_ADMM(real(1e3./sqrt(P{1}(input))),S,funProj,constraint.ADMM_options,R).^2);
            elseif strcmp(constraints.card_units,'s^2/m^2')==1
                P{counter} = @(input) Compute_projection_ADMM(input,S,funProj,constraint.ADMM_options,R);
            end
            counter                = counter+1;
            constraint.use_cardinality=0;
        end
        
        %% Transform domain bound-constraints
        if constraint.use_TD_bounds==1
            
            if freq_batch_nr==1;
                [S]=get_TD_operator(m,comp_grid,constraint);
            end;
            if freq_batch_nr>1;
                [S]=get_TD_operator(m,comp_grid,constraint);
                constraint.TD_LB = constraint.TD_LB_increment*constraint.TD_LB;
                constraint.TD_UB = constraint.TD_UB_increment*constraint.TD_UB;
            end
            funProj = @(input) boundProject(input,constraint.TD_LB.*ones(size(S,1),1),constraint.TD_UB.*ones(size(S,1),1));
            if constraint.ADMM_options.adjust_rho==0
                R=chol(speye(length(m))+constraint.ADMM_options.rho*(S'*S));
            else
                R=[];
            end
            if strcmp(constraint.TDbounds_units,'m/s')==1
                P{counter} = @(input) 1e6./(Compute_projection_ADMM(real(1e3./sqrt(P{1}(input))),S,funProj,constraint.ADMM_options,R).^2);
            elseif strcmp(constraint.TDbounds_units,'s^2/m^2')==1
                P{counter} = @(input) Compute_projection_ADMM(input,S,funProj,constraint.ADMM_options,R);
            end
            
            counter                = counter+1;
            constraint.use_TD_bounds=0;
        end
        
        %% Fourier based minimum smoothness constraints
        if constraint.use_min_smooth == 1
            cmin        = min(1e3./sqrt(m(:)));
            fmax        = max(model.freq(params.freq_index));
            k           = fmax / cmin;
            HPF         = getHPF(comp_grid, constraint, k);
            LPF         = getLPF(comp_grid, constraint, k);
            
            P{counter} = @(input) LPF*input;
            counter = counter+1;
            constraint.use_min_smooth = 0;
        end
        
        %% Projection onto a convex body (sort of...)
        if constraint.use_cvxbody == 1
            P{counter} = @(input) projection_cvxbody(P{1}(input),constraint.cvxb_threshold,comp_grid);
            counter = counter+1;
            constraint.use_cvxbody = 0;
        end
        
        %% Projection onto a subspace
        if constraint.use_subspace ==1 
           %P{counter} = @(input) vec(Proj_subspace(reshape(input,comp_grid.n),constraint.U,1));
           P{counter} = @(input) real(vec(Proj_subspace(input,constraint.U,0)));
           counter = counter+1;
           constraint.use_subspace = 0; 
        end
    end
end
