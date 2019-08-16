function P = setup_constraintsTS(constraint,model,v,freq_batch_nr,maxv)

	nr_constraints = constraint.use_bounds+constraint.use_min_smooth+constraint.use_nuclear+constraint.use_rank+constraint.use_TV  ;

	counter=1;

	if nr_constraints==0
		P{1}=@(input) 1*input;
	else

		for i=1:nr_constraints
			
			%% bound constraints
			
			if nargin==6
				v_max=maxv;
			else
				v_max=max(sqrt(1./v(:)))+constraint.v_max_plus;
			end
			v_min=min(sqrt(1./v(:)))-constraint.v_min_minus;
			v_max=1./v_max.^2;
			v_min=1./v_min.^2;
			LB = ones(length(v),1)*v_max;
			UB = ones(length(v),1)*v_min;
			
			if constraint.use_bounds == 1
				P{counter}            = @(input) boundProject(input,LB,UB);
				counter               = counter+1;
				constraint.use_bounds = 0;
			end
			
			%% minimum smoothness projector
			if constraint.use_min_smooth == 1
				cmin        = min(1e3*sqrt(1./v(:)));
				k           = 10+1e3*model.f0 / cmin;
				LPF         = getLPF(model, constraint, k);
				
				P{counter} = @(input) LPF*input;
				counter = counter+1;
				constraint.use_min_smooth = 0;
			end
			
			  %% Total-variation constraints
				if constraint.use_TV == 1
					h1=model.d(1);
                    h2=model.d(2);
					[D,~] = getDiscreteGrad(model.n(1),model.n(2),h1,h2,ones(size(v,1),1));
					
					if strcmp(constraint.TV_units,'m/s')==1
						if freq_batch_nr==1; TV = constraint.TV_initial.*sum(abs(D*(1e3./sqrt(P{1}(v(:)))))); end;
						if freq_batch_nr>1;  TV = constraint.TV_increment.*sum(abs(D*(1e3./sqrt(P{1}(v(:)))))); end;
					elseif strcmp(constraint.TV_units,'s^2/m^2')==1
						if freq_batch_nr==1; TV = constraint.TV_initial.*sum(abs(D*v(:))); end;
						if freq_batch_nr>1;  TV = constraint.TV_increment.*sum(abs(D*v(:))); end;
					elseif strcmp(constraint.TV_units,'s^2/km^2')==1
						if freq_batch_nr==1; TV = constraint.TV_initial.*sum(abs(D*v(:))); end;
						if freq_batch_nr>1;  TV = constraint.TV_increment.*sum(abs(D*v(:))); end;
					end

					if constraint.ADMM_options.adjust_rho==0
						R=chol(speye(length(v))+constraint.ADMM_options.rho*(D'*D));
					else
						R=[];
					end
					funProj = @(input) NormL1_project(input,TV);
					if strcmp(constraint.TV_units,'m/s')==1
						P{counter} = @(input) (1e6./((Compute_projection_ADMM((1e3./sqrt(P{1}(input))),D,funProj,constraint.ADMM_options,R)).^2));
					elseif strcmp(constraint.TV_units,'s^2/m^2')==1
						P{counter} = @(input)  Compute_projection_ADMM(input,D,funProj,constraint.ADMM_options,R);
					elseif strcmp(constraint.TV_units,'s^2/km^2')==1
						P{counter} = @(input)  Compute_projection_ADMM(input,D,funProj,constraint.ADMM_options,R);
					end
					
					counter                = counter+1;
					constraint.use_TV = 0;
				end
				
			
			%% Nuclear norm constraints
			if constraint.use_nuclear == 1
				if freq_batch_nr==1; PN = proj_nuclear(constraint.NN_initial*norm(svd(reshape(v,model.n)),1)); end
				if freq_batch_nr>1;  PN = proj_nuclear(constraint.NN_increment*norm(svd(reshape(v,model.n)),1)); end
				P{counter}             = @(input) vec(PN(reshape(input,model.n),1));
				counter                = counter+1;
				constraint.use_nuclear = 0;
			end
			
			%% Rank constraints
			
			if constraint.use_rank == 1
				if freq_batch_nr==1; rrank = constraint.r_initial; end;
				if freq_batch_nr>1;  rrank = floor(constraint.r_initial+freq_batch_nr*constraint.r_increment); end;
				P{counter} = @(input) vec(proj_Rank(reshape(input,model.n),rrank));
				counter                = counter+1;
				constraint.use_rank = 0;
			end
			
		end
	end
end