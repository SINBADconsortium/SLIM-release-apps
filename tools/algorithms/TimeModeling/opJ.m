classdef opJ < opSpot
    
    % opJ -- 3D Jacobian operator.
    %
    % Linear operator that (implicitely) multiplies a model perturbation 'dm'.
    % The output of a multiplication operation
    % is thus linearized data, '\delta d', measured at the level of the receivers in time
    % domain. Transpose multiplication mode will back propagate receiver
    % data residual '\delta d' to a model perturbation 'dm'.
    %
    %
    %
    % INPUTS to the constructor
    %
    %   'model'  : stucture containing the model parameters
    %       'm'  : square slowness [s^2/m^2]
    %       'q'  : source for forward propagation
    %       'd'  : True data
    %      'dens': Density
    %     'anis' : anisotropy parameters (epsilon,delta,theta in a
    %     structure)
    %
    % INPUT to multiplication
    %
    %   'x':    direct mode:  model perturbation
    %               nz*nx*nt
    %
    %           adjoint mode: time domain vectorized data residual
    %
    %
    % OUTPUT from multiplication
    %
    %   'y':    direct mode [cell array]:
    %              y : vectorized linearized data
    %
    %           adjoint mode [vector]: Model perturbation
    %
    %
    % Author: Mathias Louboutin from Xiangli's
	%         Seismic Laboratory for Imaging and Modeling
	%         Department of Earth, Ocean, and Atmosperic Sciences
	%         The University of British Columbia
    % Date: July 31, 2014
    % edited by Philipp Witte, July 2015
    %
    
    properties
        
        v=[];
        model=[];
        q=[];
        d=[];
        thomsen=[];
        dens=[];
        
    end
    
    methods
        %% Constructor
        function op=opJ(v,model,q,d,dens,thomsen)
            
            Np=length(model.mmz_loc)*length(model.mmx_loc)*length(model.mmy_loc);
            
            %Creation of the SPOT operator
            % Number of receiver and time samples for size
            nrec=length(model.rec_idx);
            if nrec==0
                nrec=1;
            end
            subT2=length(model.NyqT);
            % Operator
            op=op@opSpot('opj',nrec*subT2,Np);
            
            % Model,velocity,density, anisotropy and data
            op.model=model;
            op.v=v;
            op.q=q;
            op.d=d;
            op.thomsen = thomsen;
            op.dens=dens;
        end
    end
    
    methods ( Access = protected )
        %Multiplication function
        function y=multiply(op,x,mode)
            
            %**************************************************************
            % DIRECT MODE (mode=1)
            %
            % INPUT: perturbation in squared slowness, dm
            %
            %   Vector of size N, where N is the number of physical grid
            %   points.
            %
            % OUTPUT: data misfit (d-F) (linear approximation) where d is
            %   the receivers data and F the forward modelled receivers
            %   data.
            %
            %**************************************************************
            % TRANSPOSE MODE (mode=2)
            %
            % INPUT: data misfit, (d-F)
            %
            %   Vector of size nrec*nt.
            %
            % OUTPUT: perturbation, dm (linear approximation)
            %
            %**************************************************************
            
            %Checking of the multiplication mode
            adjoint_mode=false;
            if mode==2
                adjoint_mode=true;
            end
            
            % Number of time steps
            nt=length((0:op.model.dt:op.model.T));
            % Time axis at time steping rate
            fullT=0:op.model.dt:op.model.T;
            % Shot record time axis
            NyqT=op.model.NyqT;
            
            % Checking if simultenaous source
            if isfield(op.model,'multi')
                nsrc=op.model.multi;
                op.q=reshape(op.q,nsrc,length(0:op.model.dt:op.model.T));
            else
                nsrc=1;
            end
            
            if adjoint_mode
                
                nr=length(op.model.rlx_loc);
                % Receiver projector
                Pr	= Rec3D(op.model);
                if ~isempty(Pr)
	                % Time projector from shot record time axis to full
	                % time axis
	                Pt=kron(speye(nr),opLInterp1D(NyqT,fullT));
                    % Reshape input record projected on correct time axis
                    x=reshape(Pt*x(:),nt,nr);
                end
                % Source projector if there is one (no source with
                % model decomposition in all domains but one)
                if length(op.model.slx_loc)>0
                    Ps	= Src3D(op.model,nsrc);
                else
                    Ps=[];
                end
                
                [A1,A1_inv,A2,A3,A3_inv,dA1,dA2,dA3] = generate_time_stepping_stencil_matrix(op.v,op.model,op.dens,op.thomsen);
                % allocate place for gradient.
                y   = zeros(length(op.model.mmz_loc),length(op.model.mmx_loc),length(op.model.mmy_loc));
                
                U3	= zeros(length(op.model.mmz_loc),length(op.model.mmx_loc),length(op.model.mmy_loc));
                
                % computer adjoint wavefield
                V1	= zeros(length(op.model.mmz_loc),length(op.model.mmx_loc),length(op.model.mmy_loc));
                V2	= zeros(length(op.model.mmz_loc),length(op.model.mmx_loc),length(op.model.mmy_loc));
                V3	= zeros(length(op.model.mmz_loc),length(op.model.mmx_loc),length(op.model.mmy_loc));

                
                % Put matrices together to save time
				if ~isempty(Pr)
                	Pr = A1_inv'*Pr';
				end
                A2 = A1_inv'*A2' ;
                A3 = A1_inv'*A3';
              
                
                % Clear unnecessary matrices
                clear A3_inv A1_inv Ps A1;
                % Begin time loop
                for ti=nt:-1:1
                        % -------------- load check point per time step ---------------------
                        if length(find(op.model.nsub==ti))>0
                            if strcmp(op.model.save,'disk')
                                U3=myReadWf(ti,op.model.numsave,labindex);
                            else
                                U3=reshape(op.d(find(op.model.nsub==ti),:),size(U3));
                            end
                    	end
	                    if ~isempty(Pr)
	                        u_adj=x(ti,:)';
                        
	                        if ti == nt
	                            V3(1:end)   = full((Pr * u_adj));
	                        elseif ti == nt-1
	                            V3(1:end)   = Pr * u_adj - A2*V2(:);
	                        else
	                            V3(1:end) 	= Pr * u_adj -A2*V2(:) - A3*V1(:);
	                        end
                        
	                    else
	                        if ti == nt
	                            V3(1:end)   = 0*V3(1:end) ;
	                        elseif ti == nt-1
	                            V3(1:end)   = - A2*V2(:);
	                        else
	                            V3(1:end) 	= - A2*V2(:) - A3*V1(:);
	                        end
                        
	                    end
                    
                    % ------------ domain decomposition send data around ---------
	                    if op.model.num_model ~=1
	                        if op.model.ddcompz~=1
	                            if op.model.ddcompz_id ~= 1
	                                labSend(V3(1:op.model.np_extra*2,:,:),labindex-1)
	                            end
	                            if op.model.ddcompz_id ~= op.model.ddcompz
	                                V3(end-2*op.model.np_extra+1:end,:,:) = V3(end-2*op.model.np_extra+1:end,:,:) + labReceive(labindex + 1);
	                                labSend(V3(end-2*op.model.np_extra+1:end,:,:),labindex+1)
	                            end
	                            if op.model.ddcompz_id ~= 1
	                                V3(1:op.model.np_extra*2,:,:) = labReceive(labindex-1);
	                            end
	                        end
                        
	                        if op.model.ddcompx~=1
	                            if op.model.ddcompx_id ~= 1
	                                labSend(V3(:,1:op.model.np_extra*2,:),labindex-op.model.ddcompz)
	                            end
	                            if op.model.ddcompx_id ~= op.model.ddcompx
	                                V3(:,((end-2*op.model.np_extra+1):end),:) = V3(:,((end-2*op.model.np_extra+1):end),:) + labReceive(labindex + op.model.ddcompz);
	                                labSend(V3(:,(end-2*op.model.np_extra+1:end),:),labindex+op.model.ddcompz)
	                            end
	                            if op.model.ddcompx_id ~= 1
	                                V3(:,1:2*op.model.np_extra,:) =  labReceive(labindex-op.model.ddcompz);
	                            end
	                        end
                        
                        
	                        if op.model.ddcompy~=1
	                            if op.model.ddcompy_id ~= 1
	                                labSend(V3(:,:,1:2*op.model.np_extra),labindex-op.model.ddcompz*op.model.ddcompx)
	                            end
	                            if op.model.ddcompy_id ~= op.model.ddcompy
	                                V3(:,:,(end-2*op.model.np_extra+1:end)) = V3(:,:,(end-2*op.model.np_extra+1):end) + labReceive(labindex + op.model.ddcompz*op.model.ddcompx);
	                                labSend(V3(:,:,(end-2*op.model.np_extra+1):end),labindex+op.model.ddcompz*op.model.ddcompx)
	                            end
	                            if op.model.ddcompy_id ~= 1
	                                V3(:,:,1:2*op.model.np_extra) = labReceive(labindex-op.model.ddcompz*op.model.ddcompx);
	                            end
	                        end
                        
	                    end
                    % --------------- end of domain decomposition -----------------------
                    
                    % Compute dA^T V and update gradient when wanted
	                    if(length(find(op.model.nconv==(ti)))>0)  
	                        if ti == 1
	                            U_temp =dA3'*V1(:);
	                        elseif ti == 2
	                            U_temp = dA2'*V2(:) +dA3'*V1(:);
	                        else
	                            U_temp =dA1'*V3(:) + dA2'*V2(:) + dA3'*V1(:);
	                        end
	                        y = y - reshape(U_temp,size(V3)).*U3;
	                    end
                    % Free surface
	                    if op.model.freesurface
	                        for loc=1:op.model.space_order/2
	                            V3(loc,:,:)=-V3(op.model.space_order-loc+2,:,:);
	                        end
	                    end
                    % Swap wavefields
                    V1	= V2;
                    V2	= V3;
                    % End of time loop
                end
                
            else
                
                
                % Number of receivers and receiver projector
                nrecloc=length(op.model.rlx_loc);
                Pr	= Rec3D(op.model);
                
                % Source projector if there is one (no source with
                % model decomposition in all domains but one)
                if length(op.model.slx_loc)>0
                    Ps	= Src3D(op.model,nsrc);
                else
                    Ps=[];
                end
                
                % Time projector from full time axis to shot record
                % time axis
                Pt=kron(speye(nrecloc),opLInterp1D(fullT,NyqT));
                
                
                % generate time stepping stencil matrix
                [~,A1_inv,A2,A3,~,dA1,dA2,dA3] = generate_time_stepping_stencil_matrix(op.v,op.model,op.dens,op.thomsen);
                % combine some matrix to save time
                A2 = A1_inv * A2;
                A3 = A1_inv * A3;
                
                if ~isempty(Ps)
                    Ps = A1_inv * Ps';
                end
                
                % allocate wavefield snapshot
                U1	= zeros(length(op.model.mmz_loc),length(op.model.mmx_loc),length(op.model.mmy_loc));
                U2	= zeros(length(op.model.mmz_loc),length(op.model.mmx_loc),length(op.model.mmy_loc));
                U3	= zeros(length(op.model.mmz_loc),length(op.model.mmx_loc),length(op.model.mmy_loc));
                
                dU1	= zeros(length(op.model.mmz_loc),length(op.model.mmx_loc),length(op.model.mmy_loc));
                dU2	= zeros(length(op.model.mmz_loc),length(op.model.mmx_loc),length(op.model.mmy_loc));
                dU3	= zeros(length(op.model.mmz_loc),length(op.model.mmx_loc),length(op.model.mmy_loc));
                
                % % allocate place for loc shot record
                U_shot_record_temp	= zeros(nt,nrecloc);
                
                % Start loop over time
                for ti=1:nt
                    
                    if isempty(Ps)
                        if ti == 1
                            U3(1:end)   = 0;
                        elseif ti == 2
                            U3(1:end)   =  - A2*U2(:);
                        else
                            U3(1:end) 	= - A2*U2(:) - A3*U1(:);
                        end
                    else
                        if ti == 1
                            U3(1:end)   = full(Ps * op.q(:,ti));
                        elseif ti == 2
                            U3(1:end)   = Ps * op.q(:,ti) - A2*U2(:);
                        else
                            U3(1:end) 	= Ps * op.q(:,ti) - A2*U2(:) - A3*U1(:);
                        end
                    end
					% if labindex==1 && mod(ti,10)==0
					% 	disp([repmat(' ',1,10),'Time=',num2str(ti*op.model.dt),'s, Total=',num2str((nt-1)*op.model.dt),'s']);
					% end
                    % ---------------- domain decomposition send data around -----------------
                    if op.model.num_model ~=1
                        if op.model.ddcompz~=1
                            if op.model.ddcompz_id ~= 1
                                labSend(U3(op.model.np_extra+1:2*op.model.np_extra,:,:),labindex-1)
                            end
                            if op.model.ddcompz_id ~= op.model.ddcompz
                                U3((end-op.model.np_extra+1:end),:,:) = labReceive(labindex + 1);
                                labSend(U3((end-2*op.model.np_extra+1):(end-op.model.np_extra),:,:),labindex+1)
                            end
                            if op.model.ddcompz_id ~= 1
                                U3(1:op.model.np_extra,:,:) = labReceive(labindex-1);
                            end
                        end
                        
                        if op.model.ddcompx~=1
                            if op.model.ddcompx_id ~= 1
                                labSend(U3(:,op.model.np_extra+1:2*op.model.np_extra,:),labindex-op.model.ddcompz)
                            end
                            if op.model.ddcompx_id ~= op.model.ddcompx
                                U3(:,(end-op.model.np_extra+1:end),:) = labReceive(labindex + op.model.ddcompz);
                                labSend(U3(:,(end-2*op.model.np_extra+1):(end-op.model.np_extra),:),labindex+op.model.ddcompz)
                            end
                            if op.model.ddcompx_id ~= 1
                                U3(:,1:op.model.np_extra,:) =  labReceive(labindex-op.model.ddcompz);
                            end
                        end
                        
                        
                        if op.model.ddcompy~=1
                            if op.model.ddcompy_id ~= 1
                                labSend(U3(:,:,op.model.np_extra+1:2*op.model.np_extra),labindex-op.model.ddcompz*op.model.ddcompx)
                            end
                            if op.model.ddcompy_id ~= op.model.ddcompy
                                U3(:,:,(end-op.model.np_extra+1:end)) = labReceive(labindex + op.model.ddcompz*op.model.ddcompx);
                                labSend(U3(:,:,(end-2*op.model.np_extra+1):(end-op.model.np_extra)),labindex+op.model.ddcompz*op.model.ddcompx)
                            end
                            if op.model.ddcompy_id ~= 1
                                U3(:,:,1:op.model.np_extra) = labReceive(labindex-op.model.ddcompz*op.model.ddcompx);
                            end
                        end
                        
                    end
                    % --------------- end of domain decomposition -----------------------
                    
                    
                    if ti == 1
                        U_temp = (dA3*U3(:)).* (x(:));
                    elseif ti == 2
                        U_temp = (dA2*U2(:) + dA3*U3(:)).* (x(:));
                    else
                        U_temp = (dA1*U1(:) + dA2*U2(:) + dA3*U3(:)).* (x(:));
                    end
                    
                    % do A-1 * U_temp
                    
                    if ti == 1
                        dU3(1:end) = (A1_inv* U_temp(:));
                    elseif ti == 2
                        dU3(1:end) = (A1_inv*U_temp(:) - A2*dU2(:));
                    else
                        dU3(1:end) 	= (A1_inv*U_temp(:) - A2*dU2(:) - A3*dU1(:));
                    end

                    % ---------------- domain decomposition send data around -----------------
                    if op.model.num_model ~=1
                        if op.model.ddcompz~=1
                            if op.model.ddcompz_id ~= 1
                                labSend(dU3(1:op.model.np_extra*2,:,:),labindex-1)
                            end
                            if op.model.ddcompz_id ~= op.model.ddcompz
                                dU3(end-2*op.model.np_extra+1:end,:,:) = dU3(end-2*op.model.np_extra+1:end,:,:) + labReceive(labindex + 1);
                                labSend(dU3(end-2*op.model.np_extra+1:end,:,:),labindex+1)
                            end
                            if op.model.ddcompz_id ~= 1
                                dU3(1:op.model.np_extra*2,:,:) = labReceive(labindex-1);
                            end
                        end
                        
                        if op.model.ddcompx~=1
                            if op.model.ddcompx_id ~= 1
                                labSend(dU3(:,1:op.model.np_extra*2,:),labindex-op.model.ddcompz)
                            end
                            if op.model.ddcompx_id ~= op.model.ddcompx
                                dU3(:,((end-2*op.model.np_extra+1):end),:) = dU3(:,((end-2*op.model.np_extra+1):end),:) + labReceive(labindex + op.model.ddcompz);
                                labSend(dU3(:,(end-2*op.model.np_extra+1:end),:),labindex+op.model.ddcompz)
                            end
                            if op.model.ddcompx_id ~= 1
                                dU3(:,1:2*op.model.np_extra,:) =  labReceive(labindex-op.model.ddcompz);
                            end
                        end
                        
                        
                        if op.model.ddcompy~=1
                            if op.model.ddcompy_id ~= 1
                                labSend(dU3(:,:,1:2*op.model.np_extra),labindex-op.model.ddcompz*op.model.ddcompx)
                            end
                            if op.model.ddcompy_id ~= op.model.ddcompy
                                dU3(:,:,(end-2*op.model.np_extra+1:end)) = dU3(:,:,(end-2*op.model.np_extra+1):end) + labReceive(labindex + op.model.ddcompz*op.model.ddcompx);
                                labSend(dU3(:,:,(end-2*op.model.np_extra+1):end),labindex+op.model.ddcompz*op.model.ddcompx)
                            end
                            if op.model.ddcompy_id ~= 1
                                dU3(:,:,1:2*op.model.np_extra) = labReceive(labindex-op.model.ddcompz*op.model.ddcompx);
                            end
                        end
                        
                    end
                    % --------------- end of domain decomposition -----------------------
                    
                    % Update output
					if ~isempty(Pr)
                   	 	U_shot_record_temp(ti,:) = Pr * dU3(:);
					end
                    %
                    % figure(10);plot(U_shot_record_temp(ti,:));pause(.01);
                    U1	= U2;
                    U2	= U3;
                    dU1	= dU2;
                    dU2	= dU3;
                    
                    
                end
                
                
                y = Pt*U_shot_record_temp(:);
            end
        end
    end
end

function [U]= myReadWf(ti,src,lab)
	name=['t' num2str(ti) 'src' num2str(src) 'lab' num2str(lab)];
	fid=fopen(name,'rb');
	U=fread(fid,inf);
	fclose(fid);
	system(['rm -f ' name]);
end