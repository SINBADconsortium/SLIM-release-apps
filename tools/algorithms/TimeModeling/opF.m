classdef opF < opSpot
    
    % opF -- 3D forward modeling operator.
    %
    % Linear operator that (implicitely) multiplies a source 'q', defined
    % in time domain as a 1D array of size nt and propagates it up
    % to the receivers location. The output of a multiplication operation
    % is thus data, 'd', measured at the level of the receivers in time
    % domain. Transpose multiplication mode will back propagate receiver
    % data 'd' to a source 'q'.
    %
    %
    %
    % INPUTS to the constructor
    %
    %   'model'  : stucture containing the model parameters
    %       'm'  : square slowness [s^2/m^2]
    %     'dens' : density [g/cm^3]
    %     'anis' : anisotropy parameters (epsilon,delta,theta in a
    %     structure)
    % INPUT to multiplication
    %
    %   'x':    direct mode: time domain vectorized source vector of length
    %               nt
    %
    %           adjoint mode: time domain receiver data of length nrec*nt
    %
    %
    % OUTPUT from multiplication
    %
    %   'y':    direct mode [cell array]:
    %              y{1}: time domain receiver data
    %              y{2}: time domain source wave-field (if required by flag
    %                   'RAM'; else = [] and wavefield saved on disk)
    %
    %           adjoint mode [vector]: time domain adjoint source wave
    %               field.
    %
    %
    % Author: Mathias Louboutin from Xianli's code
	%         Seismic Laboratory for Imaging and Modeling
	%         Department of Earth, Ocean, and Atmosperic Sciences
	%         The University of British Columbia
    % Date: December, 2014
    % Edited by Philipp Witte, July 2015
    %
    
    properties
        
        model=[];
        v=[];
        den=[];
        thomsen=[];
        
    end
    
    methods
        %% Constructor
        function op=opF(v,model,dens,thomsen)
            
            nt=length((0:model.dt:model.T)); %Number of time samples
            
            %Number of receivers (size of the operator)
            nrec=length(model.rec_idx);
            if isfield(model,'multi')
                nsrc=model.multi;
            else
                nsrc=1;
            end
            if nrec==0
                nrec=1;
            end
            subT2=length(model.NyqT);
            %Creation of the SPOT operator
            op=op@opSpot('opF',nrec*subT2,nsrc*nt);
            % Model,velocity, density and anisotropy fields
            op.model=model;
            op.v=v;
            op.den=dens;
            op.thomsen = thomsen;
            
        end
    end
    
    
    methods ( Access = protected )
        %% Multiplication function
        function y=multiply(op,x,mode)
            
            %**********************************************************
            % DIRECT MODE (mode=1): forward propagation from the source
            %
            % INPUT: time domain source
            %
            %   Vector of length nt*N, where nt is the number of time
            %   samples and N is the number of grid points on the physical
            %   grid. This source vector can be built vectorizing a 3D
            %   source array of size nz*nx*nt.
            %
            % OUTPUT: time domain receiver data and, if required (flag:
            %   FULLWAVEFIELD), time domain source wave-field
            %
            %**************************************************************
            % ADJOINT MODE (mode=2): backpropagation from the receivers
            %
            % INPUT: time domain receiver data
            %
            %   Vector of size nt*nrec, where nt is the number of time
            %   samples and nrec is the number of receivers.
            %
            % OUTPUT: time domain adjoint source
            %
            %**************************************************************
            
            %Checking of the multiplication mode
            adjoint_mode=false;
            if mode==2,
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
            else
                nsrc=1;
            end

            % Time stepping ,matrices
            
            %%=========================================================================================================
            % 									Time domain forward modeling
            % ==========================================================================================================
            % 	Any wave equation can be writen as a following structure
            %		A1 U_(t-1) + A2 U_t + A3 U_(t+1) = q_t
            % 	Take acoustic wave equation with only velocity for example,
            %		A1,A3: diagnal matrix.   e.g, diag(1./(v^2 *dt^2)), if no PML
            % 		A2: symetric matrix with few off diagnals.  L +  diag(1./(v^2 *dt^2))
            %
            % 	Including PML boundary
            %		A1 = (1+yeta*dt)/(v^2*dt^2);
            %		A2 = L - (1/v^2)(-2/dt^2  + yeta^2);
            %		A3 = (1-yeta*dt)/(v^2*dt^2);
            %
            %		Forward_modeling_mode
            %			fmode: 1, 2 and -1 stands for forward from T1, forward from Tmax and backpergation
            %					from Tmax, respectively.
            % 		need to have A and its adjoint
            %
            % 		4th order in space and 2 order in time
            %
            %
            %	forward modeling solve for A U  = Q
            %
            %		[ A1                           ]  U1    Q1
            %		| A2 A1                        |  U2    Q2
            %		| A3 A2 A1                     |  U3  = Q3
            %		|             . . .            |  U4    Q4
            %		[                      A3 A2 A1]  U5    Q5
            %
            %	adjoint modeling  solve for AT V = dU
            %
            %		[A1T A2T A3T                   ]  V1    dU1
            %		|            . . .             |  V2    dU2
            %		|                   A1T A2T A3T|  V3  = dU3
            %		|                       A1T A2T|  V4    dU4
            %		[                           A1T]  V5    dU5
            %
            % Notice: Adjoint modeling is back propergation, be sure form A1T A2T.
            %		if no PML, A1 and A3 are the same.
            %
            
            
            if adjoint_mode
                %%----------------------- Back propagate --------------------%%
                nr=length(op.model.rlx_loc);
                % Receiver projector
                Pr	= Rec3D(op.model);
                % Time projector from shot record time axis to full
                % time axis
                Pt=kron(speye(nr),opLInterp1D(NyqT,fullT));
                % generate time stepping stencil matrix
                [~,A1_inv,A2,A3,~,~,~,~] = generate_time_stepping_stencil_matrix(op.v,op.model,op.den,op.thomsen);
                % combine some matrix to same time
                Pr = A1_inv'*Pr';
                A2 = A1_inv'*A2' ;
                A3 = A1_inv'*A3';
                
                % allocate wavefield snapshot
                U1	= zeros(length(op.model.mmz_loc),length(op.model.mmx_loc),length(op.model.mmy_loc));
                U2	= zeros(length(op.model.mmz_loc),length(op.model.mmx_loc),length(op.model.mmy_loc));
                U3	= zeros(length(op.model.mmz_loc),length(op.model.mmx_loc),length(op.model.mmy_loc));
                
                % Receiver wavefield
                x=reshape(Pt*x(:),nt,nr);
                
                % Output (at source position only
                y	= zeros(size(op.model.slx_loc,1),nt);
                % Source projector
                Ps	= Src3D(op.model);
                
                % Start loop over time
                for ti=nt:-1:1
                    
                    
                    if ~isempty(Pr)
                        u_adj=x(ti,:)';
                        
                        if ti == nt
                            U3(1:end)   = full((Pr* u_adj));
                        elseif ti == nt-1
                            U3(1:end)   = Pr * u_adj - A2*U2(:);
                        else
                            U3(1:end) 	= Pr * u_adj - A2*U2(:) -  A3*U1(:);
                        end
                    else
                        if ti == nt
                            U3(1:end)   = 0*U3(1:end) ;
                        elseif ti == nt-1
                            U3(1:end)   = -  A2*U2(:);
                        else
                            U3(1:end) 	= -  A2*U2(:) - A3*U1(:);
                        end
                    end
                    
                    %                         if labindex ==1 && mod(ti,50)==0
                    %                              disp([repmat(' ',1,10),'Time=',num2str(ti*op.model.dt),'ms, Total=',num2str((nt-1)*op.model.dt),'ms']);
                    %                         end
                    
                    % ------------ domain decomposition send data around ---------
                    if op.model.num_model ~=1
                        if op.model.ddcompz~=1
                            if op.model.ddcompz_id ~= 1
                                labSend(U3(1:op.model.np_extra*2,:,:),labindex-1)
                            end
                            if op.model.ddcompz_id ~= op.model.ddcompz
                                U3(end-2*op.model.np_extra+1:end,:,:) = U3(end-2*op.model.np_extra+1:end,:,:) + labReceive(labindex + 1);
                                labSend(U3(end-2*op.model.np_extra+1:end,:,:),labindex+1)
                            end
                            if op.model.ddcompz_id ~= 1
                                U3(1:op.model.np_extra*2,:,:) = labReceive(labindex-1);
                            end
                        end
                        
                        if op.model.ddcompx~=1
                            if op.model.ddcompx_id ~= 1
                                labSend(U3(:,1:op.model.np_extra*2,:),labindex-op.model.ddcompz)
                            end
                            if op.model.ddcompx_id ~= op.model.ddcompx
                                U3(:,((end-2*op.model.np_extra+1):end),:) = U3(:,((end-2*op.model.np_extra+1):end),:) + labReceive(labindex + op.model.ddcompz);
                                labSend(U3(:,(end-2*op.model.np_extra+1:end),:),labindex+op.model.ddcompz)
                            end
                            if op.model.ddcompx_id ~= 1
                                U3(:,1:2*op.model.np_extra,:) =  labReceive(labindex-op.model.ddcompz);
                            end
                        end
                        
                        
                        if op.model.ddcompy~=1
                            if op.model.ddcompy_id ~= 1
                                labSend(U3(:,:,1:2*op.model.np_extra),labindex-op.model.ddcompz*op.model.ddcompx)
                            end
                            if op.model.ddcompy_id ~= op.model.ddcompy
                                U3(:,:,(end-2*op.model.np_extra+1:end)) = U3(:,:,(end-2*op.model.np_extra+1):end) + labReceive(labindex + op.model.ddcompz*op.model.ddcompx);
                                labSend(U3(:,:,(end-2*op.model.np_extra+1):end),labindex+op.model.ddcompz*op.model.ddcompx)
                            end
                            if op.model.ddcompy_id ~= 1
                                U3(:,:,1:2*op.model.np_extra) = labReceive(labindex-op.model.ddcompz*op.model.ddcompx);
                            end
                        end
                        
                    end
                    % --------------- end of domain decomposition -----------------------
                    
                    
                    % output source wavefield
                    if ~isempty(Ps)
                        y(ti) = Ps * U3(:);
                    end
                    % swap wavefields for next time step
                    U1 = U2;
                    U2 = U3;
                    
                    % Free surface, mirror the wavefield above see
                    % level
                    if op.model.freesurface
                        for loc=1:op.model.space_order/2
                            U3(loc,:,:)=-U3(op.model.space_order-loc+2,:,:);
                        end
                        U3(op.model.space_order/2+1,:,:)=0*U3(op.model.space_order/2+1,:,:);
                    end
                    
                    % End of time loop
                end
                
            else
                %%----------------------- Forward propagate --------------------%%
                
                % Number of receivers and receiver projector
                nrecloc=length(op.model.rlx_loc);
                Pr	= Rec3D(op.model);
                
                % generate time stepping stencil matrix
                [~,A1_inv,A2,A3,~,~,~,~] = generate_time_stepping_stencil_matrix(op.v,op.model,op.den,op.thomsen);
                % combine some matrix to save time
                A2 = A1_inv * A2;
                A3 = A1_inv * A3;
                % Seup in case of continuous acquisition
                if isfield(op.model,'tfire')
                    disp('Continuous acquisition');
                    % Define a smaller shot record saved on disk every 10s
                    % of simulation
                    if nt>10000
                        chunck_size=10000;
                        recloc=zeros(chunck_size,nrecloc);
                        Tchunck=0:op.model.dt:(chunck_size-1)*op.model.dt;
                        NyqChunck=0:4:(chunck_size-1)*op.model.dt;
                        indChunk=1;
                        numChunk=1;
                        
                        % Time projector from full time axis to shot record
                        % time axis
                        Pt=kron(speye(nrecloc),opLInterp1D(Tchunck,NyqChunck));
                    else
                        Pt=kron(speye(nrecloc),opLInterp1D(fullT,NyqT));
                        recloc=zeros(nt,nrecloc);
                        indChunk=1;
                        numChunk=1;
                        chunck_size=nt;
                    end
                    k = 1;
                    % Source firing time indexes
                    tfire_ind = floor(op.model.tfire/op.model.dt);
                    tisrc = 0;
                    if tfire_ind(k)==0
                        tisrc = 1;
                        op.model.slx = op.model.sfire(k);
                        op.model.slx_loc = op.model.sfire(k);
                        Ps = Src3D(op.model);
                        Ps = A1_inv*Ps';
                        disp(['shot Idx = ',num2str(k)]);
                        k = k+1;
                    else
                        Ps=[];
                    end
                else
                    % Setup for conventional acquisition
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
                    % clera unnecessary matrix
                    if ~isempty(Ps)
                        Ps = A1_inv * Ps';
                    end
                    % Local shot record
                    recloc=zeros(nt,nrecloc);
                    clear A1_inv;
                end
                % allocate wavefield snapshot
                U1	= zeros(length(op.model.mmz_loc),length(op.model.mmx_loc),length(op.model.mmy_loc));
                U2	= zeros(length(op.model.mmz_loc),length(op.model.mmx_loc),length(op.model.mmy_loc));
                U3	= zeros(length(op.model.mmz_loc),length(op.model.mmx_loc),length(op.model.mmy_loc));
                % Wavefield for gradient (no wavefield for modelling
                % only)
                if ~strcmp(op.model.save,'disk')
                    wfloc=zeros(length(op.model.nsub),numel(U3));
                end
                
                % Reshape source for time stepping
                x=reshape(x,nsrc,nt);
                
                % Start loop over time
                for ti=1:nt
                    
                    
                    
                    % Setup source in case of continuous acquisition
                    if isfield(op.model,'tfire')
                        if ti <= tfire_ind(end)
                            if ti < tfire_ind(1)
                                Ps = [];
                            elseif ti == tfire_ind(k)
                                tisrc = 1;
                                op.model.slx = op.model.sfire(k);
                                op.model.slx_loc = op.model.sfire(k);
                                Ps = Src3D(op.model);
                                Ps = A1_inv*Ps';
                                disp(['shot Idx = ',num2str(k)]);
                                k = k+1;
                            end
                        end
                    else
                        tisrc=ti;
                    end
                    
                    if isempty(Ps)
                        if ti == 1
                            U3(1:end)   = 0;
                        elseif ti == 2
                            U3(1:end)   =  - A2 * U2(:);
                        else
                            U3(1:end) 	= - A2 * U2(:) - A3 * U1(:);
                        end
                    else
                        if ti == 1
                            U3(1:end)   = full(Ps * x(:,tisrc));
                        elseif ti == 2
                            U3(1:end)   = Ps * x(:,tisrc) - A2 * U2(:);
                        else
                            U3(1:end) 	= Ps * x(:,tisrc) - A2 * U2(:) - A3 * U1(:);
                        end
                    end
					
                    tisrc=tisrc+1;
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
%                     if labindex==1 && mod(ti,50)==0
%                         disp([repmat(' ',1,10),'Time=',num2str(ti*op.model.dt),'ms, Total=',num2str((nt-1)*op.model.dt),'ms']);
%                         figure(1);imagesc(U3(:,:,1),[-1 1]);title([num2str(ti)]);colormap(jet);colorbar;pause(.001);
%                     end


                    % ---------------- Save wavefield for gradient--------%

                    if strcmp(op.model.save,'disk')
                        if(length(find(op.model.nsub==ti))>0)
							mySaveWf(ti,op.model.numsave,labindex,U3);
                        end
                    else
                        if(length(find(op.model.nsub==ti))>0)
                            wfloc(find(op.model.nsub==ti),:)=U3(:);
                        end
                    end

                    %--------------- Save receiver data -----------------%
                    if isfield(op.model,'tfire')
                        if ~isempty(Pr)
                            recloc(indChunk,:)=Pr*U3(:);
                            indChunk=indChunk+1;
                        end
                        if (indChunk==chunck_size+1 || ti==nt) && nt>10000
                            indChunk=1;
                            recloc=reshape(Pt*recloc(:),length(NyqChunck),nrecloc);
                            name=['shot_rec_' num2str(numChunk) '.mat'];
                            save(name, '-v7.3', 'recloc','chunck_size');
                            numChunk=numChunk+1;
                            recloc=zeros(chunck_size,nrecloc);
                        end
                    else
                        %Receiver data
                        if ~isempty(Pr)
                            recloc(ti,:)=Pr*U3(:);
                        end
                    end
                    
                    % Free surface, mirror the wavefield above see
                    % level
                    if op.model.freesurface
                        for loc=1:op.model.space_order/2
                            U3(loc,:,:)=-U3(op.model.space_order-loc+2,:,:);
                        end
                        U3(op.model.space_order/2+1,:,:)=0*U3(op.model.space_order/2+1,:,:);
                    end
                    
                    % Swap wavefields
                    U1 = U2;
                    U2 = U3;
                    
                    % End of time loop
                end
                
                % Resize output correctly
                if ~isfield(op.model,'tfire') || nt<10001
                    recloc=Pt*recloc(:);
                else
                    recloc=[];
                end
                if strcmp(op.model.save,'disk')
                    y=recloc;
                else
                    y{1}=recloc;
                    y{2}=wfloc;
                end
                
            end
            
        end
    end
end

function []= mySaveWf(ti,src,lab,U)
	name=['t' num2str(ti) 'src' num2str(src) 'lab' num2str(lab)];
	fid=fopen(name,'wb');
	fwrite(fid,U);
	fclose(fid);
end