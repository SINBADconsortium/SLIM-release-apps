function [output_data,wavefield_chk_point,loc_rec_idx] = XiangLi_Time(model,src,fwdpara,wavefield_chk_point,input_data);
%% [output_data,J] = XiangLi_Time(vel,den,src,fwdpara);
% 2D/3D time domain forward modeling, 
% output arguements
%	output_data: distributed recording data, formed into a long vector, size of which is Ns x Nr x Nt.
%				 If use domain decomposition or multiple shots for one core. The index of receiver will be not well orgnized.
%	wavefield_chk_point: wavefield checkpoints for backward modeling
%	loc_rec_idx: If use domain decomposition or multiple shots for one core. this is the index of receivers.
%				
% input arguements
%	model				: structure of model parameters
%			vel					: vectorized velocity model m/s, 
%			den					: vectorized density mdoel g/cm^3
%%		TTI anisotrophy setting
%			epsilon				: Thomsen's parameter
%			delta				: Thomsen's parameter
%			theta				: angle of the symmetry axis with respect to the x and z-axis 
%			phi					: angle of the symmetry axis with respect to the y and z-axis (=0 for 2d)
%
%		-----------------------------------------------------------------------------------
%		| NOTICE:     vel should already include boundaries.                              |
%		| always need 'model.vel', other parameters depends on which wave equation to use,| 
%		| check fwdpara.wave_equ for more information.                                    |
%		-----------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------
%	src					: SEE 'fwdpara' setting (source location parameters);
% ---------------------------------------------------------------------------------------------------------
%	wavefield_chk_point	: wavefield checkpoint, To use this you need to do the modeling first 
%						  (fwdpara.fmode=1), and save checkpoint. Then you can reconstruct wavefeild
%						  from Tmax to T1.
% ---------------------------------------------------------------------------------------------------------
%	input_data			: input data for jacobian of FWI objective function;
%						  for forward mode (fwdpara.fmode= 3), input data is model perturbation.
%						  for adjoint mode (fwdpara.fmode=-3), input data is residual wavefield.
% ---------------------------------------------------------------------------------------------------------
%	fwdpara				: structure of modeling parameters
%		log out put:
%			fid: output log. =1, show infomation in command window
%							 =fopen('filename.txt','w'). output to a 'txt' file
%							 =0, silence  
%		wave equations
%			wave_equ:	1, acoustic with only velocity
%						2, acoustic with variable density 
%						3, TTI acoustic without density
%		finite difference
%			space_order: 4, 6, 8, 16. 4th, 6th, 8th order in space (high order for density in progess)
%		model geometry parameters, 
%			mmx: row vector (1 x n), model coordinate along x. (unit: meters)
%			mmy: row vector (1 x n), model coordinate along y. (unit: meters); length(mmy)=1 
%				 stands for 2D model.
%			mmz: row vector (1 x n), model coordinate along z. (unit: meters)
%		source location parameters
%			stype:	'seq_same', 'src' needs to be a vector (source wavelet,length = Nt),
%								all shot has same source wavelet.
%					'seq_diff', 'src' needs to be a matrix (size of it = Ns x Nt, Ns=No of 
%														  total shot locations);
%								each columne is a difference source wavelet
%					'sim_same', 'src' needs to be a matrix (size of it = Ns x Nt )
%								each slice along y is one sim. source.
%					'sim_diff', 'src' needs to be a 3D cube (size of it = Ns x Nt x Nshot )
%								each slice along y is one sim. source. each shot has
%								different source locations (like back propergation of 
%								marine acquisition). 
%								BUT EACH SHOT NEEDS TO HAVE SAME # OF SOURCE LOCATIONS
%								(put 0s if they are not).
%			slz: source location coordinate of z. (unit: meters)
%			slx: source location coordinate of x. (unit: meters)
%			sly: source location coordinate of y. (unit: meters); length=1 stands for 2D model.
%					For 'seq_same', 'seq_diff', they NEED TO BE transposed vectors (1xN), each 
%					element is a source point. 
%					For 'sim_same', they are vectors. 
%					For 'sim_diff', they are matrix, each column of them is source location
%					coordinate for one sim. shot.
%			shot_id: shot id for the first shot, e.g, 50, the first shot saved as shot_50.mat 
%					(default,=1).
%			saved: how to output simulated shot record. 
%					fwdpara.saved=0; output shot record as Composite format
%							NOTICE: f use domain decomp, this matrix will be NOT arranged shot by 
%									shot, if arrange data worker by worker. be sure do not mess up 
%									receiver index.
%					fwdpara.saved=1, output shot record as distributed matrix, [Nr*Ns, Nt]
%							NOTICE: if use domain decomp, this matrix will be NOT arranged shot by 
%									shot, if arrange data worker by worker. be sure do not mess up 
%									receiver index.
%					fwdpara.saved='filename', file name, save to to disk as ('./filename/shot_1.mat', ... ,'./filename/shot_n.mat');
%		receiver location parameter
%			rtype: 'full' or 'marine' stands for full acquisition and marine acq., repectively
%			rlz: receiver location coordinate of z. (unit: meters)
%			rlx: receiver location coordinate of x. (unit: meters)
%			rly: receiver location coordinate of y. (unit: meters); length=1 stands for 2D model.
%					For 'full', they are vectors. all shots share same receivers
%					For 'marine', they are matrixs. each column is receiver location for one
%					single shots. NOTE: SIZE(rl*,2) = NS;		
%		time axis & central frequency of source
%			taxis: time axis for wavefield and src. Ensure Nt of 'src' = length(taxis);
%			dt   : time step length, Ensure time step length satisfy the stability condition.
%			fcent: central frequency of source wavelet. for stability checking.
%		Boundary (unit: meters)
%			abx: number of grids for boundary along x.
%			aby: number of grids for boundary along y.length=1 stands for 2D model.
%			abz: number of grids for boundary along z.
%		freesurface
%			frees: 1 and 0 stands for freesurface and NO freesurface, respectively.
%		Modeling_mode
%			fmode:
%				 1 : solve forward modeling from T1 to Tmax.
%				 2 : solve forward modeling from Tmax to T1. You need to have U_tmax and U_tmax-1 and maybe checkpoints.
%				-1 : solve adjoint wavefield from Tmax to T1.
%				 3 : linear born modeling.              (forward mode of jacobian for FWI objective function)
%				-3 : gradient of FWIobjective function. (adjoint mode of jacobian for FWI objective function)
%		domain decomposition setup
%			ddcompx: split model along x dimension
%			ddcompz: split model along z dimension
%			ddcompy: split model along y dimension
%		choose model range. Sometimes each shot only sees a small part of the whole model. 
%			cut_model: =0, alway uses the whole model; =1, for each shot, uses a subset of the model.
%			cut_rest : []
%		tricks of computing gradient
%			chkp_space: space of saving instantaneous wavefield, e.g, fwdpara.chkp_space = 10, if you 
%						have Nt time step,You will save 2 instantaneous wavefields every 10 wavefields,
%						fwdpara.chkp_space=1 means save all wavefield. 
%			chkp_save: 	chkp_save=1, put all checkpoint wavefield in to a big Composite data cube. 
%						chkp_save=2, save all checkpoint wavefields of one local shot to disk. shot by shot
%						chkp_save=3, save each instantaneous checkpoint wavefield of each shot to disk.
%						chkp_save=0, do not save checkpoint; in this case, you can NOT compute gradient.
%									 Use it, if you only want to simulate some shots (e.g, observed data). 
%		ploting wavefield snapshot
%			disp_snapshot: =1, display wavefeild snapshot. 0, silent. 
%							NOTICE: you can not display while using domain decomposition.
%		
%		----------------------------------------------------------------------------
%         parameters for Jacobian and its adjoint of FWI least-squares objective
%		----------------------------------------------------------------------------
%		Unit of updates
%			v_up_type: 'slowness' or 'velocity'
% ---------------------------------------------------------------------------------------------------------
%
%
%   put a example here;
% 
%
%
% ===========================================================================================================
%	Author information:
%			Xiang Li, Phd candidate.
%			Seismic Laboratory of Imaging and Modeling
%			Earth and Ocean Science, University of British Columbia
%			Tel: +17789917197
%			Website: https://slim.gatech.edu
%			Email: JW.Li.Xiang@gmail.com
%      
% 	Date: July, 2014
%
% DO NOT Modify this code, you don't know how I coded it up. A simple change will cause crucial mistake.
% I am not responsiable for any problem that caused by using this code. You may not use it for commercial.
% If you need technique support or have any question, concern, sugguestions, please feel free to contact me.
% ===========================================================================================================

%% Note for intermediate structure. 
%
%   fwdpara_dis:  this is a structure that contains all setting for each worker. 
%                      x_id_loc: index of coordinate along x of local model of mmx
%                      y_id_loc: index of coordinate along y of local model of mmy
%                      z_id_loc: index of coordinate along z of local model of mmz
%                         u_idx: index of modeling area, for worker that not on boundary, it should start from 3 to end-2
%                         u_idy: index of modeling area, for worker that not on boundary, it should start from 3 to end-2
%                         u_idz: index of modeling area, for worker that not on boundary, it should start from 3 to end-2
%                      np_extra: depends on finite difference scheme. 2 point for 4th order.
%                       mmx_loc: coordinate along x for local model split
%                       mmy_loc: coordinate along y for local model split
%                       mmz_loc: coordinate along z for local model split
%                    ddcompx_id: model split index along x
%                    ddcompy_id: model split index along y
%                    ddcompz_id: model split index along z
%                      slx_band: source band locations along x, for local model splits.
%                      sly_band: source band locations along y, for local model splits.
%                      slz_band: source band locations along z, for local model splits.
%                shots_parition: source split partition. 
% max_number_of_shot_all_worker: num of all source
%                           rlx: receiver location along x for local model split band
%                           rly: receiver location along y for local model split band
%                           rlz: receiver location along z for local model split band
%                       rlx_loc: receiver location along x that on the local grid
%                       rly_loc: receiver location along y that on the local grid
%                       rlz_loc: receiver location along z that on the local grid
%                       rec_idx: index of local receiver to global receiver position
%                           slx: source location along x for local model split band
%                           sly: source location along y for local model split band
%                           slz: source location along z for local model split band
%                       slx_loc: source location along x that on the local grid
%                       sly_loc: source location along y that on the local grid
%                       slz_loc: source location along z that on the local grid
%                       src_idx: index of local source to global source position
%
% -----------------------------------------------------------------------------------------------------------
	% global vel_dis den_dis dvel_dis fwdpara_dis src_dis

	if nargin < 4, wavefield_chk_point	= 0;end
	if nargin < 5, input_data			= 0;end

	fwdpara = check_fwdpara_struct(fwdpara,model);

	% ------------------------------------------ check modeling parameters  -------------------------------------
	check_stability_condition(model,fwdpara);

	% ------------------------------------------ domain decomposition setup -------------------------------------
	[fwdpara_dis] = setup_for_domain_decomp(model,fwdpara,input_data);
	clear model
	
	% --------------------------------------------- do forward mdoeling -----------------------------------------
	[output_data,fwdpara_dis,wavefield_chk_point,loc_rec_idx] = do_Forward_modeling(fwdpara_dis,src,fwdpara,wavefield_chk_point,input_data);

end



%% do foward modeling
function  [output_data,fwdpara_dis,wavefield_chk_point,loc_rec_idx] = do_Forward_modeling(fwdpara_dis,src,fwdpara,wavefield_chk_point,input_data)
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
	%%
	% global vel_dis den_dis dvel_dis fwdpara_dis

	
	if fwdpara.fmode == 1;
		% ------------------------------------------------------
		% do forward wavefield for T1 to Tmax
		% ------------------------------------------------------
		if fwdpara.fid && labindex == 1
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n'])
			fprintf(fwdpara.fid,[repmat(' ',1,10),'Solving forward modeling from T1 to Tmax\n'])
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n'])
		end
	

		% check point wavfield
		wavefield_chk_point = 0;
		
		spmd		
			%
			ntrace_loc 			= 0; % number of loc traces, Nr_loc * Ns_loc; 
			Num_of_shots_loc	= size(fwdpara_dis.slx_band,2);
			Num_of_time_step	= length(fwdpara.taxis);
			loc_rec_idx			= [];
			
			
			% locate receiver position (all shot has same receiver)
			if strcmp(fwdpara.rtype,'full')
				fwdpara_dis.rlx	= fwdpara.rlx;
				fwdpara_dis.rly	= fwdpara.rly;
				fwdpara_dis.rlz	= fwdpara.rlz;
				fwdpara_dis		= locate_rcv_position(fwdpara_dis);
				Pr				= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
												fwdpara_dis.rlz_loc,fwdpara_dis.rlx_loc,fwdpara_dis.rly_loc);
				rec_idx			= fwdpara_dis.rec_idx;
			elseif strcmp(fwdpara.rtype,'marine')
				rec_idx			= cell(Num_of_shots_loc,1);
			end

			% Num_of_rec_global	= size(fwdpara_dis.rlx,1);
				
			% create place for shot record
			output_data		= cell(1,Num_of_shots_loc);
			
			% create place for checkpoint wavefield 
			if fwdpara.chkp_save == 1
				wavefield_chk_point	= cell(Num_of_shots_loc,1);
			end
			[chk_point_idx,N_chkp]	= generate_check_point_idx(Num_of_time_step,fwdpara.chkp_space);
		

			% generate time stepping stencil matrix
			% A1 = [];A1_inv = []; A2 = [];A3 = [];A3_inv = [];	
			[~,A1_inv,A2,A3,~] = generate_time_stepping_stencil_matrix(fwdpara_dis,fwdpara);
			% [A3,~,A2,~,A1_inv] = generate_time_stepping_stencil_matrix(fwdpara_dis,fwdpara);
			% combine some matrix to save time
			A2 = A1_inv * A2; 
			A3 = A1_inv * A3; 

			
			% ---------------------------------- loop over shots ---------------------------------------
			for si  = 1:Num_of_shots_loc
				
				if fwdpara.fid && labindex == 1
					fprintf(fwdpara.fid,[repmat(' ',1,5),repmat('-',1,5),' Simulating shot ',num2str(si),' ',repmat('-',1,5),'\n']);
				end
				
				% locate src position and receiver position (each shot has diff receivers) 
				fwdpara_dis.slx		= vec(fwdpara_dis.slx_band(:,si));
				fwdpara_dis.sly		= vec(fwdpara_dis.sly_band(:,si));
				fwdpara_dis.slz		= vec(fwdpara_dis.slz_band(:,si));
				fwdpara_dis 		= locate_shot_position(fwdpara_dis);
				Ps					= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
												fwdpara_dis.slz_loc,fwdpara_dis.slx_loc,fwdpara_dis.sly_loc);
												
				if strcmp(fwdpara.rtype,'marine')
					fwdpara_dis.rlx	= fwdpara.rlx(:,fwdpara_dis.shots_id(si));
					fwdpara_dis.rly	= fwdpara.rly(:,fwdpara_dis.shots_id(si));
					fwdpara_dis.rlz	= fwdpara.rlz(:,fwdpara_dis.shots_id(si));
					fwdpara_dis 	= locate_rcv_position(fwdpara_dis);
					Pr				= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
													fwdpara_dis.rlz_loc,fwdpara_dis.rlx_loc,fwdpara_dis.rly_loc);
				end
				
				loc_rec_idx = [loc_rec_idx(:);fwdpara_dis.rec_idx + size(fwdpara.rlx,1).*(fwdpara_dis.shots_id(si) - 1)];
				
				% generate source wavelet
				Q  					= prepare_source_wavelet(src,fwdpara,fwdpara_dis,si);
	
				% allocate place for loc shot record
				U_shot_record_temp	= zeros(length(fwdpara_dis.rec_idx),Num_of_time_step);

				
				% combine some matrix to save time
				Ps = A1_inv * Ps';
	
				% allocate wavefield snapshot
				U1	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				U2	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				U3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				
				
				% allocate space for checkpoint
				if fwdpara.chkp_save == 1 || fwdpara.chkp_save == 2
					chk_point_wavefeild_temp	= zeros(numel(U3),N_chkp);
				end
				chk_point_count				= 1;			


				% ---------------------------- loop over time step ---------------------------------
				for ti = 1:Num_of_time_step

					%  forward modeling solve for A U  = Q
					%
					% [ A1                           ]  U1    Q1
					% | A2 A1                        |  U2    Q2
					% | A3 A2 A1                     |  U3  = Q3
					% |             . . .            |  U4    Q4
					% [                      A3 A2 A1]  U5    Q5
					if ti == 1
						U3(1:end)   = full(Ps * Q(:,ti));
					elseif ti == 2
						U3(1:end)   = Ps * Q(:,ti) - A2 * U2(:);
					else
						U3(1:end) 	= Ps * Q(:,ti) - A2 * U2(:) - A3 * U1(:);
					end
					% U3(:) = A1_inv * U3(:);
				
					
					
					if labindex ==1 && mod(ti,50)==0 && fwdpara.fid
						fprintf(fwdpara.fid,[repmat(' ',1,10),'Time=',num2str(fwdpara.taxis(ti)),'s, Total=',num2str(fwdpara.taxis(end)),'s.\n'])
	
					end
				

					% ---------------- domain decomposition send data around -----------------
					if fwdpara_dis.num_model ~=1
						if fwdpara.ddcompz~=1
							if fwdpara_dis.ddcompz_id ~= 1
								labSend(U3(fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra,:,:),labindex-1)
							end
							if fwdpara_dis.ddcompz_id ~= fwdpara.ddcompz
								U3((end-fwdpara_dis.np_extra+1:end),:,:) = labReceive(labindex + 1);
								labSend(U3((end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra),:,:),labindex+1)
							end
							if fwdpara_dis.ddcompz_id ~= 1
								U3(1:fwdpara_dis.np_extra,:,:) = labReceive(labindex-1);
							end
						end
						labBarrier
						if fwdpara.ddcompx~=1
							if fwdpara_dis.ddcompx_id ~= 1
								labSend(U3(:,fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra,:),labindex-fwdpara.ddcompz)
							end
							if fwdpara_dis.ddcompx_id ~= fwdpara.ddcompx
								U3(:,(end-fwdpara_dis.np_extra+1:end),:) = labReceive(labindex + fwdpara.ddcompz);
								labSend(U3(:,(end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra),:),labindex+fwdpara.ddcompz)
							end
							if fwdpara_dis.ddcompx_id ~= 1
								U3(:,1:fwdpara_dis.np_extra,:) =  labReceive(labindex-fwdpara.ddcompz);
							end
						end
						labBarrier

						if fwdpara.ddcompy~=1
							if fwdpara_dis.ddcompy_id ~= 1
								labSend(U3(:,:,fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra),labindex-fwdpara.ddcompz*fwdpara.ddcompx)
							end
							if fwdpara_dis.ddcompy_id ~= fwdpara.ddcompy
								U3(:,:,(end-fwdpara_dis.np_extra+1:end)) = labReceive(labindex + fwdpara.ddcompz*fwdpara.ddcompx);
								labSend(U3(:,:,(end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra)),labindex+fwdpara.ddcompz*fwdpara.ddcompx)
							end
							if fwdpara_dis.ddcompy_id ~= 1
								U3(:,:,1:fwdpara_dis.np_extra) = labReceive(labindex-fwdpara.ddcompz*fwdpara.ddcompx);
							end
						end
						labBarrier
					end										
					% --------------- end of domain decomposition -----------------------


					 

					
 					% ---------------------------- display snapshot----------------------------
 					if fwdpara.disp_snapshot == 1 & mod(ti,20)==0
 						if length(fwdpara.mmy) ==1
 							% 2D ploting
 							figure(1);imagesc(U3);title([num2str(ti)]);caxis([-1 1]);colorbar;pause(.01);

 						else
 							% 3D ploting 
 							figure(1);sliceview(reshape(U3,length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc)));
 							title([num2str(ti)]);pause(.01)
 						end
 					end
					
					
 					% -------------- save check point per time step ---------------------
 					if chk_point_idx(ti) == 1
 						if fwdpara.chkp_save == 1
 							chk_point_wavefeild_temp(:,chk_point_count)	= U3(:);
 							chk_point_count 	= chk_point_count + 1;
 						elseif fwdpara.chkp_save == 2
 							chk_point_wavefeild_temp(:,chk_point_count)	= U3(:);
 							chk_point_count 	= chk_point_count + 1;
 						elseif fwdpara.chkp_save == 3
 							save_chkpoint_to_disk(U3,labindex,si,fwdpara.chkp_save,chk_point_count);
 							chk_point_count 	= chk_point_count + 1;
 						end
 					end
					

					
					% output receiver wavefield
					U_shot_record_temp(:,ti) = Pr * U3(:);
	
					
					U1 = U2;
					U2 = U3;
					U3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				
				end
				% save the wavefield at the last time step
				% U_shot_record_temp(:,ti)	= Pr * U2(:);
				output_data{si}				= U_shot_record_temp.';
				if strcmp(fwdpara.rtype,'marine')
					rec_idx{si}					= fwdpara_dis.rec_idx;
				end
				ntrace_loc = ntrace_loc + length(fwdpara_dis.rec_idx);
				
				
				% -------------- save check point per shot---------------------
				if fwdpara.chkp_save == 1
					wavefield_chk_point{si,1} = chk_point_wavefeild_temp;
				elseif fwdpara.chkp_save == 2
					
					save_chkpoint_to_disk(chk_point_wavefeild_temp,labindex,si,2,chk_point_count)
				end

			end
			
			% different workers may simulated different number of shots (e.g, worker 1 simulate 15 shots, worker 2 simulate 16 shots),
			% 'labBarrier' can not work while other workers are idle. 
			% this make sure all worker finish simulation. You can arrange workers to save one shot time. 
			if Num_of_shots_loc < fwdpara_dis.max_number_of_shot_all_worker
				for si = Num_of_shots_loc+1: fwdpara_dis.max_number_of_shot_all_worker
					for ti = 1:Num_of_time_step 
						labBarrier
						labBarrier
						labBarrier
					end
				end
			end
			
			% put all data together
			codistr			= codistributor1d(2,codistributor1d.unsetPartition,[1,numlabs]);
			ntrace_loc		= codistributed.build(ntrace_loc,codistr,'noCommunication'); 
			
			
		end
		

		[output_data,loc_rec_idx] = output_shot_record(output_data,rec_idx,fwdpara,fwdpara_dis,ntrace_loc,loc_rec_idx);
		
		
	elseif fwdpara.fmode == 2;
		if fwdpara.fid	&& labindex == 1	
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n'])
			fprintf(fwdpara.fid,[repmat(' ',1,10),'Solving forward modeling from Tmax to T1.\n'])
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n'])
		end
%				 2 : solve forward modeling from Tmax to T1. You need to have U_tmax and U_tmax-1 and maybe checkpoints.
%				-1 : solve adjoint wavefield from Tmax to T1.
%				 3 : linear born modeling.              (forward mode of jacobian for FWI objective function)
%				-3 : gradient of FWIobjective function. (adjoint mode of jacobian for FWI objective function)
		
		
		
	% ------------------------------------------------------
	% do forward wavefield from Tmax to T1
	% ------------------------------------------------------
	% NOTE: at least need u_n and u_n-1 and check point.
		spmd		
			%
			
			ntrace_loc 			= 0; % number of loc traces, Nr_loc * Ns_loc; 
			Num_of_shots_loc	= size(fwdpara_dis.slx_band,2);
			Num_of_time_step	= length(fwdpara.taxis);
			loc_rec_idx 		= 0;
			
			% locate receiver position (all shot has same receiver)
			if strcmp(fwdpara.rtype,'full')
				fwdpara_dis.rlx	= fwdpara.rlx;
				fwdpara_dis.rly	= fwdpara.rly;
				fwdpara_dis.rlz	= fwdpara.rlz;
				fwdpara_dis		= locate_rcv_position(fwdpara_dis);
				Pr				= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
												fwdpara_dis.rlz_loc,fwdpara_dis.rlx_loc,fwdpara_dis.rly_loc);
				rec_idx			= fwdpara_dis.rec_idx;
			elseif strcmp(fwdpara.rtype,'marine')
				rec_idx			= cell(Num_of_shots_loc,1);
			end

			% Num_of_rec_global	= size(fwdpara_dis.rlx,1);
		
			% create place for shot record
			output_data		= cell(1,Num_of_shots_loc);
		
			% create place for checkpoint wavefield 
			% if fwdpara.chkp_save == 1
			% 	wavefield_chk_point	= cell(Num_of_shots_loc,1);
			% end
			
			[chk_point_idx,N_chkp]	= generate_check_point_idx(Num_of_time_step,fwdpara.chkp_space);
		
			
			% generate time stepping stencil matrix
			% A1 = [];A1_inv = []; A2 = [];A3 = [];A3_inv = [];	
			[A1,~,A2,~,A3_inv] = generate_time_stepping_stencil_matrix(fwdpara_dis,fwdpara);
			% % combine some matrix to save time
			% Ps = A3_inv * Ps';

			% combine some matrix to save time
			A2	= A3_inv * A2;
			A1	= A3_inv * A1;
			
			% ---------------------------- loop over shots ---------------------------------
			for si  = 1:Num_of_shots_loc
				if fwdpara.fid && labindex == 1
					fprintf(fwdpara.fid,[repmat(' ',1,5),repmat('-',1,5),' Simulating shot ',num2str(si),' ',repmat('-',1,5),'\n']);
				end
				% locate src position and receiver position (each shot has diff receivers) 
				fwdpara_dis.slx		= vec(fwdpara_dis.slx_band(:,si));
				fwdpara_dis.sly		= vec(fwdpara_dis.sly_band(:,si));
				fwdpara_dis.slz		= vec(fwdpara_dis.slz_band(:,si));
				fwdpara_dis 		= locate_shot_position(fwdpara_dis);
				Ps					= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
												fwdpara_dis.slz_loc,fwdpara_dis.slx_loc,fwdpara_dis.sly_loc);
			
				if strcmp(fwdpara.rtype,'marine')
					fwdpara_dis.rlx = fwdpara.rlx(:,fwdpara_dis.shots_id(si));
					fwdpara_dis.rly = fwdpara.rly(:,fwdpara_dis.shots_id(si));
					fwdpara_dis.rlz = fwdpara.rlz(:,fwdpara_dis.shots_id(si));
					fwdpara_dis 	= locate_rcv_position(fwdpara_dis);
					Pr				= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
													fwdpara_dis.rlz_loc,fwdpara_dis.rlx_loc,fwdpara_dis.rly_loc);
				end
				
				loc_rec_idx = [loc_rec_idx(:);fwdpara_dis.rec_idx + size(fwdpara.rlx,1).*(fwdpara_dis.shots_id(si) - 1)];
				
				% generate source wavelet
				
				Q  = prepare_source_wavelet(src,fwdpara,fwdpara_dis,si);
		
				% allocate place for loc shot record
				U_shot_record_temp	= zeros(length(fwdpara_dis.rec_idx),Num_of_time_step);

				% combine some matrix to save time
				Ps = A3_inv * Ps';

				
				
			

				% allocate wavefield snapshot
				U1	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				U2	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				U3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
			
			
				% allocate space for checkpoint
				if fwdpara.chkp_save == 1 
					chk_point_wavefeild_temp = wavefield_chk_point{si,1};
				elseif fwdpara.chkp_save == 2
					chk_point_wavefeild_temp = load_chkpoint_from_disk(labindex,si,fwdpara.chkp_save);
				end
				chk_point_count				 = N_chkp;			
		
		
				% ---------------------------- loop over time step ---------------------------------
				for ti = Num_of_time_step:-1:1

					%  forward modeling solve for A U  = Q
					%
					% [ A1                           ]  U1    Q1
					% | A2 A1                        |  U2    Q2
					% | A3 A2 A1                     |  U3  = Q3
					% |             . . .            |  U4    Q4
					% [                      A3 A2 A1]  U5    Q5
					
					
					% -------------- load check point per time step ---------------------
					if chk_point_idx(ti) == 1
						if fwdpara.chkp_save == 1
							U3(1:end) 				= reshape(chk_point_wavefeild_temp(:,chk_point_count),size(U2));
							chk_point_count = chk_point_count - 1;
						elseif fwdpara.chkp_save == 2
							U3(1:end)  				= reshape(chk_point_wavefeild_temp(:,chk_point_count),size(U2));
							chk_point_count 		= chk_point_count - 1;
						elseif fwdpara.chkp_save == 3
							U3(1:end) 			    = load_chkpoint_from_disk(labindex,si,fwdpara.chkp_save,chk_point_count);
							chk_point_count			= chk_point_count - 1;
						end
						
					else
						U3(1:end)	=	Ps*Q(:,ti+2) - A2*U2(:) - A1 * U1(:);	
						% U3(1:end)	=	A3_inv * U3(:);
						% ---------------- domain decomposition send data around -----------------
						if fwdpara_dis.num_model ~=1
							if fwdpara.ddcompz~=1
								if fwdpara_dis.ddcompz_id ~= 1
									labSend(U3(fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra,:,:),labindex-1)
								end
								if fwdpara_dis.ddcompz_id ~= fwdpara.ddcompz
									U3((end-fwdpara_dis.np_extra+1:end),:,:) = labReceive(labindex + 1);
									labSend(U3((end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra),:,:),labindex+1)
								end
								if fwdpara_dis.ddcompz_id ~= 1
									U3(1:fwdpara_dis.np_extra,:,:) = labReceive(labindex-1);
								end
							end
							labBarrier
							if fwdpara.ddcompx~=1
								if fwdpara_dis.ddcompx_id ~= 1
									labSend(U3(:,fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra,:),labindex-fwdpara.ddcompz)
								end
								if fwdpara_dis.ddcompx_id ~= fwdpara.ddcompx
									U3(:,(end-fwdpara_dis.np_extra+1:end),:) = labReceive(labindex + fwdpara.ddcompz);
									labSend(U3(:,(end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra),:),labindex+fwdpara.ddcompz)
								end
								if fwdpara_dis.ddcompx_id ~= 1
									U3(:,1:fwdpara_dis.np_extra,:) =  labReceive(labindex-fwdpara.ddcompz);
								end
							end
							labBarrier

							if fwdpara.ddcompy~=1
								if fwdpara_dis.ddcompy_id ~= 1
									labSend(U3(:,:,fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra),labindex-fwdpara.ddcompz*fwdpara.ddcompx)
								end
								if fwdpara_dis.ddcompy_id ~= fwdpara.ddcompy
									U3(:,:,(end-fwdpara_dis.np_extra+1:end)) = labReceive(labindex + fwdpara.ddcompz*fwdpara.ddcompx);
									labSend(U3(:,:,(end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra)),labindex+fwdpara.ddcompz*fwdpara.ddcompx)
								end
								if fwdpara_dis.ddcompy_id ~= 1
									U3(:,:,1:fwdpara_dis.np_extra) = labReceive(labindex-fwdpara.ddcompz*fwdpara.ddcompx);
								end
							end
							labBarrier		
						end								
						% --------------- end of domain decomposition -----------------------


					end
			
					if labindex ==1 && mod(ti,50)==0 && fwdpara.fid;
						fprintf(fwdpara.fid,[repmat(' ',1,10),'Time=',num2str(fwdpara.taxis(ti)),'s, Total=',num2str(fwdpara.taxis(end)),'s.\n'])
					end

					% ---------------------------- display snapshot----------------------------
					if fwdpara.disp_snapshot == 1 & mod(ti,20)==0
						if length(fwdpara.mmy) ==1
							% 2D ploting
							figure(1);imagesc(U3);title([num2str(ti)]);colorbar;pause(.01)
						else
							% 3D ploting 
							figure(1);sliceview(reshape(U3,length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc)));
							title([num2str(ti)]);pause(.01)
						end
					end


					% shot record
					U_shot_record_temp(:,ti) = Pr * U3(:);
					
				
					% -------------- save check point per time step ---------------------
					% if chk_point_idx(ti) == 1
					% 	if fwdpara.chkp_save == 1
					% 		chk_point_wavefeild_temp(:,chk_point_count)	= U3(:);
					% 		chk_point_count = chk_point_count + 1;
					% 	elseif fwdpara.chkp_save == 2
					% 		chk_point_wavefeild_temp(:,chk_point_count)	= U3(:);
					% 		chk_point_count = chk_point_count + 1;
					% 	elseif fwdpara.chkp_save == 3
					% 		save_chkpoint_to_disk(U3,labindex,si,fwdpara.chkp_save,chk_point_count);
					% 		chk_point_count = chk_point_count + 1;
					% 	end
					% end
				
				
				
					U1	= U2;
					U2	= U3;
					U3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
					%
					% u1 = u2;
					% u2 = u3;
					% u3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
					%
				
				end
				% save the wavefield at the last time step
				U_shot_record_temp(:,ti)		= Pr * U2(:);
				output_data{si}				= U_shot_record_temp.';
				if strcmp(fwdpara.rtype,'marine')
					rec_idx{si}					= fwdpara_dis.rec_idx;
				end
				ntrace_loc = ntrace_loc + length(fwdpara_dis.rec_idx);
			
			
				% -------------- save check point per shot---------------------
				%
				% if fwdpara.chkp_save == 1
				% 	wavefield_chk_point{si,1} = chk_point_wavefeild_temp;
				% elseif fwdpara.chkp_save == 2
				%
				% 	save_chkpoint_to_disk(chk_point_wavefeild_temp,labindex,si,2,chk_point_count)
				% end

			
			end
		
			% different workers may simulated different number of shots (e.g, worker 1 simulate 15 shots, worker 2 simulate 16 shots),
			% 'labBarrier' can not work while other workers are idle. 
			% this make sure all worker finish simulation. You can arrange workers to save one shot time. 
			if Num_of_shots_loc < fwdpara_dis.max_number_of_shot_all_worker
				for si = Num_of_shots_loc+1: fwdpara_dis.max_number_of_shot_all_worker
					for ti = 1:Num_of_time_step 
						labBarrier
						labBarrier
						labBarrier
					end
				end
			end
		
			% put all data together
			codistr			= codistributor1d(2,codistributor1d.unsetPartition,[1,numlabs]);
			ntrace_loc		= codistributed.build(ntrace_loc,codistr,'noCommunication'); 
		
		
		end

		[output_data,loc_rec_idx] = output_shot_record(output_data,rec_idx,fwdpara,fwdpara_dis,ntrace_loc,loc_rec_idx);
	
	

	elseif fwdpara.fmode == -1
	% do adjoint backpropergation
	% loop from Nt to 1
	% do forward wavefield
	% u_n-2 = do_A1_and_adjoint(u_n-1,u_n,mode=2) + do_A2_and_adjoint(u_n,u_n-1,mode=2) + q_n-1;
	%
	% end of loop
	%
		if fwdpara.fid && labindex == 1
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n'])
			fprintf(fwdpara.fid,[repmat(' ',1,10),'Solving adjoint wavefield from Tmax to T1.\n'])
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n'])
		end
		
		
		spmd		
			%
			Num_of_shots_loc	= size(fwdpara_dis.slx_band,2);
			Num_of_time_step	= length(fwdpara.taxis);
			loc_rec_idx			= [];
			% locate receiver position (all shot has same receiver)
			Pr = 0;
			if strcmp(fwdpara.rtype,'full')
				fwdpara_dis.rlx	= fwdpara.rlx;
				fwdpara_dis.rly	= fwdpara.rly;
				fwdpara_dis.rlz	= fwdpara.rlz;
				fwdpara_dis		= locate_rcv_position_adj(fwdpara_dis);
				Pr				= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
												fwdpara_dis.rlz_loc,fwdpara_dis.rlx_loc,fwdpara_dis.rly_loc);
				rec_idx			= fwdpara_dis.rec_idx;
			elseif strcmp(fwdpara.rtype,'marine')
				rec_idx			= cell(Num_of_shots_loc,1);
			end
		


			% Num_of_rec_global	= size(fwdpara_dis.rlx,1);

			% allocate shot record
			output_data		= cell(1,Num_of_shots_loc);
		
		
			% generate time stepping stencil matrix
			% A1 = [];A1_inv = []; A2 = [];A3 = [];A3_inv = [];	
			[~,A1_inv,A2,A3,~] = generate_time_stepping_stencil_matrix(fwdpara_dis,fwdpara);
			% combine some matrix to same time
			Pr = A1_inv' *  Pr';
			A2 = A1_inv' *  A2';
			A3 = A1_inv' *  A3';

	
			% count total number of receivers for all loc shot
			ntrace_loc			 = 0;
			fwdpara_dis.Nr_count = 1;
			% ---------------------------- loop over shots ---------------------------------
			for si  = 1:Num_of_shots_loc
				% locate src position and receiver position (each shot has diff receivers) 
				if fwdpara.fid && labindex == 1
					fprintf(fwdpara.fid,[repmat(' ',1,5),repmat('-',1,5),' Simulating shot ',num2str(si),' ',repmat('-',1,5),'\n']);
				end
				fwdpara_dis.slx		= vec(fwdpara_dis.slx_band(:,si));
				fwdpara_dis.sly		= vec(fwdpara_dis.sly_band(:,si));
				fwdpara_dis.slz		= vec(fwdpara_dis.slz_band(:,si));
				fwdpara_dis 		= locate_shot_position_adj(fwdpara_dis);
				Ps					= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
												fwdpara_dis.slz_loc,fwdpara_dis.slx_loc,fwdpara_dis.sly_loc);
			
				if strcmp(fwdpara.rtype,'marine')
					fwdpara_dis.rlx = fwdpara.rlx(:,fwdpara_dis.shots_id(si));
					fwdpara_dis.rly = fwdpara.rly(:,fwdpara_dis.shots_id(si));
					fwdpara_dis.rlz = fwdpara.rlz(:,fwdpara_dis.shots_id(si));
					fwdpara_dis 	= locate_rcv_position_adj(fwdpara_dis);
					Pr				= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
													fwdpara_dis.rlz_loc,fwdpara_dis.rlx_loc,fwdpara_dis.rly_loc);
					% combine some matrix to same time
					Pr = A1_inv' *  Pr';
				end
			
				% generate source wavelet
				Q  					= prepare_source_wavelet(src,fwdpara,fwdpara_dis,si);
				
				% allocate place for loc shot record
				% U_shot_record_temp	= zeros(length(fwdpara_dis.rec_idx),Num_of_time_step);
				U_shot_record_temp	= zeros(size(fwdpara_dis.slx_loc,1),Num_of_time_step);

				



				
				% allocate wavefield snapshot 
				U1	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				U2	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				U3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
			
				% ---------------------------- loop over time step ---------------------------------
				for ti = Num_of_time_step:-1:1

					% adjoint modeling  solve for AT V = dU
					%
					% [A1T A2T A3T                   ]  V1    dU1
					% |            . . .             |  V2    dU2
					% |                   A1T A2T A3T|  V3  = dU3
					% |                       A1T A2T|  V4    dU4
					% [                           A1T]  V5    dU5
					%

					if ti == Num_of_time_step
						U3(1:end)   = full((Pr * Q(:,ti)));
					elseif ti == Num_of_time_step-1
						U3(1:end)   = Pr * Q(:,ti) - A2 * U2(:);
					else
						U3(1:end) 	= Pr * Q(:,ti) - A2 * U2(:) - A3 * U1(:);
					end
					% U3(1:end) = A1_inv' *  U3(:);
					if labindex ==1 && mod(ti,50)==0 && fwdpara.fid;
						fprintf(fwdpara.fid,[repmat(' ',1,10),'Time=',num2str(fwdpara.taxis(ti)),'s, Total=',num2str(fwdpara.taxis(end)),'s.\n'])
					end
					% norm(U3(:))
				
					% display snapshot
					if fwdpara.disp_snapshot == 1 & mod(ti,20)==0
						if length(fwdpara.mmy) ==1
							% 2D ploting
							figure(1);imagesc(U3);title([num2str(ti)]);pause(.01)
						else
							% 3D ploting 
							figure(1);sliceview(reshape(U3,length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc)));
							title(num2str(ti));pause(.01)
						end
					end

					% ------------ domain decomposition send data around ---------
					if fwdpara_dis.num_model ~=1
						if fwdpara.ddcompz~=1
							if fwdpara_dis.ddcompz_id ~= 1
								labSend(U3(1:fwdpara_dis.np_extra*2,:,:),labindex-1)
							end
							if fwdpara_dis.ddcompz_id ~= fwdpara.ddcompz
								U3(end-2*fwdpara_dis.np_extra+1:end,:,:) = U3(end-2*fwdpara_dis.np_extra+1:end,:,:) + labReceive(labindex + 1);
								labSend(U3(end-2*fwdpara_dis.np_extra+1:end,:,:),labindex+1)
							end
							if fwdpara_dis.ddcompz_id ~= 1
								U3(1:fwdpara_dis.np_extra*2,:,:) = labReceive(labindex-1);
							end
						end
						labBarrier
						if fwdpara.ddcompx~=1
							if fwdpara_dis.ddcompx_id ~= 1
								labSend(U3(:,1:fwdpara_dis.np_extra*2,:),labindex-fwdpara.ddcompz)
							end
							if fwdpara_dis.ddcompx_id ~= fwdpara.ddcompx
								U3(:,((end-2*fwdpara_dis.np_extra+1):end),:) = U3(:,((end-2*fwdpara_dis.np_extra+1):end),:) + labReceive(labindex + fwdpara.ddcompz);
								labSend(U3(:,(end-2*fwdpara_dis.np_extra+1:end),:),labindex+fwdpara.ddcompz)
							end
							if fwdpara_dis.ddcompx_id ~= 1
								U3(:,1:2*fwdpara_dis.np_extra,:) =  labReceive(labindex-fwdpara.ddcompz);
							end
						end
						labBarrier

						if fwdpara.ddcompy~=1
							if fwdpara_dis.ddcompy_id ~= 1
								labSend(U3(:,:,1:2*fwdpara_dis.np_extra),labindex-fwdpara.ddcompz*fwdpara.ddcompx)
							end
							if fwdpara_dis.ddcompy_id ~= fwdpara.ddcompy
								U3(:,:,(end-2*fwdpara_dis.np_extra+1:end)) = U3(:,:,(end-2*fwdpara_dis.np_extra+1):end) + labReceive(labindex + fwdpara.ddcompz*fwdpara.ddcompx);
								labSend(U3(:,:,(end-2*fwdpara_dis.np_extra+1):end),labindex+fwdpara.ddcompz*fwdpara.ddcompx)
							end
							if fwdpara_dis.ddcompy_id ~= 1
								U3(:,:,1:2*fwdpara_dis.np_extra) = labReceive(labindex-fwdpara.ddcompz*fwdpara.ddcompx);
							end
						end
						labBarrier		
					end								
					% --------------- end of domain decomposition -----------------------
					
				
					% output receiver wavefield
					U_shot_record_temp(:,ti) = Ps * U3(:);

					U1 = U2;
					U2 = U3;
					U3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));

				end
			
				% save the wavefield at the last time step
				% U_shot_record_temp(:,ti)	= Ps * U2(:);
				output_data{si}				= U_shot_record_temp.';
				if strcmp(fwdpara.rtype,'marine')
					rec_idx{si}					= fwdpara_dis.rec_idx;
				end
			
				fwdpara_dis.Nr_count 			= fwdpara_dis.Nr_count + length(fwdpara_dis.rec_idx);
				ntrace_loc = ntrace_loc + size(U_shot_record_temp,1);
			end
		
			% different workers may simulated different number of shots (e.g, worker 1 simulate 15 shots, worker 2 simulate 16 shots),
			% 'labBarrier' can not work while other workers are idle. 
			% this make sure all worker finish simulation. You can arrange workers to save one shot time. 
			if Num_of_shots_loc < fwdpara_dis.max_number_of_shot_all_worker
				for si = Num_of_shots_loc+1: fwdpara_dis.max_number_of_shot_all_worker
					for ti = 1:Num_of_time_step 
						labBarrier
						labBarrier
						labBarrier
					end
				end
			end
		
			% put all data together
		
			% put all data together
			codistr			= codistributor1d(2,codistributor1d.unsetPartition,[1,numlabs]);
			ntrace_loc		= codistributed.build(ntrace_loc,codistr,'noCommunication'); 
	
		
		end
		[output_data,~] = output_shot_record(output_data,rec_idx,fwdpara,fwdpara_dis,ntrace_loc,0);
	
	
	
	%%==========================================================================================================
	%					Jacobian of FWI objective function ||D-U||_2
	% ==========================================================================================================
	%
	% 	Time domain stepping modeling can be form as 
	% 		AU = Q;
	%	Take the derivative 
	%      dA            dU
	%    ---- U  =  -A ----
	%    dm_i          dm_i
	% ----------------------------------------------------------------------------------------------------------
	%	Forward mode: Born modeling operator
	% 
	%	       N      -1   dA                 -1       dA    -1
	%	 dU = sum [ -A   (---- U) dm_i ]  = -A  {diag[---- (A  Q)] dm }
	%	       i          dm_i                         dm
	% 
	%	Note: notice background wavefield can be solved time step by step as virtul source
	%
	% ----------------------------------------------------------------------------------------------------------
	%	Adjoint mode: RTM operator
	% 
	%	        src and rec               dA    -1         -1 T          
	%	 g  =       sum      - conj{diag[---- (A  Q)]}.*[(A  )  dU];
	%	                                  dm                   
	%	         -1 T      T -1
	%	 Note: (A  )  =  (A )  ; 
	%
	% 	Note: Background wavefield, forward modeling, solve for A U  = Q
	% ----------------------------------------------------------------------------------------------------------
	%
	% [ A1                           ]  U1      Q1
	% | A2 A1                        |  U2      Q2
	% | A3 A2 A1                     |  U3      Q3
	% |             . . .            |  .     =  .
	% |                   A3 A2 A1   |  Un-1    Q4
	% [                      A3 A2 A1]  Un      Q5
	% 
	%
	% Back propergation wavefield, adjoint modeling,  solve for AT V = dU
	%
	% [A1T A2T A3T                   ]  V1      dU1
	% |    A1T A2T A3T               |  V2      dU2
	% |            . . .             |  .        .
	% |                   A1T A2T A3T|  Vn-2  = dUn-2
	% |                       A1T A2T|  Vn-1    dUn-1
	% [                           A1T]  Vn      dUn
	% 
	% do checkpoints, becarefull with the memory.
	% Usually, forward wavefield is solved from T1 to Tmax
	% But if we know "Un-1" and "U_n", Then we can from Tmax to T1
	%
	elseif fwdpara.fmode == 3
		%	Forward mode: Born modeling operator
		% ---------------------------------------------------------------------- 
		%	       N      -1   dA                 -1       dA    -1
		%	 dU = sum [ -A   (---- U) dm_i ]  = -A  {diag[---- (A  Q)] dm }
		%	       i          dm_i                         dm
		% 
		% ----------------------------------------------------------------------
		if fwdpara.fid && labindex == 1
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n'])
			fprintf(fwdpara.fid,[repmat(' ',1,10),'Born modeling.\n'])
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n'])
		end
	
		spmd		
			%
			ntrace_loc 			= 0; % number of loc traces, Nr_loc * Ns_loc; 
			Num_of_shots_loc	= size(fwdpara_dis.slx_band,2);
			Num_of_time_step	= length(fwdpara.taxis);
			loc_rec_idx			= [];
			
			% locate receiver position (all shot has same receiver)
			if strcmp(fwdpara.rtype,'full')
				fwdpara_dis.rlx	= fwdpara.rlx;
				fwdpara_dis.rly	= fwdpara.rly;
				fwdpara_dis.rlz	= fwdpara.rlz;
				fwdpara_dis		= locate_rcv_position(fwdpara_dis);
				Pr				= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
												fwdpara_dis.rlz_loc,fwdpara_dis.rlx_loc,fwdpara_dis.rly_loc);
				rec_idx			= fwdpara_dis.rec_idx;
			elseif strcmp(fwdpara.rtype,'marine')
				rec_idx			= cell(Num_of_shots_loc,1);
			end

			% Num_of_rec_global	= size(fwdpara_dis.rlx,1);
			
			% create place for shot record
			output_data		= cell(1,Num_of_shots_loc);
			
	
			% generate time stepping stencil matrix
			% A1 = [];A1_inv = []; A2 = [];A3 = [];A3_inv = [];
			[~,A1_inv,A2,A3,~,dA1,dA2,dA3] = generate_time_stepping_stencil_matrix(fwdpara_dis,fwdpara);
			% combine some matrix to save time
			A2	= A1_inv * A2;
			A3	= A1_inv * A3;
			dA1	= A1_inv * dA1;
			dA2	= A1_inv * dA2;
			dA3	= A1_inv * dA3;

			
			
			% ---------------------------------- loop over shots ---------------------------------------
			for si  = 1:Num_of_shots_loc
				if fwdpara.fid && labindex == 1
					fprintf(fwdpara.fid,[repmat(' ',1,5),repmat('-',1,5),' Simulating shot ',num2str(si),' ',repmat('-',1,5),'\n']);
				end				
				% locate src position and receiver position (each shot has diff receivers) 
				fwdpara_dis.slx		= vec(fwdpara_dis.slx_band(:,si));
				fwdpara_dis.sly		= vec(fwdpara_dis.sly_band(:,si));
				fwdpara_dis.slz		= vec(fwdpara_dis.slz_band(:,si));
				fwdpara_dis 		= locate_shot_position(fwdpara_dis);
				Ps					= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
												fwdpara_dis.slz_loc,fwdpara_dis.slx_loc,fwdpara_dis.sly_loc);
				
				if strcmp(fwdpara.rtype,'marine')
					fwdpara_dis.rlx	= fwdpara.rlx(:,fwdpara_dis.shots_id(si));
					fwdpara_dis.rly	= fwdpara.rly(:,fwdpara_dis.shots_id(si));
					fwdpara_dis.rlz	= fwdpara.rlz(:,fwdpara_dis.shots_id(si));
					fwdpara_dis 	= locate_rcv_position(fwdpara_dis);
					Pr				= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
													fwdpara_dis.rlz_loc,fwdpara_dis.rlx_loc,fwdpara_dis.rly_loc);
				end
				
				loc_rec_idx = [loc_rec_idx(:);fwdpara_dis.rec_idx + size(fwdpara.rlx,1).*(fwdpara_dis.shots_id(si) - 1)];
				% generate source wavelet
				Q  					= prepare_source_wavelet(src,fwdpara,fwdpara_dis,si);
			
				% % allocate place for loc shot record
				U_shot_record_temp	= zeros(length(fwdpara_dis.rec_idx),Num_of_time_step);
				%
				
				
				% combine some matrix to save time
				Ps = A1_inv * Ps';


				% allocate wavefield snapshot
				U1	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				U2	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				U3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				
				dU1	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				dU2	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				dU3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				
								
				% % allocate space for checkpoint
				% if fwdpara.chkp_save == 1 | fwdpara.chkp_save == 2
				% 	chk_point_wavefeild_temp	= zeros(numel(U3),N_chkp);
				% end
				% chk_point_count				= 1;			


				% ---------------------------- loop over time step ---------------------------------
				for ti = 1:Num_of_time_step

					%  forward modeling solve for A U  = Q
					%
					% [ A1                           ]  U1    Q1
					% | A2 A1                        |  U2    Q2
					% | A3 A2 A1                     |  U3  = Q3
					% |             . . .            |  U4    Q4
					% [                      A3 A2 A1]  U5    Q5
					if ti == 1
						U3(1:end)   = full(Ps * Q(:,ti));
					elseif ti == 2
						U3(1:end)   = Ps * Q(:,ti) - A2 * U2(:);
					else
						U3(1:end) 	= Ps * Q(:,ti) - A2 * U2(:) - A3 * U1(:);
					end
					% U3(:) = A1_inv * U3(:);
					if labindex ==1 && mod(ti,50)==0 && fwdpara.fid;
						fprintf(fwdpara.fid,[repmat(' ',1,10),'Time=',num2str(fwdpara.taxis(ti)),'s, Total=',num2str(fwdpara.taxis(end)),'s.\n'])
					end
				

					% ---------------- domain decomposition send data around -----------------
					if fwdpara_dis.num_model ~=1
						if fwdpara.ddcompz~=1
							if fwdpara_dis.ddcompz_id ~= 1
								labSend(U3(fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra,:,:),labindex-1)
							end
							if fwdpara_dis.ddcompz_id ~= fwdpara.ddcompz
								U3((end-fwdpara_dis.np_extra+1:end),:,:) = labReceive(labindex + 1);
								labSend(U3((end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra),:,:),labindex+1)
							end
							if fwdpara_dis.ddcompz_id ~= 1
								U3(1:fwdpara_dis.np_extra,:,:) = labReceive(labindex-1);
							end
						end
						labBarrier
						if fwdpara.ddcompx~=1
							if fwdpara_dis.ddcompx_id ~= 1
								labSend(U3(:,fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra,:),labindex-fwdpara.ddcompz)
							end
							if fwdpara_dis.ddcompx_id ~= fwdpara.ddcompx
								U3(:,(end-fwdpara_dis.np_extra+1:end),:) = labReceive(labindex + fwdpara.ddcompz);
								labSend(U3(:,(end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra),:),labindex+fwdpara.ddcompz)
							end
							if fwdpara_dis.ddcompx_id ~= 1
								U3(:,1:fwdpara_dis.np_extra,:) =  labReceive(labindex-fwdpara.ddcompz);
							end
						end
						labBarrier

						if fwdpara.ddcompy~=1
							if fwdpara_dis.ddcompy_id ~= 1
								labSend(U3(:,:,fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra),labindex-fwdpara.ddcompz*fwdpara.ddcompx)
							end
							if fwdpara_dis.ddcompy_id ~= fwdpara.ddcompy
								U3(:,:,(end-fwdpara_dis.np_extra+1:end)) = labReceive(labindex + fwdpara.ddcompz*fwdpara.ddcompx);
								labSend(U3(:,:,(end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra)),labindex+fwdpara.ddcompz*fwdpara.ddcompx)
							end
							if fwdpara_dis.ddcompy_id ~= 1
								U3(:,:,1:fwdpara_dis.np_extra) = labReceive(labindex-fwdpara.ddcompz*fwdpara.ddcompx);
							end
						end
						labBarrier	
					end									
					% --------------- end of domain decomposition -----------------------

					% do  U_temp =  diag(dA U)*dm
					if ti == 1 
						U_temp = (dA1 * U3(:)).* (-fwdpara_dis.dvel_dis(:));
					elseif ti == 2
						U_temp = (dA2 * U2(:) + dA1 * U3(:)).* (-fwdpara_dis.dvel_dis(:));
					else
						U_temp = (dA3 * U1(:) + dA2 * U2(:) + dA1 * U3(:)).* (-fwdpara_dis.dvel_dis(:));
					end
					
					% do A-1 * U_temp
					
					if ti == 1
						dU3(1:end) = ( U_temp);
					elseif ti == 2
						dU3(1:end) = ((U_temp - A2 * dU2(:)));
					else
						dU3(1:end) 	= ((U_temp - A2 * dU2(:) - A3 * dU1(:)));
					end
					
					% ---------------- domain decomposition send data around -----------------
					if fwdpara_dis.num_model ~=1
						if fwdpara.ddcompz~=1
							if fwdpara_dis.ddcompz_id ~= 1
								labSend(dU3(1:fwdpara_dis.np_extra*2,:,:),labindex-1)
							end
							if fwdpara_dis.ddcompz_id ~= fwdpara.ddcompz
								dU3(end-2*fwdpara_dis.np_extra+1:end,:,:) = dU3(end-2*fwdpara_dis.np_extra+1:end,:,:) + labReceive(labindex + 1);
								labSend(dU3(end-2*fwdpara_dis.np_extra+1:end,:,:),labindex+1)
							end
							if fwdpara_dis.ddcompz_id ~= 1
								dU3(1:fwdpara_dis.np_extra*2,:,:) = labReceive(labindex-1);
							end
						end
						labBarrier
						if fwdpara.ddcompx~=1
							if fwdpara_dis.ddcompx_id ~= 1
								labSend(dU3(:,1:fwdpara_dis.np_extra*2,:),labindex-fwdpara.ddcompz)
							end
							if fwdpara_dis.ddcompx_id ~= fwdpara.ddcompx
								dU3(:,((end-2*fwdpara_dis.np_extra+1):end),:) = dU3(:,((end-2*fwdpara_dis.np_extra+1):end),:) + labReceive(labindex + fwdpara.ddcompz);
								labSend(dU3(:,(end-2*fwdpara_dis.np_extra+1:end),:),labindex+fwdpara.ddcompz)
							end
							if fwdpara_dis.ddcompx_id ~= 1
								dU3(:,1:2*fwdpara_dis.np_extra,:) =  labReceive(labindex-fwdpara.ddcompz);
							end
						end
						labBarrier

						if fwdpara.ddcompy~=1
							if fwdpara_dis.ddcompy_id ~= 1
								labSend(dU3(:,:,1:2*fwdpara_dis.np_extra),labindex-fwdpara.ddcompz*fwdpara.ddcompx)
							end
							if fwdpara_dis.ddcompy_id ~= fwdpara.ddcompy
								dU3(:,:,(end-2*fwdpara_dis.np_extra+1:end)) = dU3(:,:,(end-2*fwdpara_dis.np_extra+1):end) + labReceive(labindex + fwdpara.ddcompz*fwdpara.ddcompx);
								labSend(dU3(:,:,(end-2*fwdpara_dis.np_extra+1):end),labindex+fwdpara.ddcompz*fwdpara.ddcompx)
							end
							if fwdpara_dis.ddcompy_id ~= 1
								dU3(:,:,1:2*fwdpara_dis.np_extra) = labReceive(labindex-fwdpara.ddcompz*fwdpara.ddcompx);
							end
						end
						labBarrier	
					end									
					% --------------- end of domain decomposition -----------------------


					
 					% ---------------------------- display snapshot----------------------------
 					if fwdpara.disp_snapshot == 1 & mod(ti,20)==0
 						if length(fwdpara.mmy) ==1
 							% 2D ploting
 							figure(1);imagesc(dU3);title([num2str(ti)]);colorbar; pause(.01)
 						else
 							% 3D ploting 
 							figure(1);sliceview(reshape(dU3,length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc)));
 							title([num2str(ti)]);pause(.01)
 						end
 					end
					
				
			 					%
 					% % -------------- save check point per time step ---------------------
 					% if chk_point_idx(ti) == 1
 					% 	if fwdpara.chkp_save == 1
 					% 		chk_point_wavefeild_temp(:,chk_point_count)	= U3(:);
 					% 		chk_point_count 	= chk_point_count + 1;
 					% 	elseif fwdpara.chkp_save == 2
 					% 		chk_point_wavefeild_temp(:,chk_point_count)	= U3(:);
 					% 		chk_point_count 	= chk_point_count + 1;
 					% 	elseif fwdpara.chkp_save == 3
 					% 		save_chkpoint_to_disk(U3,labindex,si,fwdpara.chkp_save,chk_point_count);
 					% 		chk_point_count 	= chk_point_count + 1;
 					% 	end
 					% end
 					%

					
					% % output receiver wavefield
					U_shot_record_temp(:,ti) = Pr * dU3(:);
			%
					
					
				% reallocate wavefield.	
					U1	= U2;
					U2	= U3;
					U3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
					dU1	= dU2;
					dU2	= dU3;
					dU3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));

					
				end
		
				% save the wavefield at the last time step
				% U_shot_record_temp(:,ti)	= Pr * U2(:);
				output_data{si}				= U_shot_record_temp.';
				if strcmp(fwdpara.rtype,'marine')
					rec_idx{si}					= fwdpara_dis.rec_idx;
				end
				ntrace_loc = ntrace_loc + length(fwdpara_dis.rec_idx);
				
				

				% % -------------- save check point per shot---------------------
				% if fwdpara.chkp_save == 1
				% 	wavefield_chk_point{si,1} = chk_point_wavefeild_temp;
				% elseif fwdpara.chkp_save == 2
				%
				% 	save_chkpoint_to_disk(chk_point_wavefeild_temp,labindex,si,2,chk_point_count)
				% end

			end
			
			% different workers may simulated different number of shots (e.g, worker 1 simulate 15 shots, worker 2 simulate 16 shots),
			% 'labBarrier' can not work while other workers are idle. 
			% this make sure all worker finish simulation. You can arrange workers to save one shot time. 
			
			
			if Num_of_shots_loc < fwdpara_dis.max_number_of_shot_all_worker
				for si = Num_of_shots_loc+1: fwdpara_dis.max_number_of_shot_all_worker
					for ti = 1:Num_of_time_step 
						labBarrier
						labBarrier
						labBarrier
						labBarrier
						labBarrier
						labBarrier
					end
				end
			end
			
			% put all data together
			codistr			= codistributor1d(2,codistributor1d.unsetPartition,[1,numlabs]);
			ntrace_loc		= codistributed.build(ntrace_loc,codistr,'noCommunication'); 
			
		
		end
	
		[output_data,loc_rec_idx] = output_shot_record(output_data,rec_idx,fwdpara,fwdpara_dis,ntrace_loc,loc_rec_idx);
		
	elseif fwdpara.fmode == -3
		%	Adjoint mode: RTM operator
		% ---------------------------------------------------------------------- 
		%	        src and rec               dA    -1         -1 T          
		%	 g  =       sum      - conj{diag[---- (A  Q)]}.*[(A  )  dU];
		%	                                  dm                   
		%	         -1 T      T -1
		%	 Note: (A  )  =  (A )  ; 
		% ---------------------------------------------------------------------- 
		if fwdpara.fid && labindex == 1
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n'])
			fprintf(fwdpara.fid,[repmat(' ',1,10),'Gradient of FWI objective function.\n'])
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n'])
		end
		
		
		spmd		
			%
			
			ntrace_loc 			= 0; % number of loc traces, Nr_loc * Ns_loc; 
			Num_of_shots_loc	= size(fwdpara_dis.slx_band,2);
			Num_of_time_step	= length(fwdpara.taxis);
			loc_rec_idx			= [];
			% locate receiver position (all shot has same receiver)
			if strcmp(fwdpara.rtype,'full')
				fwdpara_dis.rlx	= fwdpara.rlx;
				fwdpara_dis.rly	= fwdpara.rly;
				fwdpara_dis.rlz	= fwdpara.rlz;
				fwdpara_dis		= locate_rcv_position(fwdpara_dis);
				Pr				= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
												fwdpara_dis.rlz_loc,fwdpara_dis.rlx_loc,fwdpara_dis.rly_loc);
				rec_idx			= fwdpara_dis.rec_idx;
			elseif strcmp(fwdpara.rtype,'marine')
				rec_idx			= cell(Num_of_shots_loc,1);
			end

			% Num_of_rec_global	= size(fwdpara_dis.rlx,1);
		
			% create place for shot record
			% U_shot_record		= cell(Num_of_shots_loc,1);
		
			% create place for checkpoint wavefield 
			% if fwdpara.chkp_save == 1
			% 	wavefield_chk_point	= cell(Num_of_shots_loc,1);
			% end
			
			[chk_point_idx,N_chkp]	= generate_check_point_idx(Num_of_time_step,fwdpara.chkp_space);
		
			
			% generate time stepping stencil matrix
			% A1 = [];A1_inv = []; A2 = [];A3 = [];A3_inv = [];	
			[A1,A1_inv,A2,A3,A3_inv,dA1,dA2,dA3] = generate_time_stepping_stencil_matrix(fwdpara_dis,fwdpara);
			
			
			% allocate place for gradient. 
			g   = zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc)); 
			
			% count total number of receivers for all loc shot
			% ntrace_loc			 = 0;
			fwdpara_dis.Nr_count = 1;
			
			% ---------------------------- loop over shots ---------------------------------
			for si  = 1:Num_of_shots_loc
				if fwdpara.fid && labindex == 1
					fprintf(fwdpara.fid,[repmat(' ',1,5),repmat('-',1,5),' Simulating shot ',num2str(si),' ',repmat('-',1,5),'\n']);
				end			
				% locate src position and receiver position (each shot has diff receivers) 
				fwdpara_dis.slx		= vec(fwdpara_dis.slx_band(:,si));
				fwdpara_dis.sly		= vec(fwdpara_dis.sly_band(:,si));
				fwdpara_dis.slz		= vec(fwdpara_dis.slz_band(:,si));
				fwdpara_dis 		= locate_shot_position(fwdpara_dis);
				Ps					= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
												fwdpara_dis.slz_loc,fwdpara_dis.slx_loc,fwdpara_dis.sly_loc);
			
				if strcmp(fwdpara.rtype,'marine')
					fwdpara_dis.rlx = fwdpara.rlx(:,fwdpara_dis.shots_id(si));
					fwdpara_dis.rly = fwdpara.rly(:,fwdpara_dis.shots_id(si));
					fwdpara_dis.rlz = fwdpara.rlz(:,fwdpara_dis.shots_id(si));
					fwdpara_dis 	= locate_rcv_position(fwdpara_dis);
					Pr				= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
													fwdpara_dis.rlz_loc,fwdpara_dis.rlx_loc,fwdpara_dis.rly_loc);
				end
			
				% generate source wavelet
				Q 					 = prepare_source_wavelet(src,fwdpara,fwdpara_dis,si);
				labbarrier_count	 = 0;
				% allocate place for loc shot record
				% U_shot_record_temp	= zeros(length(fwdpara_dis.rec_idx),Num_of_time_step);
				
				% combine some matrix to save time
				
				
			

				% allocate wavefield snapshot
				U1	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				U2	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				U3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				
				% computer adjoint wavefield
				V1	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				V2	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				V3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
		
				
				if fwdpara.chkp_save == 1 
					chk_point_wavefeild_temp = wavefield_chk_point{si,1};
				elseif fwdpara.chkp_save == 2
					chk_point_wavefeild_temp = load_chkpoint_from_disk(labindex,si,fwdpara.chkp_save);
				end
				chk_point_count				 = N_chkp;			
		
				% initialize the last two wavefield
				if fwdpara.chkp_save == 1
					U1(1:end) 				= reshape(chk_point_wavefeild_temp(:,chk_point_count),size(U2));
					chk_point_count			= chk_point_count - 1;
					U2(1:end) 				= reshape(chk_point_wavefeild_temp(:,chk_point_count),size(U2));
					chk_point_count			= chk_point_count - 1;
				elseif fwdpara.chkp_save == 2
					U1(1:end)  				= reshape(chk_point_wavefeild_temp(:,chk_point_count),size(U2));
					chk_point_count 		= chk_point_count - 1;
					U2(1:end)  				= reshape(chk_point_wavefeild_temp(:,chk_point_count),size(U2));
					chk_point_count 		= chk_point_count - 1;
				elseif fwdpara.chkp_save == 3
					U1(1:end) 			    = load_chkpoint_from_disk(labindex,si,fwdpara.chkp_save,chk_point_count);
					chk_point_count			= chk_point_count - 1;
					U2(1:end) 			    = load_chkpoint_from_disk(labindex,si,fwdpara.chkp_save,chk_point_count);
					chk_point_count			= chk_point_count - 1;
				end

		
				% ---------------------------- loop over time step ---------------------------------
				for ti = Num_of_time_step:-1:1

					%  forward modeling solve for A U  = Q
					%
					% [ A1                           ]  U1    Q1
					% | A2 A1                        |  U2    Q2
					% | A3 A2 A1                     |  U3  = Q3
					% |             . . .            |  U4    Q4
					% [                      A3 A2 A1]  U5    Q5
					

					% wavefield completed while ti = 3
					if ti > 2
					% -------------- load check point per time step ---------------------
						if chk_point_idx(ti-2) == 1
							if fwdpara.chkp_save == 1
								U3(1:end) 				= reshape(chk_point_wavefeild_temp(:,chk_point_count),size(U2));
								chk_point_count			= chk_point_count - 1;
							elseif fwdpara.chkp_save == 2
								U3(1:end)  				= reshape(chk_point_wavefeild_temp(:,chk_point_count),size(U2));
								chk_point_count 		= chk_point_count - 1;
							elseif fwdpara.chkp_save == 3
								U3(1:end) 			    = load_chkpoint_from_disk(labindex,si,fwdpara.chkp_save,chk_point_count);
								chk_point_count			= chk_point_count - 1;
							end
						
						else
							U3(1:end)	 =	Ps'*Q(:,ti) - A2*U2(:) - A1 * U1(:);		
							U3(1:end)	=	A3_inv * U3(:);	
							% ---------------- domain decomposition send data around -----------------
							if fwdpara_dis.num_model ~=1
								if fwdpara.ddcompz~=1
									if fwdpara_dis.ddcompz_id ~= 1
										labSend(U3(fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra,:,:),labindex-1)
									end
									if fwdpara_dis.ddcompz_id ~= fwdpara.ddcompz
										U3((end-fwdpara_dis.np_extra+1:end),:,:) = labReceive(labindex + 1);
										labSend(U3((end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra),:,:),labindex+1)
									end
									if fwdpara_dis.ddcompz_id ~= 1
										U3(1:fwdpara_dis.np_extra,:,:) = labReceive(labindex-1);
									end
								end
								labBarrier
								if fwdpara.ddcompx~=1
									if fwdpara_dis.ddcompx_id ~= 1
										labSend(U3(:,fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra,:),labindex-fwdpara.ddcompz)
									end
									if fwdpara_dis.ddcompx_id ~= fwdpara.ddcompx
										U3(:,(end-fwdpara_dis.np_extra+1:end),:) = labReceive(labindex + fwdpara.ddcompz);
										labSend(U3(:,(end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra),:),labindex+fwdpara.ddcompz)
									end
									if fwdpara_dis.ddcompx_id ~= 1
										U3(:,1:fwdpara_dis.np_extra,:) =  labReceive(labindex-fwdpara.ddcompz);
									end
								end
								labBarrier

								if fwdpara.ddcompy~=1
									if fwdpara_dis.ddcompy_id ~= 1
										labSend(U3(:,:,fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra),labindex-fwdpara.ddcompz*fwdpara.ddcompx)
									end
									if fwdpara_dis.ddcompy_id ~= fwdpara.ddcompy
										U3(:,:,(end-fwdpara_dis.np_extra+1:end)) = labReceive(labindex + fwdpara.ddcompz*fwdpara.ddcompx);
										labSend(U3(:,:,(end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra)),labindex+fwdpara.ddcompz*fwdpara.ddcompx)
									end
									if fwdpara_dis.ddcompy_id ~= 1
										U3(:,:,1:fwdpara_dis.np_extra) = labReceive(labindex-fwdpara.ddcompz*fwdpara.ddcompx);
									end
								end
								labBarrier										
								labbarrier_count  = labbarrier_count+1;	
							end								
							% --------------- end of domain decomposition -----------------------
						end
					end % end of if
					% do  U_temp =  diag(dA U)
					if ti == 1
						U_temp = -(dA1 * U1(:));
					elseif ti == 2
						U_temp = -(dA2 * U2(:) + dA1 * U1(:));
					else
						U_temp = -(dA3 * U3(:) + dA2 * U2(:) + dA1 * U1(:));
					end
			
					
					if labindex ==1 && mod(ti,50)==0 && fwdpara.fid;
						fprintf(fwdpara.fid,[repmat(' ',1,10),'Time=',num2str(fwdpara.taxis(ti)),'s, Total=',num2str(fwdpara.taxis(end)),'s.\n'])
					end


					% shot record
					% U_shot_record_temp(:,ti) = Pr * U3(:);
					
					
					% generate source wavelet
					du  					= prepare_adjoint_source(input_data,fwdpara,fwdpara_dis,si);

					% adjoint modeling  solve for AT V = dU
					%
					% [A1T A2T A3T                   ]  V1    dU1
					% |            . . .             |  V2    dU2
					% |                   A1T A2T A3T|  V3  = dU3
					% |                       A1T A2T|  V4    dU4
					% [                           A1T]  V5    dU5
					%

					if ti == Num_of_time_step
						V3(1:end)   = full( (Pr' * du(:,ti)));
					elseif ti == Num_of_time_step-1
						V3(1:end)   = Pr' * du(:,ti) - A2' * V2(:);
					else
						V3(1:end) 	= Pr' * du(:,ti) - A2' * V2(:) - A3' * V1(:);
					end 
					V3(1:end) = A1_inv' *  V3(:);					
					% ------------ domain decomposition send data around ---------
					if fwdpara_dis.num_model ~=1
						if fwdpara.ddcompz~=1
							if fwdpara_dis.ddcompz_id ~= 1
								labSend(V3(1:fwdpara_dis.np_extra*2,:,:),labindex-1)
							end
							if fwdpara_dis.ddcompz_id ~= fwdpara.ddcompz
								V3(end-2*fwdpara_dis.np_extra+1:end,:,:) = V3(end-2*fwdpara_dis.np_extra+1:end,:,:) + labReceive(labindex + 1);
								labSend(V3(end-2*fwdpara_dis.np_extra+1:end,:,:),labindex+1)
							end
							if fwdpara_dis.ddcompz_id ~= 1
								V3(1:fwdpara_dis.np_extra*2,:,:) = labReceive(labindex-1);
							end
						end
						labBarrier
						if fwdpara.ddcompx~=1
							if fwdpara_dis.ddcompx_id ~= 1
								labSend(V3(:,1:fwdpara_dis.np_extra*2,:),labindex-fwdpara.ddcompz)
							end
							if fwdpara_dis.ddcompx_id ~= fwdpara.ddcompx
								V3(:,((end-2*fwdpara_dis.np_extra+1):end),:) = V3(:,((end-2*fwdpara_dis.np_extra+1):end),:) + labReceive(labindex + fwdpara.ddcompz);
								labSend(V3(:,(end-2*fwdpara_dis.np_extra+1:end),:),labindex+fwdpara.ddcompz)
							end
							if fwdpara_dis.ddcompx_id ~= 1
								V3(:,1:2*fwdpara_dis.np_extra,:) =  labReceive(labindex-fwdpara.ddcompz);
							end
						end
						labBarrier

						if fwdpara.ddcompy~=1
							if fwdpara_dis.ddcompy_id ~= 1
								labSend(V3(:,:,1:2*fwdpara_dis.np_extra),labindex-fwdpara.ddcompz*fwdpara.ddcompx)
							end
							if fwdpara_dis.ddcompy_id ~= fwdpara.ddcompy
								V3(:,:,(end-2*fwdpara_dis.np_extra+1:end)) = V3(:,:,(end-2*fwdpara_dis.np_extra+1):end) + labReceive(labindex + fwdpara.ddcompz*fwdpara.ddcompx);
								labSend(V3(:,:,(end-2*fwdpara_dis.np_extra+1):end),labindex+fwdpara.ddcompz*fwdpara.ddcompx)
							end
							if fwdpara_dis.ddcompy_id ~= 1
								V3(:,:,1:2*fwdpara_dis.np_extra) = labReceive(labindex-fwdpara.ddcompz*fwdpara.ddcompx);
							end
						end
						labBarrier
						labbarrier_count  = labbarrier_count+1;	
					end							
					% --------------- end of domain decomposition -----------------------
					
					g = g + reshape(U_temp,size(V3)).*V3;
					% ---------------------------- display snapshot----------------------------
					if fwdpara.disp_snapshot == 1 & mod(ti,10)==0 
						if length(fwdpara.mmy) ==1
							% 2D ploting
							figure(1);%colormap(promax)
							subplot(3,1,1);imagesc(U3);title([num2str(ti)]);caxis(1e0*[-1 1])
							subplot(3,1,2);imagesc(V3);title([num2str(ti)]);caxis([-20 20])
							subplot(3,1,3);imagesc(g);title([num2str(ti)]);caxis(1e3*[-3 3])
						else
							% 3D ploting 
							figure(1);%colormap(promax)
							subplot(3,1,1);sliceview(reshape(-U3,length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc)));
							subplot(3,1,2);sliceview(reshape(-V3,length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc)));
							subplot(3,1,3);sliceview(reshape(-g,length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc)));
							title([num2str(ti)]);pause(.01)
						end
					end
				
					V1	= V2;
					V2	= V3;
					V3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
					
					
					
					if ti == 1
						U1	= U2;
					elseif ti == 2
						U1	= U2;
						U2	= U3;
					else					
						U1	= U2;
						U2	= U3;
						U3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
					end
					%
					% u1 = u2;
					% u2 = u3;
					% u3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
					%
				
				end % end of time loop
				% save the wavefield at the last time step
				% U_shot_record_temp(:,ti)		= Pr * U2(:);
				% U_shot_record{si}				= U_shot_record_temp;
				if strcmp(fwdpara.rtype,'marine')
					rec_idx{si}					= fwdpara_dis.rec_idx;
				end
				% ntrace_loc = ntrace_loc + length(fwdpara_dis.rec_idx);
				
				
				fwdpara_dis.Nr_count 			= fwdpara_dis.Nr_count + length(fwdpara_dis.rec_idx);
				% ntrace_loc = ntrace_loc + size(U_shot_record_temp,1);
			
				
				% -------------- save check point per shot---------------------
				%
				% if fwdpara.chkp_save == 1
				% 	wavefield_chk_point{si,1} = chk_point_wavefeild_temp;
				% elseif fwdpara.chkp_save == 2
				%
				% 	save_chkpoint_to_disk(chk_point_wavefeild_temp,labindex,si,2,chk_point_count)
				% end

			
			end % end of shot loop

			% different workers may simulated different number of shots (e.g, worker 1 simulate 15 shots, worker 2 simulate 16 shots),
			% 'labBarrier' can not work while other workers are idle. 
			% this make sure all worker finish simulation. You can arrange workers to save one shot time. 
			if Num_of_shots_loc < fwdpara_dis.max_number_of_shot_all_worker
				for si = Num_of_shots_loc+1: fwdpara_dis.max_number_of_shot_all_worker
					for ti = 1:labbarrier_count*3;
						labBarrier
					end
				end
			end
		
			% put all data together
			% codistr			= codistributor1d(2,codistributor1d.unsetPartition,[1,numlabs]);
			% ntrace_loc		= codistributed.build(ntrace_loc,codistr,'noCommunication');
			%
			
			g				= g(fwdpara_dis.u_idz,fwdpara_dis.u_idx,fwdpara_dis.u_idy);

		end
		
		output_data 			= gather_gradient(g,fwdpara_dis,fwdpara);
		% [U_shot_record,shot_idx_all] = output_shot_record(U_shot_record,rec_idx,fwdpara,fwdpara_dis,ntrace_loc);
	

		%% =========================================================================
		%			TTI acoustic wave equation
		%% =========================================================================
		%   1    d^2 U                            d^2 U
		% ----- ------- - (1 + 2 sigma) H U - H0 -------- = (1+2 sigma) H V + Q
		%  c^2   d t^2                            d z^2
		%
		%   1    d^2 V
		% ----- ------- - 2(epsilo-sigma)H V  = 2 (epsilo-sigma) H U
		%  c^2   d t^2
		%
		% ---------------------------------- 2D case -------------------------------
		%
		%                      d^2                     d^2                    d^2
		%  H = cos^2 (theta) -------  + sin^2(theta) ------- - sin(2 theta) ------
		%                     d x^2                   d z^2                  dx dz
		%                      d^2                    d^2                     d^2
		%  H0 = sin^2(theta) -------  + cos^2(theta) -------- + sin(2 theta) ------
		%                     d x^2                   d z^2                   dx dz
		%
		% ---------------------------------- 3D case --------------------------------
		%         d^2                 d^2                d^2         d^2        d^2 
		% H = A ------- + (B + A D) ------- + (A - DB) ------- - CH ----- - CG -----
		%        d x^2               d y^2              d z^2        dxdy       dxdz
		%           d^2
		%     - AF ------
		%           dydz
		%         d^2           d^2         d^2         d^2         d^2         d^2    
		% H0 = B ------  + AE ------- + AD ------ + CH ------ + CG ------ + AF -----
		%         dx^2         dy^2         dz^2        dxdy        dxdz        dydz
		%       
		%
		% A = cos^2(theta); B = sin^2(theta); C = sin(2 theta); D = cos^2(phi);
		%
		% E = sin^2(phi); F = sin(2 phi); G = cos(phi); H = sin(phi).
		%
		% ---------------------------------------------------------------------------
		% Thus, one time step linear algerba form of TTI wave equation can be writen
		%		A1 U1 + A2 U2 + A3 U3 = B V2 + Q_t
		%		A1 V1 + B2 V2 + A3 V3 = C U2
		% 
		% ALL time step of it can be form as 
		%		A_all U + B V = Q; V = B_all^-1 C U;
		%	
		%			[ A1                           ] 
		%			| A2 A1                        | 
		%	A_all =	| A3 A2 A1                     | 
		%			|             . . .            | 
		%			[                      A3 A2 A1] 
		%
		%			[ A1                           ] 
		%			| B2 A1                        | 
		%	B_all =	| A3 B2 A1                     | 
		%			|             . . .            | 
		%			[                      A3 B2 A1] 
		%	
	elseif fwdpara.fmode == 11;
		% ------------------------------------------------------
		% do forward wavefield for T1 to Tmax
		% ------------------------------------------------------
		if fwdpara.fid && labindex == 1
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n'])
			fprintf(fwdpara.fid,[repmat(' ',1,10),'Solving forward modeling from T1 to Tmax for TTI\n'])
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n'])
		end


		% check point wavfield
		wavefield_chk_point = 0;
	
		spmd		
			%
			ntrace_loc 			= 0; % number of loc traces, Nr_loc * Ns_loc; 
			Num_of_shots_loc	= size(fwdpara_dis.slx_band,2);
			Num_of_time_step	= length(fwdpara.taxis);
		
		
			% locate receiver position (all shot has same receiver)
			if strcmp(fwdpara.rtype,'full')
				fwdpara_dis.rlx	= fwdpara.rlx;
				fwdpara_dis.rly	= fwdpara.rly;
				fwdpara_dis.rlz	= fwdpara.rlz;
				fwdpara_dis		= locate_rcv_position(fwdpara_dis);
				Pr				= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
												fwdpara_dis.rlz_loc,fwdpara_dis.rlx_loc,fwdpara_dis.rly_loc);
				rec_idx			= fwdpara_dis.rec_idx;
			elseif strcmp(fwdpara.rtype,'marine')
				rec_idx			= cell(Num_of_shots_loc,1);
			end

			% Num_of_rec_global	= size(fwdpara_dis.rlx,1);
			
			% create place for shot record
			output_data		= cell(1,Num_of_shots_loc);
		
			% create place for checkpoint wavefield 
			if fwdpara.chkp_save == 1
				wavefield_chk_point	= cell(Num_of_shots_loc,1);
			end
			[chk_point_idx,N_chkp]	= generate_check_point_idx(Num_of_time_step,fwdpara.chkp_space);
	

			% generate time stepping stencil matrix
			% A1 = [];A1_inv = []; A2 = [];A3 = [];A3_inv = [];	A1,A1_inv,A2,A3,A3_inv,dA1,dA2,dA3,B,B2,C
			[~,A1_inv,A2,A3,~,~,~,~,B,B2,C] = generate_time_stepping_stencil_matrix(fwdpara_dis,fwdpara);
			% combine some matrix to save time
			A2 = A1_inv * A2; 
			A3 = A1_inv * A3; 
			B	= A1_inv * B;
			B2	= A1_inv * B2;
			C	= A1_inv * C;
			
				
			
		
			% ---------------------------------- loop over shots ---------------------------------------
			for si  = 1:Num_of_shots_loc
			
				if fwdpara.fid && labindex == 1
					fprintf(fwdpara.fid,[repmat(' ',1,5),repmat('-',1,5),' Simulating shot ',num2str(si),' ',repmat('-',1,5),'\n']);
				end
			
				% locate src position and receiver position (each shot has diff receivers) 
				fwdpara_dis.slx		= vec(fwdpara_dis.slx_band(:,si));
				fwdpara_dis.sly		= vec(fwdpara_dis.sly_band(:,si));
				fwdpara_dis.slz		= vec(fwdpara_dis.slz_band(:,si));
				fwdpara_dis 		= locate_shot_position(fwdpara_dis);
				Ps					= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
												fwdpara_dis.slz_loc,fwdpara_dis.slx_loc,fwdpara_dis.sly_loc);
											
				if strcmp(fwdpara.rtype,'marine')
					fwdpara_dis.rlx	= fwdpara.rlx(:,fwdpara_dis.shots_id(si));
					fwdpara_dis.rly	= fwdpara.rly(:,fwdpara_dis.shots_id(si));
					fwdpara_dis.rlz	= fwdpara.rlz(:,fwdpara_dis.shots_id(si));
					fwdpara_dis 	= locate_rcv_position(fwdpara_dis);
					Pr				= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
													fwdpara_dis.rlz_loc,fwdpara_dis.rlx_loc,fwdpara_dis.rly_loc);
				end
			
				% generate source wavelet
				Q  					= prepare_source_wavelet(src,fwdpara,fwdpara_dis,si);
		
				% allocate place for loc shot record
				U_shot_record_temp	= zeros(length(fwdpara_dis.rec_idx),Num_of_time_step);

			
				% combine some matrix to save time
				Ps = A1_inv * Ps';

				% allocate wavefield snapshot
				U1	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				U2	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				U3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
			
				V1	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				V2	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				V3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
			
			
				
				% allocate space for checkpoint
				if fwdpara.chkp_save == 1 || fwdpara.chkp_save == 2
					chk_point_wavefeild_temp	= zeros(numel(U3).*2 ,N_chkp);
				end
				chk_point_count				= 1;			


				% ---------------------------- loop over time step ---------------------------------
				for ti = 1:Num_of_time_step

					%  forward modeling solve for A U  = Q
					%
					% [ A1                           ]  U1    Q1
					% | A2 A1                        |  U2    Q2
					% | A3 A2 A1                     |  U3  = Q3
					% |             . . .            |  U4    Q4
					% [                      A3 A2 A1]  U5    Q5
					if ti == 1
						U3(1:end)   =	B * V2(:)  +  Ps * Q(:,ti);
					elseif ti == 2
						U3(1:end)   =	B * V2(:)  +  Ps * Q(:,ti) - A2 * U2(:);
					else
						U3(1:end) 	=	B * V2(:)  +  Ps * Q(:,ti) - A2 * U2(:) - A3 * U1(:);
					end
					% U3(:) = A1_inv * U3(:);
			
				
				
					if labindex ==1 && mod(ti,50)==0 && fwdpara.fid
						fprintf(fwdpara.fid,[repmat(' ',1,10),'Time=',num2str(fwdpara.taxis(ti)),'s, Total=',num2str(fwdpara.taxis(end)),'s.\n'])

					end
			

					% ---------------- domain decomposition send data around -----------------
					if fwdpara_dis.num_model ~=1
						if fwdpara.ddcompz~=1
							if fwdpara_dis.ddcompz_id ~= 1
								labSend(U3(fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra,:,:),labindex-1)
							end
							if fwdpara_dis.ddcompz_id ~= fwdpara.ddcompz
								U3((end-fwdpara_dis.np_extra+1:end),:,:) = labReceive(labindex + 1);
								labSend(U3((end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra),:,:),labindex+1)
							end
							if fwdpara_dis.ddcompz_id ~= 1
								U3(1:fwdpara_dis.np_extra,:,:) = labReceive(labindex-1);
							end
						end
						labBarrier
						if fwdpara.ddcompx~=1
							if fwdpara_dis.ddcompx_id ~= 1
								labSend(U3(:,fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra,:),labindex-fwdpara.ddcompz)
							end
							if fwdpara_dis.ddcompx_id ~= fwdpara.ddcompx
								U3(:,(end-fwdpara_dis.np_extra+1:end),:) = labReceive(labindex + fwdpara.ddcompz);
								labSend(U3(:,(end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra),:),labindex+fwdpara.ddcompz)
							end
							if fwdpara_dis.ddcompx_id ~= 1
								U3(:,1:fwdpara_dis.np_extra,:) =  labReceive(labindex-fwdpara.ddcompz);
							end
						end
						labBarrier

						if fwdpara.ddcompy~=1
							if fwdpara_dis.ddcompy_id ~= 1
								labSend(U3(:,:,fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra),labindex-fwdpara.ddcompz*fwdpara.ddcompx)
							end
							if fwdpara_dis.ddcompy_id ~= fwdpara.ddcompy
								U3(:,:,(end-fwdpara_dis.np_extra+1:end)) = labReceive(labindex + fwdpara.ddcompz*fwdpara.ddcompx);
								labSend(U3(:,:,(end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra)),labindex+fwdpara.ddcompz*fwdpara.ddcompx)
							end
							if fwdpara_dis.ddcompy_id ~= 1
								U3(:,:,1:fwdpara_dis.np_extra) = labReceive(labindex-fwdpara.ddcompz*fwdpara.ddcompx);
							end
						end
						labBarrier
					end										
					% --------------- end of domain decomposition -----------------------

					if ti == 1
						
					elseif ti == 2
						V3(1:end)	= C * U2(:) - B2 * V2(:);
					else
						V3(1:end)	= C * U2(:) - B2 * V2(:) - A3 * V1(:);
					end
					% ---------------- domain decomposition send data around -----------------
					if fwdpara_dis.num_model ~=1
						if fwdpara.ddcompz~=1
							if fwdpara_dis.ddcompz_id ~= 1
								labSend(V3(fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra,:,:),labindex-1)
							end
							if fwdpara_dis.ddcompz_id ~= fwdpara.ddcompz
								V3((end-fwdpara_dis.np_extra+1:end),:,:) = labReceive(labindex + 1);
								labSend(V3((end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra),:,:),labindex+1)
							end
							if fwdpara_dis.ddcompz_id ~= 1
								V3(1:fwdpara_dis.np_extra,:,:) = labReceive(labindex-1);
							end
						end
						labBarrier
						if fwdpara.ddcompx~=1
							if fwdpara_dis.ddcompx_id ~= 1
								labSend(V3(:,fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra,:),labindex-fwdpara.ddcompz)
							end
							if fwdpara_dis.ddcompx_id ~= fwdpara.ddcompx
								V3(:,(end-fwdpara_dis.np_extra+1:end),:) = labReceive(labindex + fwdpara.ddcompz);
								labSend(V3(:,(end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra),:),labindex+fwdpara.ddcompz)
							end
							if fwdpara_dis.ddcompx_id ~= 1
								V3(:,1:fwdpara_dis.np_extra,:) =  labReceive(labindex-fwdpara.ddcompz);
							end
						end
						labBarrier

						if fwdpara.ddcompy~=1
							if fwdpara_dis.ddcompy_id ~= 1
								labSend(V3(:,:,fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra),labindex-fwdpara.ddcompz*fwdpara.ddcompx)
							end
							if fwdpara_dis.ddcompy_id ~= fwdpara.ddcompy
								V3(:,:,(end-fwdpara_dis.np_extra+1:end)) = labReceive(labindex + fwdpara.ddcompz*fwdpara.ddcompx);
								labSend(V3(:,:,(end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra)),labindex+fwdpara.ddcompz*fwdpara.ddcompx)
							end
							if fwdpara_dis.ddcompy_id ~= 1
								V3(:,:,1:fwdpara_dis.np_extra) = labReceive(labindex-fwdpara.ddcompz*fwdpara.ddcompx);
							end
						end
						labBarrier
					end										
					% --------------- end of domain decomposition -----------------------				 



				
 					% ---------------------------- display snapshot----------------------------
 					if fwdpara.disp_snapshot == 1 & mod(ti,20)==0
 						if length(fwdpara.mmy) ==1
 							% 2D ploting
 							figure(1);imagesc(U3);title([num2str(ti)]);colorbar;pause(.01)
							% figure(2);imagesc(V3);title([num2str(ti)]);colorbar;pause(.01)

 						else
 							% 3D ploting 
 							figure(1);sliceview(reshape(U3,length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc)));
 							title([num2str(ti)]);pause(.01)
 						end
 					end
				
				
 					% -------------- save check point per time step ---------------------
 					if chk_point_idx(ti) == 1
 						if fwdpara.chkp_save == 1
 							chk_point_wavefeild_temp(:,chk_point_count)	= [U3(:);V3(:)];
 							chk_point_count 	= chk_point_count + 1;
 						elseif fwdpara.chkp_save == 2
 							chk_point_wavefeild_temp(:,chk_point_count)	= [U3(:);V3(:)];
 							chk_point_count 	= chk_point_count + 1;
 						elseif fwdpara.chkp_save == 3
 							save_chkpoint_to_disk([U3(:);V3(:)],labindex,si,fwdpara.chkp_save,chk_point_count);
 							chk_point_count 	= chk_point_count + 1;
 						end
 					end
				

				
					% output receiver wavefield
					U_shot_record_temp(:,ti) = Pr * U3(:);

				
					U1 = U2;
					U2 = U3;
					U3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
			
					V1 = V2;
					V2 = V3;
					V3 = zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				end
				% save the wavefield at the last time step
				% U_shot_record_temp(:,ti)	= Pr * U2(:);
				output_data{si}				= U_shot_record_temp.';
				if strcmp(fwdpara.rtype,'marine')
					rec_idx{si}					= fwdpara_dis.rec_idx;
				end
				ntrace_loc = ntrace_loc + length(fwdpara_dis.rec_idx);
			
			
				% -------------- save check point per shot---------------------
				if fwdpara.chkp_save == 1
					wavefield_chk_point{si,1} = chk_point_wavefeild_temp;
				elseif fwdpara.chkp_save == 2
					save_chkpoint_to_disk(chk_point_wavefeild_temp,labindex,si,2,chk_point_count)
				end

			end
		
			% different workers may simulated different number of shots (e.g, worker 1 simulate 15 shots, worker 2 simulate 16 shots),
			% 'labBarrier' can not work while other workers are idle. 
			% this make sure all worker finish simulation. You can arrange workers to save one shot time. 
			if Num_of_shots_loc < fwdpara_dis.max_number_of_shot_all_worker
				for si = Num_of_shots_loc+1: fwdpara_dis.max_number_of_shot_all_worker
					for ti = 1:Num_of_time_step 
						labBarrier
						labBarrier
						labBarrier
					end
				end
			end
		
			% put all data together
			codistr			= codistributor1d(2,codistributor1d.unsetPartition,[1,numlabs]);
			ntrace_loc		= codistributed.build(ntrace_loc,codistr,'noCommunication'); 
		
		
		end

		[output_data,loc_rec_idx] = output_shot_record(output_data,rec_idx,fwdpara,fwdpara_dis,ntrace_loc,0);
	
	
	elseif fwdpara.fmode == 12; 
		if fwdpara.fid	&& labindex == 1	
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n'])
			fprintf(fwdpara.fid,[repmat(' ',1,10),'Solving forward modeling from Tmax to T1 for TTI.\n'])
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n'])
		end
		
	% ------------------------------------------------------
	% do forward wavefield from Tmax to T1 for TTI media
	% ------------------------------------------------------
	% NOTE: at least need u_n and u_n-1 and check point.
		spmd		
			%
			
			ntrace_loc 			= 0; % number of loc traces, Nr_loc * Ns_loc; 
			Num_of_shots_loc	= size(fwdpara_dis.slx_band,2);
			Num_of_time_step	= length(fwdpara.taxis);
			
			% locate receiver position (all shot has same receiver)
			if strcmp(fwdpara.rtype,'full')
				fwdpara_dis.rlx	= fwdpara.rlx;
				fwdpara_dis.rly	= fwdpara.rly;
				fwdpara_dis.rlz	= fwdpara.rlz;
				fwdpara_dis		= locate_rcv_position(fwdpara_dis);
				Pr				= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
												fwdpara_dis.rlz_loc,fwdpara_dis.rlx_loc,fwdpara_dis.rly_loc);
				rec_idx			= fwdpara_dis.rec_idx;
			elseif strcmp(fwdpara.rtype,'marine')
				rec_idx			= cell(Num_of_shots_loc,1);
			end

			% Num_of_rec_global	= size(fwdpara_dis.rlx,1);
		
			% create place for shot record
			output_data		= cell(1,Num_of_shots_loc);
		
			% create place for checkpoint wavefield 
			% if fwdpara.chkp_save == 1
			% 	wavefield_chk_point	= cell(Num_of_shots_loc,1);
			% end
			
			[chk_point_idx,N_chkp]	= generate_check_point_idx(Num_of_time_step,fwdpara.chkp_space);
		
			
			% generate time stepping stencil matrix
			% A1 = [];A1_inv = []; A2 = [];A3 = [];A3_inv = [];	
			[A1,~,A2,~,A3_inv,~,~,~,B,B2,C] = generate_time_stepping_stencil_matrix(fwdpara_dis,fwdpara);
			% [A1,~,A2,~,A3_inv] = generate_time_stepping_stencil_matrix(fwdpara_dis,fwdpara);

			% combine some matrix to save time
			A2	= A3_inv * A2; 
			A1	= A3_inv * A1; 
			B	= A3_inv * B;
			B2	= A3_inv * B2;
			C	= A3_inv * C;
			
			
			% ---------------------------- loop over shots ---------------------------------
			for si  = 1:Num_of_shots_loc
				if fwdpara.fid && labindex == 1
					fprintf(fwdpara.fid,[repmat(' ',1,5),repmat('-',1,5),' Simulating shot ',num2str(si),' ',repmat('-',1,5),'\n']);
				end
				% locate src position and receiver position (each shot has diff receivers) 
				fwdpara_dis.slx		= vec(fwdpara_dis.slx_band(:,si));
				fwdpara_dis.sly		= vec(fwdpara_dis.sly_band(:,si));
				fwdpara_dis.slz		= vec(fwdpara_dis.slz_band(:,si));
				fwdpara_dis 		= locate_shot_position(fwdpara_dis);
				Ps					= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
												fwdpara_dis.slz_loc,fwdpara_dis.slx_loc,fwdpara_dis.sly_loc);
			
				if strcmp(fwdpara.rtype,'marine')
					fwdpara_dis.rlx = fwdpara.rlx(:,fwdpara_dis.shots_id(si));
					fwdpara_dis.rly = fwdpara.rly(:,fwdpara_dis.shots_id(si));
					fwdpara_dis.rlz = fwdpara.rlz(:,fwdpara_dis.shots_id(si));
					fwdpara_dis 	= locate_rcv_position(fwdpara_dis);
					Pr				= Src_Rcv_interp(fwdpara_dis.mmz_loc,fwdpara_dis.mmx_loc,fwdpara_dis.mmy_loc, ...
													fwdpara_dis.rlz_loc,fwdpara_dis.rlx_loc,fwdpara_dis.rly_loc);
				end
			
				% generate source wavelet
				
				Q  = prepare_source_wavelet(src,fwdpara,fwdpara_dis,si);
		
				% allocate place for loc shot record
				U_shot_record_temp	= zeros(length(fwdpara_dis.rec_idx),Num_of_time_step);

				% combine some matrix to save time
				Ps = A3_inv * Ps';



				% allocate wavefield snapshot
				U1	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				U2	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				U3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
			
				V1	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				V2	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
				V3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));

			
				% allocate space for checkpoint
				if fwdpara.chkp_save == 1 
					chk_point_wavefeild_temp = wavefield_chk_point{si,1};
				elseif fwdpara.chkp_save == 2
					chk_point_wavefeild_temp = load_chkpoint_from_disk(labindex,si,fwdpara.chkp_save);
				end
				chk_point_count				 = N_chkp;			
		
		
				% ---------------------------- loop over time step ---------------------------------
				for ti = Num_of_time_step:-1:1

					%  forward modeling solve for A U  = Q
					%
					% [ A1                           ]  U1    Q1
					% | A2 A1                        |  U2    Q2
					% | A3 A2 A1                     |  U3  = Q3
					% |             . . .            |  U4    Q4
					% [                      A3 A2 A1]  U5    Q5
					
					
					% -------------- load check point per time step ---------------------
					if chk_point_idx(ti) == 1
						if fwdpara.chkp_save == 1
							U3(1:end) 				= reshape(chk_point_wavefeild_temp(1:numel(U3)    ,chk_point_count),size(U2));
							V3(1:end) 				= reshape(chk_point_wavefeild_temp(numel(U3)+1:end,chk_point_count),size(V2));
							chk_point_count = chk_point_count - 1;
						elseif fwdpara.chkp_save == 2
							U3(1:end)  				= reshape(chk_point_wavefeild_temp(1:numel(U3)    ,chk_point_count),size(U2));
							V3(1:end)  				= reshape(chk_point_wavefeild_temp(numel(U3)+1:end,chk_point_count),size(V2));
							chk_point_count 		= chk_point_count - 1;
						elseif fwdpara.chkp_save == 3
							chk_point_wavefeild_temp = load_chkpoint_from_disk(labindex,si,fwdpara.chkp_save,chk_point_count);
							U3(1:end)  				= reshape(chk_point_wavefeild_temp(1:numel(U3)    ),size(U2));
							V3(1:end)  				= reshape(chk_point_wavefeild_temp(numel(U3)+1:end),size(V2));
							chk_point_count			= chk_point_count - 1;
						end
						
					else
						U3(1:end)	=	B * V2(:)  +  Ps*Q(:,ti+2) - A2*U2(:) - A1 * U1(:);	
						% U3(1:end)	=	A3_inv * U3(:);
						% ---------------- domain decomposition send data around -----------------
						if fwdpara_dis.num_model ~=1
							if fwdpara.ddcompz~=1
								if fwdpara_dis.ddcompz_id ~= 1
									labSend(U3(fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra,:,:),labindex-1)
								end
								if fwdpara_dis.ddcompz_id ~= fwdpara.ddcompz
									U3((end-fwdpara_dis.np_extra+1:end),:,:) = labReceive(labindex + 1);
									labSend(U3((end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra),:,:),labindex+1)
								end
								if fwdpara_dis.ddcompz_id ~= 1
									U3(1:fwdpara_dis.np_extra,:,:) = labReceive(labindex-1);
								end
							end
							labBarrier
							if fwdpara.ddcompx~=1
								if fwdpara_dis.ddcompx_id ~= 1
									labSend(U3(:,fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra,:),labindex-fwdpara.ddcompz)
								end
								if fwdpara_dis.ddcompx_id ~= fwdpara.ddcompx
									U3(:,(end-fwdpara_dis.np_extra+1:end),:) = labReceive(labindex + fwdpara.ddcompz);
									labSend(U3(:,(end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra),:),labindex+fwdpara.ddcompz)
								end
								if fwdpara_dis.ddcompx_id ~= 1
									U3(:,1:fwdpara_dis.np_extra,:) =  labReceive(labindex-fwdpara.ddcompz);
								end
							end
							labBarrier

							if fwdpara.ddcompy~=1
								if fwdpara_dis.ddcompy_id ~= 1
									labSend(U3(:,:,fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra),labindex-fwdpara.ddcompz*fwdpara.ddcompx)
								end
								if fwdpara_dis.ddcompy_id ~= fwdpara.ddcompy
									U3(:,:,(end-fwdpara_dis.np_extra+1:end)) = labReceive(labindex + fwdpara.ddcompz*fwdpara.ddcompx);
									labSend(U3(:,:,(end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra)),labindex+fwdpara.ddcompz*fwdpara.ddcompx)
								end
								if fwdpara_dis.ddcompy_id ~= 1
									U3(:,:,1:fwdpara_dis.np_extra) = labReceive(labindex-fwdpara.ddcompz*fwdpara.ddcompx);
								end
							end
							labBarrier		
						end								
						% --------------- end of domain decomposition -----------------------
			
						V3(1:end)	= C * U2(:) - B2 * V2(:) - A1 * V1(:);
						% ---------------- domain decomposition send data around -----------------
						if fwdpara_dis.num_model ~=1
							if fwdpara.ddcompz~=1
								if fwdpara_dis.ddcompz_id ~= 1
									labSend(V3(fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra,:,:),labindex-1)
								end
								if fwdpara_dis.ddcompz_id ~= fwdpara.ddcompz
									V3((end-fwdpara_dis.np_extra+1:end),:,:) = labReceive(labindex + 1);
									labSend(V3((end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra),:,:),labindex+1)
								end
								if fwdpara_dis.ddcompz_id ~= 1
									V3(1:fwdpara_dis.np_extra,:,:) = labReceive(labindex-1);
								end
							end
							labBarrier
							if fwdpara.ddcompx~=1
								if fwdpara_dis.ddcompx_id ~= 1
									labSend(V3(:,fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra,:),labindex-fwdpara.ddcompz)
								end
								if fwdpara_dis.ddcompx_id ~= fwdpara.ddcompx
									V3(:,(end-fwdpara_dis.np_extra+1:end),:) = labReceive(labindex + fwdpara.ddcompz);
									labSend(V3(:,(end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra),:),labindex+fwdpara.ddcompz)
								end
								if fwdpara_dis.ddcompx_id ~= 1
									V3(:,1:fwdpara_dis.np_extra,:) =  labReceive(labindex-fwdpara.ddcompz);
								end
							end
							labBarrier

							if fwdpara.ddcompy~=1
								if fwdpara_dis.ddcompy_id ~= 1
									labSend(V3(:,:,fwdpara_dis.np_extra+1:2*fwdpara_dis.np_extra),labindex-fwdpara.ddcompz*fwdpara.ddcompx)
								end
								if fwdpara_dis.ddcompy_id ~= fwdpara.ddcompy
									V3(:,:,(end-fwdpara_dis.np_extra+1:end)) = labReceive(labindex + fwdpara.ddcompz*fwdpara.ddcompx);
									labSend(V3(:,:,(end-2*fwdpara_dis.np_extra+1):(end-fwdpara_dis.np_extra)),labindex+fwdpara.ddcompz*fwdpara.ddcompx)
								end
								if fwdpara_dis.ddcompy_id ~= 1
									V3(:,:,1:fwdpara_dis.np_extra) = labReceive(labindex-fwdpara.ddcompz*fwdpara.ddcompx);
								end
							end
							labBarrier
						end										
						% --------------- end of domain decomposition -----------------------				 


					end
			
					if labindex ==1 && mod(ti,50)==0 && fwdpara.fid;
						fprintf(fwdpara.fid,[repmat(' ',1,10),'Time=',num2str(fwdpara.taxis(ti)),'s, Total=',num2str(fwdpara.taxis(end)),'s.\n'])
					end

					% ---------------------------- display snapshot----------------------------
					if fwdpara.disp_snapshot == 1 & mod(ti,20)==0
						if length(fwdpara.mmy) ==1
							% 2D ploting
							figure(1);imagesc(U3);title([num2str(ti)]);colorbar;pause(.01)
							% figure(1);imagesc(V3);title([num2str(ti)]);colorbar;pause(.01)
						else
							% 3D ploting 
							figure(1);sliceview(reshape(U3,length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc)));
							title([num2str(ti)]);pause(.01)
						end
					end


					% shot record
					U_shot_record_temp(:,ti) = Pr * U3(:);
					
				
					% -------------- save check point per time step ---------------------
					% if chk_point_idx(ti) == 1
					% 	if fwdpara.chkp_save == 1
					% 		chk_point_wavefeild_temp(:,chk_point_count)	= U3(:);
					% 		chk_point_count = chk_point_count + 1;
					% 	elseif fwdpara.chkp_save == 2
					% 		chk_point_wavefeild_temp(:,chk_point_count)	= U3(:);
					% 		chk_point_count = chk_point_count + 1;
					% 	elseif fwdpara.chkp_save == 3
					% 		save_chkpoint_to_disk(U3,labindex,si,fwdpara.chkp_save,chk_point_count);
					% 		chk_point_count = chk_point_count + 1;
					% 	end
					% end
				
				
				
					U1	= U2;
					U2	= U3;
					U3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
					
					
					V1 = V2;
					V2 = V3;
					V3 = zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
					%
					% u1 = u2;
					% u2 = u3;
					% u3	= zeros(length(fwdpara_dis.mmz_loc),length(fwdpara_dis.mmx_loc),length(fwdpara_dis.mmy_loc));
					%
				
				end
				% save the wavefield at the last time step
				U_shot_record_temp(:,ti)		= Pr * U2(:);
				output_data{si}				= U_shot_record_temp.';
				if strcmp(fwdpara.rtype,'marine')
					rec_idx{si}					= fwdpara_dis.rec_idx;
				end
				ntrace_loc = ntrace_loc + length(fwdpara_dis.rec_idx);
			
			
				% -------------- save check point per shot---------------------
				%
				% if fwdpara.chkp_save == 1
				% 	wavefield_chk_point{si,1} = chk_point_wavefeild_temp;
				% elseif fwdpara.chkp_save == 2
				%
				% 	save_chkpoint_to_disk(chk_point_wavefeild_temp,labindex,si,2,chk_point_count)
				% end

			
			end
		
			% different workers may simulated different number of shots (e.g, worker 1 simulate 15 shots, worker 2 simulate 16 shots),
			% 'labBarrier' can not work while other workers are idle. 
			% this make sure all worker finish simulation. You can arrange workers to save one shot time. 
			if Num_of_shots_loc < fwdpara_dis.max_number_of_shot_all_worker
				for si = Num_of_shots_loc+1: fwdpara_dis.max_number_of_shot_all_worker
					for ti = 1:Num_of_time_step 
						labBarrier
						labBarrier
						labBarrier
					end
				end
			end
		
			% put all data together
			codistr			= codistributor1d(2,codistributor1d.unsetPartition,[1,numlabs]);
			ntrace_loc		= codistributed.build(ntrace_loc,codistr,'noCommunication'); 
		
		
		end

		[output_data,loc_rec_idx] = output_shot_record(output_data,rec_idx,fwdpara,fwdpara_dis,ntrace_loc,0);
	

	end % end of fmode
end % end of main function

% !@#$%^&*()*&^%$#@!@#$%^&*(&^%$#@!@#$%^&*(*&^%$#@!@#$%^&*(*&^%$#@!~!@#$%^&*()*&^%$#@!@#$%^&*()_)(*&^%$#@!@#$*(&)))
%                                      Main function is ended here      
%_)&^%$%^&*()_)(*&^%$#$%^&*()(*&^%$#@#$%^&*()_(*&^%$#@!@#$%^&*()_+_)(*&^%$#@!@#$%^&*()_(*&^%$#@!~!@#$%^&*()_)())))



%%=========================================================================================================
% =========================================================================================================
%                                                                                                         
%                                                                                                        
%                                             Subroutines                                                                          
%                                                                                                        
%                                                                                                        
% =========================================================================================================
% =========================================================================================================

function fwdpara = check_fwdpara_struct(fwdpara,model)
	% check forward structure file parameter 
	if  ~isfield(fwdpara,'slz')   | ~isfield(fwdpara,'slx')   | ~isfield(fwdpara,'rlx') | ~isfield(fwdpara,'rlz') | ...
		~isfield(fwdpara,'fcent') | ~isfield(fwdpara,'taxis') | ~isfield(fwdpara,'mmz') | ~isfield(fwdpara,'mmx') | ...
		~isfield(fwdpara,'fmode')
		error('missing modeling parameter struct fields, please check ');
	end
	if ~isfield(fwdpara,'stype'			), fwdpara.stype		 = 'seq_same'						;end
	if ~isfield(fwdpara,'rly'			), fwdpara.rly			 = ones(size(fwdpara.rlx))			;end
	if ~isfield(fwdpara,'sly'			), fwdpara.sly			 = ones(size(fwdpara.slx))			;end
	if ~isfield(fwdpara,'dt'			), fwdpara.dt			 = fwdpara.taxis(2)-fwdpara.taxis(1);end			
	if ~isfield(fwdpara,'abx'			), fwdpara.abx			 = 50								;end		
	if ~isfield(fwdpara,'abz'			), fwdpara.abz			 = 50								;end		
	if ~isfield(fwdpara,'mmy'			), fwdpara.mmy			 = 1								;end		
	if ~isfield(fwdpara,'v_up_type'		), fwdpara.v_up_type	 = 'slowness'						;end		
	if ~isfield(fwdpara,'chkp_space'	), fwdpara.chkp_space	 = 50								;end		
	if ~isfield(fwdpara,'chkp_save'		), fwdpara.chkp_save	 = 0								;end		
	if ~isfield(fwdpara,'disp_snapshot'	), fwdpara.disp_snapshot = 1								;end			
	if ~isfield(fwdpara,'fid'			), fwdpara.fid			 = 1								;end		
	if ~isfield(fwdpara,'free'			), fwdpara.free			 = 0								;end		
	if ~isfield(fwdpara,'wave_equ'		), fwdpara.wave_equ		 = 1								;end			
	if ~isfield(fwdpara,'space_order'	), fwdpara.space_order	 = 4								;end			
	if ~isfield(fwdpara,'free'			), fwdpara.free			 = 0								;end	
	if ~isfield(fwdpara,'ddcompx'		), fwdpara.ddcompx		 = 1								;end
	if ~isfield(fwdpara,'ddcompz'		), fwdpara.ddcompz		 = 1								;end
	if ~isfield(fwdpara,'ddcompy'		), fwdpara.ddcompy		 = 1								;end
	if ~isfield(fwdpara,'saved'			), fwdpara.saved		 = 1								;end			
	if ~isfield(fwdpara,'shot_id'		), fwdpara.shot_id		 = 1 								;end		
	if ~isfield(fwdpara,'aby'			), fwdpara.aby			 = length(fwdpara.mmy)				;end
	if ~isfield(fwdpara,'cut_model'		), fwdpara.cut_model	 = 0								;end
	if ~isfield(fwdpara,'cut_extra'		), fwdpara.cut_extra	 = 1000								;end			
	if ~isfield(fwdpara,'rtype'			)
		if size(fwdpara.rlx,2) == 1
			fwdpara.rtype = 'full';
		else
			fwdpara.rtype = 'marine';
		end
	end
	if fwdpara.wave_equ == 3 %TTI
		fwdpara.fmode = fwdpara.fmode + 10; % TTI using different mode.
		if  ~isfield(model,'epsilon') |  ~isfield(model,'delta') | ~isfield(model,'theta')
			error('fwdpara is missing TTI parameter struct fields, please check ');
		end
		if  ~isfield(model,'phi') & length(fwdpara.mmy) ~=1
			error('fwdpara is missing 3D TTI parameter struct fields (''phi''), please check ');
		end
	end	
end





function check_stability_condition(model,fwdpara)
	% function check_stability_condition(vel,fwdpara)
	% check stabilitiy condition.
	% -------- stability condition ---------------
	% dh < vmin/10/fcent,
	% Vmax * dt / dh  <= a. a=sqrt(3/8) for 2D, a=1/2 for 3D.
	% See Lines 1999, A recipe for stability of finite-difference wave-equation computations
	% --------------------------------------------
	dx = fwdpara.mmx(2) - fwdpara.mmx(1);
	dz = fwdpara.mmz(2) - fwdpara.mmz(1);
	dy = 1;
	if length(fwdpara.mmy)~=1
		dy = fwdpara.mmy(2) - fwdpara.mmy(1);
	end
	dt 		= fwdpara.taxis(2) - fwdpara.taxis(1);
	dh 		= max([dx,dy,dz]);
	vmin 	= min(model.vel);
	vmax 	= max(model.vel);

	if fwdpara.space_order == 2
		a		=sqrt(1/2);
		if length(fwdpara.mmy)~=1, a	= sqrt(1/3); end
	elseif fwdpara.space_order == 4 
		a		=sqrt(3/8);
		if length(fwdpara.mmy)~=1, a	= .5; end
	end

	dh_opt = vmin/10/fwdpara.fcent;
	fcent_opt =vmin/10/dh;
	if fwdpara.fid
		if dh > dh_opt;
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n']);
			fprintf(fwdpara.fid,['  Max space interval is ',num2str(dh),'m, which is larger than optinal ',num2str(dh_opt),'m.\n']);
			fprintf(fwdpara.fid,['  Modify it to optinal or lower the centrel frequency of your wavelet to ',num2str(fcent_opt),'Hz\n']);
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n']);
		end
		if dh < .5*dh_opt;
			fprintf(fwdpara.fid,[repmat('-',1,80),'.\n']);
			fprintf(fwdpara.fid,['  Max space interval is ',num2str(dh),'m, which is 2 time smaller than optinal ',num2str(dh_opt),'m.\n']);
			fprintf(fwdpara.fid,['  You could modify it to optinal for shorter modeling time or enlarger the centrel frequency of your\n']);
			fprintf(fwdpara.fid,['  wavelet to   ',num2str(fcent_opt),'Hz. All up to you...\n']);
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n']);
		end


		dt_opt = a * dh /vmax;
		if dt > dt_opt 
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n']);
			fprintf(fwdpara.fid,['  Time step length is ',num2str(dt),'s, which is larger than optinal ',num2str(dt_opt),'s\n']);
			fprintf(fwdpara.fid,['  Modify it to optinal. \n']);
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n']);

		end
		if dt < .5 * dt_opt
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n']);
			fprintf(fwdpara.fid,['  Time step length is ',num2str(dt),'s, which is too smaller than optinal ',num2str(dt_opt),'s\n']);
			fprintf(fwdpara.fid,['  You could enlarge it for shorter modeling time, unless you really want this. \n']);
			fprintf(fwdpara.fid,[repmat('-',1,80),'\n']);
		end
	end
end

%%
function  [fwdpara_dis] = setup_for_domain_decomp(model,fwdpara,input_data)
	% this function will setup parameters for demain decomperzation.

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
	Num_of_worker=parpool_size(); if ~Num_of_worker,Num_of_worker=1;end% number of labs
	Num_of_model_split = fwdpara.ddcompx * fwdpara.ddcompz * fwdpara.ddcompy;
	if  Num_of_model_split >  Num_of_worker
		error(['Number of works (',num2str( Num_of_worker),') can not be smaller ' ...
		'than number of model splits (',num2str( Num_of_worker),'); Need to be N ' ...
		'times larger, then you do not waste works.']);
	end
	 Num_of_shot_band = floor( Num_of_worker/ Num_of_model_split);
	 Num_of_all_shots = size(fwdpara.slx,2);
	 
	 shots_parition	  = round(linspace(0,Num_of_all_shots,Num_of_shot_band+1));
	 
	 
	 
	 nx		= length(fwdpara.mmx);
	 ny		= length(fwdpara.mmy);
	 nz		= length(fwdpara.mmz);
	 
	 ngx	= round(linspace(0,nx,fwdpara.ddcompx+1));
	 ngz	= round(linspace(0,nz,fwdpara.ddcompz+1));
 	 ngy	= round(linspace(0,ny,fwdpara.ddcompy+1));

	 
	mmx_cut_idx = 1:nx;
	mmy_cut_idx = 1:ny;
	mmz_cut_idx = 1:nz;  %% 

	 % create some datas for difference works
	fwdpara_dis = Composite();
	% vel_dis     = Composite();
	% den_dis     = Composite();
	% dvel_dis	= Composite();
	% mmx mmy mmz for each model split
	model.vel = reshape(model.vel,nz,nx,ny);
	if fwdpara.wave_equ == 2
		model.den = reshape(model.den,nz,nx,ny);
	elseif fwdpara.wave_equ == 3
		model.epsilon = reshape(model.epsilon,nz,nx,ny);
		model.delta	= reshape(model.delta,nz,nx,ny);
		model.theta	= reshape(model.theta,nz,nx,ny);
		
		% epsilon_dis	=	Composite();
		% delta_dis	=	Composite();
		% theta_dis	=	Composite();
		if length(fwdpara.mmy) ~= 1
			fwdpara.phi	=	reshape(fwdpara.phi,nz,nx,ny);
			phi_dis			=	Composite();
		end
	end

	if fwdpara.fmode   	== 3
		input_data	= reshape(input_data,nz,nx,ny);
	end
	
	for labi = 1:Num_of_worker
		model_parition_index	= mod(labi-1,Num_of_model_split)+1;
		[ddcompz_id,ddcompx_id,ddcompy_id]	= ind2sub([fwdpara.ddcompz,fwdpara.ddcompx,fwdpara.ddcompy],model_parition_index);
		shot_band_index			= ceil(labi/Num_of_model_split);
		num_of_extra_point		= max([2,fwdpara.space_order/2]); % extra points needs for domain decompzation

		% find x,y,z grid coordinate of model that on this worker.
		% find src locations on this work
		% We split all shots in to difference worker bands.
		% this is subset of the shots location for the worker band that this worker in
		shot_id		= shots_parition(shot_band_index)+1:shots_parition(shot_band_index+1);
		slx_dis		= fwdpara.slx(:,shot_id);
		sly_dis     = fwdpara.sly(:,shot_id);
		slz_dis     = fwdpara.slz(:,shot_id);
	
		%local receiver position 
		rlx_dis		= fwdpara.rlx(:,shot_id);
		rly_dis     = fwdpara.rly(:,shot_id);
		rlz_dis     = fwdpara.rlz(:,shot_id);
		
		

		
		if fwdpara.cut_model
			
			min_x   = min([slx_dis(:);rlx_dis(:)])-fwdpara.cut_extra;
			max_x   = max([slx_dis(:);rlx_dis(:)])+fwdpara.cut_extra;
			min_y   = min([sly_dis(:);rly_dis(:)])-fwdpara.cut_extra;
			max_y   = max([sly_dis(:);rly_dis(:)])+fwdpara.cut_extra;
			
			mmx_cut_idx = find(fwdpara.mmx>=min_x & fwdpara.mmx<=max_x);
			mmy_cut_idx = find(fwdpara.mmy>=min_y & fwdpara.mmy<=max_y);
			mmz_cut_idx = 1:length(fwdpara.mmz);  %% 

			
	   		ngx	= round(linspace(mmx_cut_idx(1)-1,mmx_cut_idx(end),fwdpara.ddcompx+1));
	   		ngy	= round(linspace(mmy_cut_idx(1)-1,mmy_cut_idx(end),fwdpara.ddcompy+1));
	    	ngz	= round(linspace(mmz_cut_idx(1)-1,mmz_cut_idx(end),fwdpara.ddcompz+1));

		end
		
		
		if fwdpara.ddcompx == 1
			mmx_loc = fwdpara.mmx(ngx(1)+1:ngx(2));
			u_idx   = 1:length(mmx_loc);
			x_id	= ngx(1)+1:ngx(2);
		else
			switch ddcompx_id
				case 1
					x_id	= ngx(ddcompx_id)+1:ngx(ddcompx_id+1)+num_of_extra_point;
					mmx_loc = fwdpara.mmx(x_id);
					u_idx   = 1:(length(mmx_loc)-num_of_extra_point);
				case fwdpara.ddcompx
					x_id	= ngx(ddcompx_id)-num_of_extra_point+1:ngx(ddcompx_id+1);
					mmx_loc = fwdpara.mmx(x_id);
					u_idx   = (num_of_extra_point+1):length(mmx_loc);
				otherwise
					x_id	= ngx(ddcompx_id)-num_of_extra_point+1:ngx(ddcompx_id+1)+num_of_extra_point;
					mmx_loc = fwdpara.mmx(x_id);
					u_idx   = num_of_extra_point+1:(length(mmx_loc)-num_of_extra_point);
			end
		end
		if fwdpara.ddcompz == 1
			mmz_loc = fwdpara.mmz(ngz(1)+1:ngz(2));
			u_idz   = 1:length(mmz_loc);
			z_id	= ngz(1)+1:ngz(2);
		else
			switch ddcompz_id
				case 1
					z_id	= ngz(ddcompz_id)+1:ngz(ddcompz_id+1)+num_of_extra_point;
					mmz_loc = fwdpara.mmz(z_id);
					u_idz   = 1:(length(mmz_loc)-num_of_extra_point);
				case fwdpara.ddcompz
					z_id	= ngz(ddcompz_id)-num_of_extra_point+1:ngz(ddcompz_id+1);
					mmz_loc = fwdpara.mmz(z_id);
					u_idz   = (num_of_extra_point+1):length(mmz_loc);
				otherwise
					z_id	= ngz(ddcompz_id)-num_of_extra_point+1:ngz(ddcompz_id+1)+num_of_extra_point;
					mmz_loc = fwdpara.mmz(z_id);
					u_idz   = (num_of_extra_point+1):(length(mmz_loc)-num_of_extra_point);
			end	
		end
		if fwdpara.ddcompy == 1
			mmy_loc  = fwdpara.mmy(ngy(1)+1:ngy(2));
			u_idy   = 1:length(mmy_loc);
			y_id	= ngy(1)+1:ngy(2);
		else
			switch ddcompy_id
				case 1
					y_id	= ngy(ddcompy_id)+1:ngy(ddcompy_id+1)+num_of_extra_point;
					mmy_loc = fwdpara.mmy(y_id);
					u_idy   = 1:(length(mmy_loc)-num_of_extra_point);
				case fwdpara.ddcompy
					y_id	= ngy(ddcompy_id)-num_of_extra_point+1:ngy(ddcompy_id+1);
					mmy_loc = fwdpara.mmy(y_id);
					u_idy   = (num_of_extra_point+1):length(mmy_loc);
				otherwise
					y_id	= ngy(ddcompy_id)-num_of_extra_point+1:ngy(ddcompy_id+1)+num_of_extra_point;
					mmy_loc = fwdpara.mmy(y_id);
					u_idy   = (num_of_extra_point+1):(length(mmy_loc)-num_of_extra_point);
			end
		end


		% find rec locations on this work
		
		
		% find model partition
		% x_id = ngx(ddcompx_id)+1:ngx(ddcompx_id+1);
		% z_id = ngz(ddcompz_id)+1:ngz(ddcompz_id+1);
		% y_id = ngy(ddcompy_id)+1:ngy(ddcompy_id+1);
	

		% find rec location on this work
		fwdpara_temp		= struct('vel_dis' ,model.vel(z_id,x_id,y_id)		...
									,'x_id_loc',x_id					...	
									,'y_id_loc',y_id					...
									,'z_id_loc',z_id					...
									,'u_idx',u_idx						...
									,'u_idy',u_idy						...
									,'u_idz',u_idz						...
									,'np_extra',num_of_extra_point 		...
									,'mmx_cut_idx',mmx_cut_idx			...
									,'mmy_cut_idx',mmy_cut_idx			...
									,'mmz_cut_idx',mmz_cut_idx			...
									,'mmx_loc',mmx_loc					...
									,'mmy_loc',mmy_loc					...
									,'mmz_loc',mmz_loc					...
									,'num_model',Num_of_model_split     ...
									,'ddcompx_id', ddcompx_id			...
									,'ddcompy_id', ddcompy_id			...
									,'ddcompz_id', ddcompz_id			...
									,'slx_band',slx_dis 				...
									,'sly_band',sly_dis 				...
									,'slz_band',slz_dis					...
									,'shots_id',shot_id					...
									,'max_number_of_shot_all_worker', max(diff(shots_parition)));



		% fwdpara_temp.vel_dis{labi} 		= vel(z_id,x_id,y_id);
		% fwdpara_temp.den_dis{labi}		= 1;
		if fwdpara.wave_equ == 2
			fwdpara_temp.den_dis	= model.den(z_id,x_id,y_id);
		elseif fwdpara.wave_equ == 3
			fwdpara_temp.epsilon_dis	=	model.epsilon(z_id,x_id,y_id);
			fwdpara_temp.delta_dis		=	model.delta(z_id,x_id,y_id);
			fwdpara_temp.theta_dis		=	model.theta(z_id,x_id,y_id);
			if length(fwdpara.mmy) ~= 1
				fwdpara_temp.phi_dis	=	model.phi(z_id,x_id,y_id);
			end
		end
		
		if fwdpara.fmode   	== 3
			fwdpara_temp.dvel_dis		= input_data(z_id,x_id,y_id);
		end

		
		fwdpara_dis{labi} = fwdpara_temp;
	end
	

end
	
	
function fwdpara_dis = locate_rcv_position(fwdpara_dis)
	% 
	ny = length(fwdpara_dis.mmy_loc);
	if ny == 1
		rec_idx  = find(fwdpara_dis.rlx>=fwdpara_dis.mmx_loc(fwdpara_dis.np_extra) & fwdpara_dis.rlx<fwdpara_dis.mmx_loc(end-fwdpara_dis.np_extra) & ...
						fwdpara_dis.rlz>=fwdpara_dis.mmz_loc(fwdpara_dis.np_extra) & fwdpara_dis.rlz<fwdpara_dis.mmz_loc(end-fwdpara_dis.np_extra));
						
		fwdpara_dis.rlx_loc =  fwdpara_dis.rlx(rec_idx);
		fwdpara_dis.rly_loc =  1;
		fwdpara_dis.rlz_loc =  fwdpara_dis.rlz(rec_idx);
		fwdpara_dis.rec_idx =  rec_idx;
	else
		rec_idx  = find(fwdpara_dis.rlx>=fwdpara_dis.mmx_loc(fwdpara_dis.np_extra) & fwdpara_dis.rlx<fwdpara_dis.mmx_loc(end-fwdpara_dis.np_extra) & ...
						fwdpara_dis.rly>=fwdpara_dis.mmy_loc(fwdpara_dis.np_extra) & fwdpara_dis.rly<fwdpara_dis.mmy_loc(end-fwdpara_dis.np_extra) & ...
						fwdpara_dis.rlz>=fwdpara_dis.mmz_loc(fwdpara_dis.np_extra) & fwdpara_dis.rlz<fwdpara_dis.mmz_loc(end-fwdpara_dis.np_extra));
		fwdpara_dis.rlx_loc =  fwdpara_dis.rlx(rec_idx);
		fwdpara_dis.rly_loc =  fwdpara_dis.rly(rec_idx);
		fwdpara_dis.rlz_loc =  fwdpara_dis.rlz(rec_idx);
		fwdpara_dis.rec_idx =  rec_idx;
	end	
end


function fwdpara_dis = locate_shot_position(fwdpara_dis)
	
	ny = length(fwdpara_dis.mmy_loc);
	
	if ny == 1
		src_idx  = find(fwdpara_dis.slx>fwdpara_dis.mmx_loc(1) & fwdpara_dis.slx<fwdpara_dis.mmx_loc(end) & ...
						fwdpara_dis.slz>fwdpara_dis.mmz_loc(1) & fwdpara_dis.slz<fwdpara_dis.mmz_loc(end));
		fwdpara_dis.slx_loc =  fwdpara_dis.slx(src_idx);
		fwdpara_dis.sly_loc =  fwdpara_dis.sly(src_idx);
		fwdpara_dis.slz_loc =  fwdpara_dis.slz(src_idx);		
	else
		src_idx  = find(fwdpara_dis.slx>fwdpara_dis.mmx_loc(1) & fwdpara_dis.slx<fwdpara_dis.mmx_loc(end) & ...
						fwdpara_dis.sly>fwdpara_dis.mmy_loc(1) & fwdpara_dis.sly<fwdpara_dis.mmy_loc(end) & ...
						fwdpara_dis.slz>fwdpara_dis.mmz_loc(1) & fwdpara_dis.slz<fwdpara_dis.mmz_loc(end));
		fwdpara_dis.slx_loc =  fwdpara_dis.slx(src_idx);
		fwdpara_dis.sly_loc =  fwdpara_dis.sly(src_idx);
		fwdpara_dis.slz_loc =  fwdpara_dis.slz(src_idx);
	end
	fwdpara_dis.src_idx = src_idx;
end

function fwdpara_dis = locate_rcv_position_adj(fwdpara_dis)
	% 
	ny = length(fwdpara_dis.mmy_loc);
	if ny == 1
		rec_idx  = find(fwdpara_dis.rlx>=fwdpara_dis.mmx_loc(fwdpara_dis.np_extra) & fwdpara_dis.rlx<fwdpara_dis.mmx_loc(end-fwdpara_dis.np_extra) & ...
						fwdpara_dis.rlz>=fwdpara_dis.mmz_loc(fwdpara_dis.np_extra) & fwdpara_dis.rlz<fwdpara_dis.mmz_loc(end-fwdpara_dis.np_extra));
						
		fwdpara_dis.rlx_loc =  fwdpara_dis.rlx(rec_idx);
		fwdpara_dis.rly_loc =  fwdpara_dis.rly(rec_idx);
		fwdpara_dis.rlz_loc =  fwdpara_dis.rlz(rec_idx);
		fwdpara_dis.rec_idx =  rec_idx;
	else
		rec_idx  = find(fwdpara_dis.rlx>=fwdpara_dis.mmx_loc(fwdpara_dis.np_extra) & fwdpara_dis.rlx<fwdpara_dis.mmx_loc(end-fwdpara_dis.np_extra) & ...
						fwdpara_dis.rly>=fwdpara_dis.mmy_loc(fwdpara_dis.np_extra) & fwdpara_dis.rly<fwdpara_dis.mmy_loc(end-fwdpara_dis.np_extra) & ...
						fwdpara_dis.rlz>=fwdpara_dis.mmz_loc(fwdpara_dis.np_extra) & fwdpara_dis.rlz<fwdpara_dis.mmz_loc(end-fwdpara_dis.np_extra));
		fwdpara_dis.rlx_loc =  fwdpara_dis.rlx(rec_idx);
		fwdpara_dis.rly_loc =  fwdpara_dis.rly(rec_idx);
		fwdpara_dis.rlz_loc =  fwdpara_dis.rlz(rec_idx);
		fwdpara_dis.rec_idx =  rec_idx;
	end	
end


function fwdpara_dis = locate_shot_position_adj(fwdpara_dis)
	
	ny = length(fwdpara_dis.mmy_loc);
	
	if ny == 1
		src_idx  = find(fwdpara_dis.slx>=fwdpara_dis.mmx_loc(fwdpara_dis.np_extra) & fwdpara_dis.slx<fwdpara_dis.mmx_loc(end-fwdpara_dis.np_extra) & ...
						fwdpara_dis.slz>=fwdpara_dis.mmz_loc(fwdpara_dis.np_extra) & fwdpara_dis.slz<fwdpara_dis.mmz_loc(end-fwdpara_dis.np_extra));
		fwdpara_dis.slx_loc =  fwdpara_dis.slx(src_idx);
		fwdpara_dis.sly_loc =  fwdpara_dis.sly(src_idx);
		fwdpara_dis.slz_loc =  fwdpara_dis.slz(src_idx);		
	else
		src_idx  = find(fwdpara_dis.slx>=fwdpara_dis.mmx_loc(fwdpara_dis.np_extra) & fwdpara_dis.slx<fwdpara_dis.mmx_loc(end-fwdpara_dis.np_extra) & ...
						fwdpara_dis.sly>=fwdpara_dis.mmy_loc(fwdpara_dis.np_extra) & fwdpara_dis.sly<fwdpara_dis.mmy_loc(end-fwdpara_dis.np_extra) & ...
						fwdpara_dis.slz>=fwdpara_dis.mmz_loc(fwdpara_dis.np_extra) & fwdpara_dis.slz<fwdpara_dis.mmz_loc(end-fwdpara_dis.np_extra));
		fwdpara_dis.slx_loc =  fwdpara_dis.slx(src_idx);
		fwdpara_dis.sly_loc =  fwdpara_dis.sly(src_idx);
		fwdpara_dis.slz_loc =  fwdpara_dis.slz(src_idx);
	end
	fwdpara_dis.src_idx = src_idx;
end



%% function can generate source term, it is very important when you save data to disk
function Q  = prepare_source_wavelet(src,fwdpara,fwdpara_dis,si)
	% this function 
	% Q need to be a matrix (number of source points x number of time step)
	if fwdpara.fmode == 1 || fwdpara.fmode ==2 || fwdpara.fmode ==3 ||  fwdpara.fmode ==-3 || fwdpara.fmode ==11 || fwdpara.fmode ==12
		
		%			stype:	'seq_same', 'src' needs to be a vector (source wavelet,length = Nt),
		%								all shot has same source wavelet.
		%					'seq_diff', 'src' needs to be a matrix (size of it = Ns x Nt, Ns=No of 
		%														  total shot locations);
		%								each columne is a difference source wavelet
		%					'sim_same', 'src' needs to be a matrix (size of it = Ns x Nt )
		%								each slice along y is one sim. source.
		%					'sim_diff', 'src' needs to be a 3D cube (size of it = Ns x Nt x Nshot )
		%								each slice along y is one sim. source. each shot has
		%								different source locations (like back propergation of 
		%								marine acquisition). 
		%								BUT EACH SHOT NEEDS TO HAVE SAME # OF SOURCE LOCATIONS
		%								(put 0s if they are not).	
		switch fwdpara.stype
		
			case 'seq_same'
				Ns = length(fwdpara_dis.slx_loc);
				Q  = repmat(src(:)',Ns,1);
			case 'seq_diff'
				Ns = length(fwdpara_dis.slx_loc);
				Q  = repmat(src(si,:),Ns,1);

			case 'sim_same'
				Ns = length(fwdpara_dis.slx_loc);
				% Q  = src(fwdpara_dis.src_idx,:);
				
				Q = repmat(vec(src)',Ns,1);
				% Q  = repmat(src,1,Ns);
				% Q  = reshape(Q,size(src),Ns);
			case 'sim_diff'
				Ns = length(fwdpara_dis.slx_loc);
				Q  = src(fwdpara_dis.src_idx,:,si);
				% Q  = repmat(src(:,:,si),1,Ns); 
				% Q  = reshape(Q,size(src(:,:,si)),Ns);
			end
	
	elseif fwdpara.fmode == -1
		switch fwdpara.saved
			case 0
				Q = src{si};
				Q = Q.';
			case 1
				Q = getLocalPart(src);
				n = length(Q);
				nt = length(fwdpara.taxis);
				nr = n/nt;
				Q = reshape(Q,nt,nr).';
				Q = Q(fwdpara_dis.Nr_count:fwdpara_dis.Nr_count+length(fwdpara_dis.rec_idx)-1,:);
			end

	end



end


function Q  = prepare_adjoint_source(in,fwdpara,fwdpara_dis,si)

	switch fwdpara.saved
		case 0
			Q = in{si};
			Q = Q.';
		case 1
			Q = getLocalPart(in);
			n = length(Q);
			nt = length(fwdpara.taxis);
			nr = n/nt;
			Q = reshape(Q,nt,nr).';
			Q = Q(fwdpara_dis.Nr_count:fwdpara_dis.Nr_count+length(fwdpara_dis.rec_idx)-1,:);
	end
end
%
function [U,loc_rec_idx] = output_shot_record(U,rec_idx,fwdpara,fwdpara_dis,ntrace_loc,loc_rec_idx)
	% how to output shot record, into a long distribute array or save to disk shot by shot
	
	if fwdpara.saved == 0
		
		%% output U and rec_idx in composite 

	elseif fwdpara.saved == 1
		if fwdpara.fmode == 1 | fwdpara.fmode == 11
			nr = size(fwdpara.rlx,1);
		elseif fwdpara.fmode == 2  | fwdpara.fmode == 12
			nr = size(fwdpara.rlx,1);
		elseif fwdpara.fmode == 3
			nr = size(fwdpara.rlx,1);
		elseif fwdpara.fmode == -1
			nr = size(fwdpara.slx,1);
		end
		ns = size(fwdpara.slx,2);
		nt = length(fwdpara.taxis);
		
		spmd
			U		= cell2mat(U);
			codistr	= codistributor1d(1,[ntrace_loc.*nt],[ns*nr*nt,1]);
			U		= codistributed.build(U(:),codistr,'noCommunication');
			codistr1	= codistributor1d(1,[ntrace_loc],[ns*nr,1]);
			loc_rec_idx = codistributed.build(loc_rec_idx(:),codistr1,'noCommunication');
		end
	elseif isstr(fwdpara.saved)
		Num_of_worker		= length(U);
		Num_of_model_split	= fwdpara.ddcompx * fwdpara.ddcompz * fwdpara.ddcompy;
		Num_of_shot_band	= floor( Num_of_worker/ Num_of_model_split);
		Num_of_all_shots	= size(fwdpara.slx,2);
	 	shots_parition		= round(linspace(0,Num_of_all_shots,Num_of_shot_band+1));
	 	
		if ~isdir(['./',fwdpara.saved]),mkdir (['./',fwdpara.saved]);end
		shot_id_count		= 1;
		 for shot_band_id 	= 1:Num_of_shot_band
			 idx			= (shot_band_id-1) * Num_of_model_split+1:shot_band_id*Num_of_model_split;
			 ns_loc			= length(U{idx(1)});
			 
			 for m = 1:ns_loc % loc shot number
				 
				 shot = zeros(length(fwdpara.taxis),size(fwdpara.rlx,1));
				 
				 for n = idx
					 if strcmp(fwdpara.rtype,'marine')
						 rec_id = rec_idx{n};
						 rec_id = rec_id{m};
						 rlx 	= fwdpara.rlx{:,si};
						 rly 	= fwdpara.rly{:,si};
						 rlz 	= fwdpara.rlz{:,si};
					 elseif  strcmp(fwdpara.rtype,'full')
						 rec_id	= rec_idx{n};
						 rlx 	= fwdpara.rlx;
						 rly 	= fwdpara.rly;
						 rlz 	= fwdpara.rlz; 
					 end
					 shot_temp = U{n};
					 shot(rec_id,:)  = shot_temp{m};

				 end

				 shot_id_count	= fwdpara.shot_id - 1 + shot_id_count;
				 taxis			= fwdpara.taxis;
				 eval(['save ./',fwdpara.saved,'/shot',num2str(shot_id_count),' shot taxis rlx rlz rly'])
				 shot_id_count	= shot_id_count + 1;
			 end
		 end
	end
		
end

function 	[chk_point_idx,N_chkp]	= generate_check_point_idx(Num_of_time_step,chkp_space)
	
	idx					= sort([Num_of_time_step:-chkp_space:1,Num_of_time_step-1:-chkp_space:1],'descend');
	N_chkp				= length(idx);
	chk_point_idx		= spalloc(1,Num_of_time_step,length(idx));
	chk_point_idx(idx)	= 1;

	
end


%% save instanteous checkpoint wavefield to disk
function save_chkpoint_to_disk(U,labind,si,smode,chk_point_count);
	if nargin < 5, chk_point_count = 1;end
	
	
	switch smode
		case 2
			tmppath = [getenv('TMPDIR'),'/check_point_wavefield/'];
			if ~isdir(tmppath),mkdir(tmppath);end
			eval(['save ',tmppath,'labind',num2str(labind),'shotind', ...
					num2str(si),' U']);		
		case 3
			tmppath = [getenv('TMPDIR'),'/check_point_wavefield/',num2str(labind),'shotind',num2str(si)];
			if ~isdir(tmppath)
				mkdir(tmppath);
			end
			eval(['save ',tmppath,'/chkpoint',num2str(chk_point_count),' U']);
	end

end

%% load instanteous checkpoint wavefield from disk
function U = load_chkpoint_from_disk(labind,si,smode,chk_point_count);
	if nargin < 4, chk_point_count = 1;end
	
	switch smode
		case 2
			tmppath = [getenv('TMPDIR'),'/check_point_wavefield/'];
			eval(['load ',tmppath,'labind',num2str(labind),'shotind', ...
					num2str(si)]);	
		case 3
			tmppath = [getenv('TMPDIR'),'/check_point_wavefield/',num2str(labind),'shotind',num2str(si)];
			eval(['load ',tmppath,'/chkpoint',num2str(chk_point_count)]);
	end

end

function g_all = gather_gradient(g,fwdpara_dis,fwdpara)
	% gather composite gradient into a long vec	
	
	
	
	
	
	nx		=	length(fwdpara.mmx);
	ny		=	length(fwdpara.mmy);
	nz		=	length(fwdpara.mmz);
	% ngx		=   round(linspace(0,nx,fwdpara.ddcompx+1));
	% ngz		=	round(linspace(0,nz,fwdpara.ddcompz+1));
	% ngy		=	round(linspace(0,ny,fwdpara.ddcompy+1));
	% g_all	= 	zeros(nz*nx*ny,1);
	g_all	= 	zeros(nz,nx,ny);
	
	Num_of_worker		= length(g);
	% Num_of_model_split	= fwdpara.ddcompx * fwdpara.ddcompz * fwdpara.ddcompy;
	% Num_of_shot_band	= floor( Num_of_worker/ Num_of_model_split);
	for labi = 1:Num_of_worker
		% model_parition_index	= mod(labi-1,Num_of_model_split)+1;
		% [ddcompz_id,ddcompx_id,ddcompy_id]	= ind2sub([fwdpara.ddcompz,fwdpara.ddcompx,fwdpara.ddcompy],model_parition_index);
		% shot_band_index			= ceil(labi/Num_of_model_split);
		% num_of_extra_point		= fwdpara.space_order/2;; % extra points needs for domain decompzation


		% find model partition
		% x_id = ngx(ddcompx_id)+1:ngx(ddcompx_id+1);
		% z_id = ngz(ddcompz_id)+1:ngz(ddcompz_id+1);
		% y_id = ngy(ddcompy_id)+1:ngy(ddcompy_id+1);
	
		% [XX,ZZ,YY] = meshgrid(x_id,z_id,y_id);
		
		% Model_idx  = sub2ind([nz,nx,ny],ZZ(:),XX(:),YY(:));
		
		fwdpara_temp = fwdpara_dis{labi};
		% g_all(Model_idx)    = g_all(Model_idx)  + vec(g{labi});
		g_all(fwdpara_temp.z_id_loc(fwdpara_temp.u_idz),fwdpara_temp.x_id_loc(fwdpara_temp.u_idx),fwdpara_temp.y_id_loc(fwdpara_temp.u_idy)) = ...
		g_all(fwdpara_temp.z_id_loc(fwdpara_temp.u_idz),fwdpara_temp.x_id_loc(fwdpara_temp.u_idx),fwdpara_temp.y_id_loc(fwdpara_temp.u_idy)) + g{labi}; 
	
	end
	
	g_all = g_all(:);
	
end


%% This function can give you a stencil matrix for time stepping
%   A3 Ut-1 + A2 Ut + A1 Ut+1 = Qt 
function [A1,A1_inv,A2,A3,A3_inv,dA1,dA2,dA3,B,B2,C] = generate_time_stepping_stencil_matrix(fwdpara_dis,fwdpara)
	% code up difference mode. 
	% mode = 1, 4 th order finite difference with acoustic wave equation of constant density
		dA1 = 0;dA2 = 0;dA3 = 0;B=0;B2=0;C=0;

		nx = length(fwdpara_dis.mmx_cut_idx);
		ny = length(fwdpara_dis.mmy_cut_idx);
		nz = length(fwdpara_dis.mmz_cut_idx);
		
		% define PML parameter
		% create pml boundary tensor
		dampf       = 30; % damping factor
		powerf      = 2;
		yeta_x(fwdpara_dis.mmx_cut_idx)      = dampf.*[fliplr(linspace(0,1,fwdpara.abx)),zeros(1,nx-2*fwdpara.abx),linspace(0,1,fwdpara.abx)].^powerf;
		yeta_y(fwdpara_dis.mmy_cut_idx)      = 0;
		if ny ~= 1
			yeta_y(fwdpara_dis.mmy_cut_idx)	 = dampf.*[fliplr(linspace(0,1,fwdpara.aby)),zeros(1,ny-2*fwdpara.aby),linspace(0,1,fwdpara.aby)].^powerf;
		end
		if fwdpara.free == 0
			yeta_z(fwdpara_dis.mmz_cut_idx)	= dampf.*[fliplr(linspace(0,1,fwdpara.abz)),zeros(1,nz-2*fwdpara.abz),linspace(0,1,fwdpara.abz)].^powerf;
		else
			yeta_z(fwdpara_dis.mmz_cut_idx)	= dampf.*[zeros(1,nz-fwdpara.abz),linspace(0,1,fwdpara.abz)].^2;
		end
		
		yeta_x_loc 	= yeta_x(fwdpara_dis.x_id_loc);
		yeta_y_loc 	= yeta_y(fwdpara_dis.y_id_loc);
		yeta_z_loc 	= yeta_z(fwdpara_dis.z_id_loc);
		
		[xx,zz,yy]	= meshgrid(yeta_x_loc,yeta_z_loc,yeta_y_loc);
		yeta		= xx + zz + yy; clear xx yy zz
	
		
		% size of the local model and wavefield
		nx_loc_m 	= length(fwdpara_dis.u_idx);
		ny_loc_m 	= length(fwdpara_dis.u_idy);
		nz_loc_m 	= length(fwdpara_dis.u_idz);
		
		nx_loc_u 	= length(fwdpara_dis.mmx_loc);
		ny_loc_u 	= length(fwdpara_dis.mmy_loc);
		nz_loc_u 	= length(fwdpara_dis.mmz_loc);		
		n_loc_u  	= nx_loc_u * ny_loc_u * nz_loc_u;
		
		% find out wavefield simulating area inside the local wavefeild(exclude domain decomp)
		U_idx_x 	= 1:nx_loc_u;
		U_idx_y 	= 1:ny_loc_u;
		U_idx_z 	= 1:nz_loc_u;
	
		
		[xx,zz,yy] = meshgrid(fwdpara_dis.u_idx,fwdpara_dis.u_idz,fwdpara_dis.u_idy);
		U_index    = sub2ind([nz_loc_u,nx_loc_u,ny_loc_u],zz(:),xx(:),yy(:));
		clear xx yy zz

		switch fwdpara.wave_equ
			case 1
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
				%	A1 = -(1+yeta*dt)/(v^2*dt^2);
				%	A2 = L - (1/v^2)(-2/dt^2  + yeta^2);
				%	A3 = -(1-yeta*dt)/(v^2*dt^2);
				%
				% -----------------------------------------------
				A1 			= zeros(n_loc_u,1);
				A1(U_index)	= -(1+vec(yeta(U_index)).*fwdpara.dt)./(vec(fwdpara_dis.vel_dis(U_index)).^2*fwdpara.dt^2);
				A1			= spdiags(A1,0,n_loc_u,n_loc_u);
				
				%A1_inv		= zeros(n_loc_u,1);
				A1_inv		= -(fwdpara_dis.vel_dis(:).^2*fwdpara.dt^2)./(1+yeta(:).*fwdpara.dt);
				A1_inv		= spdiags(A1_inv,0,n_loc_u,n_loc_u);

				A3 			= zeros(n_loc_u,1);
				A3(U_index)	= -(1-vec(yeta(U_index)).*fwdpara.dt)./(vec(fwdpara_dis.vel_dis(U_index)).^2*fwdpara.dt^2);
				A3			= spdiags(A3,0,n_loc_u,n_loc_u);
				
				%A3_inv		=  zeros(n_loc_u,1);
				A3_inv		= -(fwdpara_dis.vel_dis(:).^2*fwdpara.dt^2)./(1-yeta(:).*fwdpara.dt);
				A3_inv		= spdiags(A3_inv,0,n_loc_u,n_loc_u);
				
				
				% ------------------- finite difference stencil coefficients ----------------------------
				switch fwdpara.space_order
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
				dx 									= fwdpara.mmx(2) - fwdpara.mmx(1);
				if ny~=1,dy							= fwdpara.mmy(2) - fwdpara.mmy(1);	end
				dz 									= fwdpara.mmz(2) - fwdpara.mmz(1);
				Lx									= zeros(nx_loc_u,length(stencil_ceo));
				if ny~=1,Ly							= zeros(ny_loc_u,length(stencil_ceo));end
 				Lz									= zeros(nz_loc_u,length(stencil_ceo));
				Lx(fwdpara_dis.u_idx,:)				= ones(nx_loc_m,1)*stencil_ceo./(dx^2);
				if ny~=1,Ly(fwdpara_dis.u_idy,:)	= ones(ny_loc_m,1)*stencil_ceo./(dy^2);end
				Lz(fwdpara_dis.u_idz,:)				= ones(nz_loc_m,1)*stencil_ceo./(dz^2);
				
				diag_U_ind_x						= zeros(nx_loc_u,1);
				diag_U_ind_x(fwdpara_dis.u_idx)		= 1;
				diag_U_ind_y						= zeros(ny_loc_u,1);
				diag_U_ind_y(fwdpara_dis.u_idy)		= 1;
				diag_U_ind_z						= zeros(nz_loc_u,1);
				diag_U_ind_z(fwdpara_dis.u_idz)		= 1;
								
								
				if fwdpara.ddcompx~=1 && fwdpara_dis.ddcompx_id ==1
					Lx			= spdiags(Lx,-G_np_extra:G_np_extra,nx_loc_u,nx_loc_u)';
				elseif fwdpara.ddcompx~=1 && fwdpara_dis.ddcompx_id ==fwdpara.ddcompx
					Lx			= spdiags(Lx,-G_np_extra:G_np_extra,nx_loc_u,nx_loc_u)';
				else
					Lx			= spdiags(Lx,-G_np_extra:G_np_extra,nx_loc_u,nx_loc_u)';
				end	
				Lx 		 		= kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),kron(Lx,spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
				
				
				if fwdpara.ddcompz~=1 && fwdpara_dis.ddcompz_id ==1
					Lz			= spdiags(Lz,-G_np_extra:G_np_extra,nz_loc_u,nz_loc_u)';
				elseif fwdpara.ddcompz~=1 && fwdpara_dis.ddcompz_id ==fwdpara.ddcompz
					Lz			= spdiags(Lz,-G_np_extra:G_np_extra,nz_loc_u,nz_loc_u)';
				else
					Lz			= spdiags(Lz,-G_np_extra:G_np_extra,nz_loc_u,nz_loc_u)';
				end	
				
				Lz 				= kron(kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u)),Lz);
				
				if ny ~= 1
					if fwdpara.ddcompy~=1 && fwdpara_dis.ddcompy_id ==1
						Ly		= spdiags(Ly,-G_np_extra:G_np_extra,ny_loc_u,ny_loc_u)';
					elseif fwdpara.ddcompy~=1 && fwdpara_dis.ddcompy_id ==fwdpara.ddcompy
						Ly		= spdiags(Ly,-G_np_extra:G_np_extra,ny_loc_u,ny_loc_u)';
					else
						Ly		= spdiags(Ly,-G_np_extra:G_np_extra,ny_loc_u,ny_loc_u)';
					end			
	
					Ly			= kron(Ly,kron(spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u),spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
				
				
					L			= Lx + Lz + Ly;clear Lx Ly Lz
				% L(U_index,:) 	= kron(Ly,kron(Lx,Lz));clear Lx Ly Lz
				elseif ny == 1
					L			= Lx + Lz;clear Lx Ly Lz
				end
				
				M				=  zeros(n_loc_u,1);
				M(U_index)		=  (-2/fwdpara.dt^2 + vec(yeta(U_index)).^2)./vec(fwdpara_dis.vel_dis(U_index)).^2;
				M				=  spdiags(M,0,n_loc_u,n_loc_u);

 
				 A2 			= L - M; clear L M
				 
				 % dA for jacobian.
				 if fwdpara.fmode == 3 || fwdpara.fmode == -3	 
					 switch lower(fwdpara.v_up_type)
					 	case 'slowness'
							dA1				= zeros(n_loc_u,1);
							dA1(U_index)	= -2*(1+vec(yeta(U_index)).*fwdpara.dt)./(vec(fwdpara_dis.vel_dis(U_index))*fwdpara.dt^2);
							dA1				= spdiags(dA1,0,n_loc_u,n_loc_u);
							
							dA2				= zeros(n_loc_u,1);
							dA2(U_index)	= -2*(-2/fwdpara.dt^2 + vec(yeta(U_index)).^2)./vec(fwdpara_dis.vel_dis(U_index));
							dA2				= spdiags(dA2,0,n_loc_u,n_loc_u);

							dA3				= zeros(n_loc_u,1);
							dA3(U_index)	= -2*(1-vec(yeta(U_index)).*fwdpara.dt)./(vec(fwdpara_dis.vel_dis(U_index))*fwdpara.dt^2);
							dA3				= spdiags(dA3,0,n_loc_u,n_loc_u);
						case 'velocity'				 
							dA1				= zeros(n_loc_u,1);
							dA1(U_index)	= 2*(1+vec(yeta(U_index)).*fwdpara.dt)./(vec(fwdpara_dis.vel_dis(U_index)).^3*fwdpara.dt^2);
							dA1				= spdiags(dA1,0,n_loc_u,n_loc_u);
							
							dA2				= zeros(n_loc_u,1);
							dA2(U_index)	= 2*(-2/fwdpara.dt^2 + vec(yeta(U_index)).^2)./vec(fwdpara_dis.vel_dis(U_index).^3);
							dA2				= spdiags(dA2,0,n_loc_u,n_loc_u);

							dA3				= zeros(n_loc_u,1);
							dA3(U_index)	= 2*(1-vec(yeta(U_index)).*fwdpara.dt)./(vec(fwdpara_dis.vel_dis(U_index).^3)*fwdpara.dt^2);
							dA3				= spdiags(dA3,0,n_loc_u,n_loc_u);
						end
				 end

			case 2
				%% ===============================================
				%    acoustic wave equation with variable density
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
				switch fwdpara.space_order
					case 2
						% stencil_ceo = [-1/2 0 1/2]; % 2nd order
						stencil_ceo = [-1 1]; % 2nd order
						G_np_extra    = 1;
	
						% create gradient operator 'G'
						dx 									= fwdpara.mmx(2) - fwdpara.mmx(1);
						if ny~=1,dy							= fwdpara.mmy(2) - fwdpara.mmy(1);	end
						dz 									= fwdpara.mmz(2) - fwdpara.mmz(1);
						Lx									= zeros(nx_loc_u,length(stencil_ceo));
						if ny~=1,Ly							= zeros(ny_loc_u,length(stencil_ceo));end
		 				Lz									= zeros(nz_loc_u,length(stencil_ceo));
						Lx(fwdpara_dis.u_idx,:)				= ones(nx_loc_m,1)*stencil_ceo./(dx);
						if ny~=1,Ly(fwdpara_dis.u_idy,:)	= ones(ny_loc_m,1)*stencil_ceo./(dy);end
						Lz(fwdpara_dis.u_idz,:)				= ones(nz_loc_m,1)*stencil_ceo./(dz);
						
						diag_U_ind_x						= zeros(nx_loc_u,1);
						diag_U_ind_x(fwdpara_dis.u_idx)		= 1;
						diag_U_ind_y						= zeros(ny_loc_u,1);
						diag_U_ind_y(fwdpara_dis.u_idy)		= 1;
						diag_U_ind_z						= zeros(nz_loc_u,1);
						diag_U_ind_z(fwdpara_dis.u_idz)		= 1;
						
						G				= spdiags(vec(1./fwdpara_dis.den_dis(:)),0,n_loc_u,n_loc_u);
						% lx
						Lx1				= fliplr(Lx);
						Lx1				= spdiags(Lx1, -G_np_extra:0,nx_loc_u,nx_loc_u)';
						Lx1 		 	= kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),kron(Lx1,spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
						
						Gx				= .5 * ones(nx_loc_u,length(stencil_ceo));Gx(end,1) = 1;
						Gx				= spdiags(Gx, 0:G_np_extra,nx_loc_u,nx_loc_u);
						Gx 		 		= kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),kron(Gx,spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
						Gx				= spdiags(Gx*vec(1./fwdpara_dis.den_dis(:)),0,n_loc_u,n_loc_u);
						
						Lx1				= Gx * Lx1;
						
						Lx2				= fliplr(Lx);
						Lx2				= spdiags(Lx2, 0:G_np_extra,nx_loc_u,nx_loc_u)';
						Lx2 		 	= kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),kron(Lx2,spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
						
						Gx				= .5 * ones(nx_loc_u,length(stencil_ceo));Gx(1,2) = 1;
						Gx				= spdiags(Gx, -G_np_extra:0,nx_loc_u,nx_loc_u);
						Gx 		 		= kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),kron(Gx,spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
						Gx				= spdiags(Gx*vec(1./fwdpara_dis.den_dis(:)),0,n_loc_u,n_loc_u);
						Lx2				= Gx * Lx2;

						Lx				= (-Lx2  + Lx1)./dx;	clear Lx1 Lx2 Gx				
						
						%lz 
						
						Lz1				= fliplr(Lz);
						Lz1				= spdiags(Lz1, -G_np_extra:0,nz_loc_u,nz_loc_u)';
						Lz1 		 	= kron(kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u)),Lz1);
						
						Gz				= .5 * ones(nz_loc_u,length(stencil_ceo));Gz(end,1) = 1;
						Gz				= spdiags(Gz, 0:G_np_extra,nz_loc_u,nz_loc_u);
						Gz 		 		= kron(kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u)),Gz);
						Gz				= spdiags(Gz*vec(1./fwdpara_dis.den_dis(:)),0,n_loc_u,n_loc_u);
						
						Lz1				= Gz * Lz1;
						
						Lz2				= fliplr(Lz);
						Lz2				= spdiags(Lz2, 0:G_np_extra,nz_loc_u,nz_loc_u)';
						Lz2 		 	= kron(kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u)),Lz2);
						
						Gz				= .5 * ones(nz_loc_u,length(stencil_ceo));Gz(1,2) = 1;
						Gz				= spdiags(Gz, -G_np_extra:0,nz_loc_u,nz_loc_u);
						Gz 		 		= kron(kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u)),Gz);
						Gz				= spdiags(Gz*vec(1./fwdpara_dis.den_dis(:)),0,n_loc_u,n_loc_u);
						Lz2				= Gz * Lz2;

						Lz				= (- Lz2  + Lz1)./dz;	clear Lz1 Lz2 Gz								
						
						
						if ny ~= 1
						
							Ly1				= fliplr(Ly);
							Ly1				= spdiags(Ly1, -G_np_extra:0,ny_loc_u,ny_loc_u)';
							Ly1 		 	= kron(Ly1,kron(spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u),spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
						
							Gy				= .5 * ones(ny_loc_u,length(stencil_ceo));Gy(end,1) = 1;
							Gy				= spdiags(Gy, 0:G_np_extra,ny_loc_u,ny_loc_u);
							Gy 		 		= kron(Gy,kron(spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u),spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
							Gy				= spdiags(Gy*vec(1./fwdpara_dis.den_dis(:)),0,n_loc_u,n_loc_u);
						
							Ly1				= Gy * Ly1;
						
							Ly2				= fliplr(Ly);
							Ly2				= spdiags(Ly2, 0:G_np_extra,ny_loc_u,ny_loc_u)';
							Ly2 		 	= kron(Ly2,kron(spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u),spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
						
							Gy				= .5 * ones(ny_loc_u,length(stencil_ceo));Gy(1,2) = 1;
							Gy				= spdiags(Gy, -G_np_extra:0,ny_loc_u,ny_loc_u);
							Gy 		 		= kron(Gy,kron(spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u),spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
							Gy				= spdiags(Gy*vec(1./fwdpara_dis.den_dis(:)),0,n_loc_u,n_loc_u);
							Ly2				= Gy * Ly2;

							Ly				= (- Ly2  + Ly1)./dy;	clear Lz1 Lz2 Gy	
							L			= Lx + Ly + Lz;clear Lx Ly Lz G
						else
							
							L			= Lx + Lz;clear Lx Ly G
							
						end
						
				
						M			=  zeros(n_loc_u,1);
						M(U_index)	=  (-2/fwdpara.dt^2 + vec(yeta(U_index)).^2)./vec(fwdpara_dis.den_dis(U_index).*fwdpara_dis.vel_dis(U_index).^2);
						M			=  spdiags(M,0,n_loc_u,n_loc_u);

						A2 			= L - M;clear L M

				
						A1 			= zeros(n_loc_u,1);
						A1(U_index)	= -(1+vec(yeta(U_index)).*fwdpara.dt)./(vec(fwdpara_dis.den_dis(U_index).*fwdpara_dis.vel_dis(U_index).^2)*fwdpara.dt^2);
						A1			= spdiags(A1,0,n_loc_u,n_loc_u);
				
						%A1_inv		= zeros(n_loc_u,1);
						A1_inv		= -((fwdpara_dis.den_dis(:).*fwdpara_dis.vel_dis(:).^2)*fwdpara.dt^2)./(1+yeta(:).*fwdpara.dt);
						A1_inv		= spdiags(A1_inv,0,n_loc_u,n_loc_u);

						A3 			= zeros(n_loc_u,1);
						A3(U_index)	= -(1-vec(yeta(U_index)).*fwdpara.dt)./(vec(fwdpara_dis.den_dis(U_index).*fwdpara_dis.vel_dis(U_index).^2)*fwdpara.dt^2);
						A3			= spdiags(A3,0,n_loc_u,n_loc_u);
				
						%A3_inv		=  zeros(n_loc_u,1);
						A3_inv		=  -((fwdpara_dis.den_dis(:).*fwdpara_dis.vel_dis(:).^2)*fwdpara.dt^2)./(1-yeta(:).*fwdpara.dt);
						A3_inv		= spdiags(A3_inv,0,n_loc_u,n_loc_u);
				
	
						 % dA for jacobian.
						if fwdpara.fmode == 3 || fwdpara.fmode == -3	 
							 switch lower(fwdpara.v_up_type)
							 	case 'slowness'
									dA1				= zeros(n_loc_u,1);
									dA1(U_index)	= -2*(1+vec(yeta(U_index)).*fwdpara.dt)./(vec(fwdpara_dis.den_dis(U_index).*fwdpara_dis.vel_dis(U_index))*fwdpara.dt^2);
									dA1				= spdiags(dA1,0,n_loc_u,n_loc_u);
						
									dA2				= zeros(n_loc_u,1);
									dA2(U_index)	= -2*(-2/fwdpara.dt^2 + vec(yeta(U_index)).^2)./vec(fwdpara_dis.den_dis(U_index).*fwdpara_dis.vel_dis(U_index));
									dA2				= spdiags(dA2,0,n_loc_u,n_loc_u);

									dA3				= zeros(n_loc_u,1);
									dA3(U_index)	= -2*(1-vec(yeta(U_index)).*fwdpara.dt)./(vec(fwdpara_dis.den_dis(U_index).*fwdpara_dis.vel_dis(U_index))*fwdpara.dt^2);
									dA3				= spdiags(dA3,0,n_loc_u,n_loc_u);
								case 'velocity'				 
									dA1				= zeros(n_loc_u,1);
									dA1(U_index)	= 2*(1+vec(yeta(U_index)).*fwdpara.dt)./(vec(fwdpara_dis.den_dis(U_index).*fwdpara_dis.vel_dis(U_index).^3)*fwdpara.dt^2);
									dA1				= spdiags(dA1,0,n_loc_u,n_loc_u);
						
									dA2				= zeros(n_loc_u,1);
									dA2(U_index)	= 2*(-2/fwdpara.dt^2 + vec(yeta(U_index)).^2)./vec(fwdpara_dis.den_dis(U_index).*fwdpara_dis.vel_dis(U_index).^3);
									dA2				= spdiags(dA2,0,n_loc_u,n_loc_u);

									dA3				= zeros(n_loc_u,1);
									dA3(U_index)	= 2*(1-vec(yeta(U_index)).*fwdpara.dt)./(vec(fwdpara_dis.den_dis(U_index).*fwdpara_dis.vel_dis(U_index).^3)*fwdpara.dt^2);
									dA3				= spdiags(dA3,0,n_loc_u,n_loc_u);
								end
						 end
					 
					case 4
						stencil_ceo = [1, -8,0,8,-1]./12; % 4 th order
						G_np_extra    = 2;
					case 8
						stencil_ceo = [-1/60	3/20	-3/4	0	3/4	-3/20	1/60]; % 6th order
						G_np_extra    = 3;
					case 16
						stencil_ceo = [1/280	-4/105	1/5	-4/5	0	4/5	-1/5	4/105	-1/280];
						G_np_extra    = 4;
				end
				
			% --------------------------------------------------------------------------	
			case 3
				%% =========================================================================
				%			TTI acoustic wave equation
				%% =========================================================================
				%   1    d^2 U                            
				% ----- ------- - (1 + 2 sigma) H U - H0 U = (1+2 sigma) H V + Q
				%  c^2   d t^2                            
				%
				%   1    d^2 V
				% ----- ------- - 2(epsilo-sigma)H V  = 2 (epsilo-sigma) H U
				%  c^2   d t^2
				%
				% ---------------------------------- 2D case -------------------------------
				%
				%                      d^2                     d^2                    d^2
				%  H = cos^2 (theta) -------  + sin^2(theta) ------- - sin(2 theta) ------
				%                     d x^2                   d z^2                  dx dz
				%                      d^2                    d^2                     d^2
				%  H0 = sin^2(theta) -------  + cos^2(theta) -------- + sin(2 theta) ------
				%                     d x^2                   d z^2                   dx dz
				%
				% ---------------------------------- 3D case --------------------------------
				%         d^2                 d^2                d^2         d^2        d^2 
				% H = A ------- + (B + A D) ------- + (A - DB) ------- - CH ----- - CG -----
				%        d x^2               d y^2              d z^2        dxdy       dxdz
				%           d^2
				%     - AF ------
				%           dydz
				%         d^2           d^2         d^2         d^2         d^2         d^2    
				% H0 = B ------  + AE ------- + AD ------ + CH ------ + CG ------ + AF -----
				%         dx^2         dy^2         dz^2        dxdy        dxdz        dydz
				%       
				%
				% A = cos^2(theta); B = sin^2(theta); C = sin(2 theta); D = cos^2(phi);
				%
				% E = sin^2(phi); F = sin(2 phi); G = cos(phi); H = sin(phi).
				%
				% ---------------------------------------------------------------------------
				% Thus, one time step linear algerba form of TTI wave equation can be writen
				%		A1 U1 + A2 U2 + A3 U3 = B V2 + Q_t
				%		A1 V1 + B2 V2 + A3 V3 = C U2
				% 
				% ALL time step of it can be form as 
				%		A_all U + B V = Q; V = B_all^-1 C U;
				%	
				%			[ A1                           ] 
				%			| A2 A1                        | 
				%	A_all =	| A3 A2 A1                     | 
				%			|             . . .            | 
				%			[                      A3 A2 A1] 
				%
				%			[ A1                           ] 
				%			| B2 A1                        | 
				%	B_all =	| A3 B2 A1                     | 
				%			|             . . .            | 
				%			[                      A3 B2 A1] 
				%		
				
				
				A1 			= zeros(n_loc_u,1);
				A1(U_index)	= (1+vec(yeta(U_index)).*fwdpara.dt)./(vec(fwdpara_dis.vel_dis(U_index)).^2*fwdpara.dt^2);
				A1			= spdiags(A1,0,n_loc_u,n_loc_u);
			
				%A1_inv		= zeros(n_loc_u,1);
				A1_inv		= (fwdpara_dis.vel_dis(:).^2*fwdpara.dt^2)./(1+yeta(:).*fwdpara.dt);
				A1_inv		= spdiags(A1_inv,0,n_loc_u,n_loc_u);

				A3 			= zeros(n_loc_u,1);
				A3(U_index)	= (1-vec(yeta(U_index)).*fwdpara.dt)./(vec(fwdpara_dis.vel_dis(U_index)).^2*fwdpara.dt^2);
				A3			= spdiags(A3,0,n_loc_u,n_loc_u);
			
				%A3_inv		=  zeros(n_loc_u,1);
				A3_inv		=  (fwdpara_dis.vel_dis(:).^2*fwdpara.dt^2)./(1-yeta(:).*fwdpara.dt);
				A3_inv		= spdiags(A3_inv,0,n_loc_u,n_loc_u);

				M			=  zeros(n_loc_u,1);
				M(U_index)	=  (-2/fwdpara.dt^2 + vec(yeta(U_index)).^2)./vec(fwdpara_dis.vel_dis(U_index)).^2;
				M			=  spdiags(M,0,n_loc_u,n_loc_u);



				switch fwdpara.space_order
					case 2
						stencil_ceo_2 = [0 1 -2 1 0]; % 2nd order
						stencil_ceo_1 = [0 -1 0 1 0]./2;
						G_np_extra  = 2;
					case 4
						stencil_ceo_2 = [-1, 16,-30,16,-1]./12; % 4 th order
						stencil_ceo_1 = [ 1 -8 0 8 -1]./12;
						% stencil_ceo = [-0.1,1.4,-2.6,1.4,-0.1]; % Lluis's coefficient, contact him for more info
						%  "Guasch, Lluis" <l.guasch08@imperial.ac.uk>
						G_np_extra  = 2; 
					case 6
						stencil_ceo_2 = [ 1   -14   135  -245   135   -14     1]./90; % 6th order
						stencil_ceo_1 = [-1     9   -45     0    45    -9     1]./60;
						G_np_extra  = 3;
					case 8
						stencil_ceo_2 = [-1/560	 8/315	-1/5	 8/5	-205/72	8/5	-1/5	8/315	-1/560];
						stencil_ceo_1 = [ 1/280	-4/105	 1/5	-4/5		0	4/5	-1/5	4/105	-1/280];
						G_np_extra  = 4;
				end
				dx 									= fwdpara.mmx(2) - fwdpara.mmx(1);			
				dz 									= fwdpara.mmz(2) - fwdpara.mmz(1);
				Lx									= zeros(nx_loc_u,length(stencil_ceo_2));
 				Lz									= zeros(nz_loc_u,length(stencil_ceo_2));
				Lx(fwdpara_dis.u_idx,:)				= ones(nx_loc_m,1)*stencil_ceo_2./(dx^2);
				Lz(fwdpara_dis.u_idz,:)				= ones(nz_loc_m,1)*stencil_ceo_2./(dz^2);
				diag_U_ind_x						= zeros(nx_loc_u,1);
				diag_U_ind_x(fwdpara_dis.u_idx)		= 1;
				diag_U_ind_z						= zeros(nz_loc_u,1);
				diag_U_ind_z(fwdpara_dis.u_idz)		= 1;
				diag_U_ind_y						= zeros(ny_loc_u,1);
				diag_U_ind_y(fwdpara_dis.u_idy)		= 1;			
							
				if fwdpara.ddcompx~=1 && fwdpara_dis.ddcompx_id ==1
					Lx			= spdiags(Lx,-G_np_extra:G_np_extra,nx_loc_u,nx_loc_u)';
				elseif fwdpara.ddcompx~=1 && fwdpara_dis.ddcompx_id ==fwdpara.ddcompx
					Lx			= spdiags(Lx,-G_np_extra:G_np_extra,nx_loc_u,nx_loc_u)';
				else
					Lx			= spdiags(Lx,-G_np_extra:G_np_extra,nx_loc_u,nx_loc_u)';
				end	
				Lx 		 		= kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),kron(Lx,spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
			
				if fwdpara.ddcompz~=1 && fwdpara_dis.ddcompz_id ==1
					Lz			= spdiags(Lz,-G_np_extra:G_np_extra,nz_loc_u,nz_loc_u)';
				elseif fwdpara.ddcompz~=1 && fwdpara_dis.ddcompz_id ==fwdpara.ddcompz
					Lz			= spdiags(Lz,-G_np_extra:G_np_extra,nz_loc_u,nz_loc_u)';
				else
					Lz			= spdiags(Lz,-G_np_extra:G_np_extra,nz_loc_u,nz_loc_u)';
				end	
			
				Lz 				= kron(kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u)),Lz);
			



				switch ny
				case 1
					%--------------------- 2D case ---------------------------
					H	= spdiags(cos(fwdpara_dis.theta_dis(:)).^2,0,n_loc_u,n_loc_u) * Lx + spdiags(sin(fwdpara_dis.theta_dis(:)).^2,0,n_loc_u,n_loc_u) * Lz;
					
					H0	= spdiags(sin(fwdpara_dis.theta_dis(:)).^2,0,n_loc_u,n_loc_u) * Lx + spdiags(cos(fwdpara_dis.theta_dis(:)).^2,0,n_loc_u,n_loc_u) * Lz;
					clear Lx Lz
					
					Lx									= zeros(nx_loc_u,length(stencil_ceo_1));
	 				Lz									= zeros(nz_loc_u,length(stencil_ceo_1));
					Lx(fwdpara_dis.u_idx,:)				= ones(nx_loc_m,1)*stencil_ceo_1./(dx);
					Lz(fwdpara_dis.u_idz,:)				= ones(nz_loc_m,1)*stencil_ceo_1./(dz);
					if fwdpara.ddcompx~=1 && fwdpara_dis.ddcompx_id ==1
						Lx			= spdiags(Lx,-G_np_extra:G_np_extra,nx_loc_u,nx_loc_u)';
					elseif fwdpara.ddcompx~=1 && fwdpara_dis.ddcompx_id ==fwdpara.ddcompx
						Lx			= spdiags(Lx,-G_np_extra:G_np_extra,nx_loc_u,nx_loc_u)';
					else
						Lx			= spdiags(Lx,-G_np_extra:G_np_extra,nx_loc_u,nx_loc_u)';
					end	
					Lx 		 		= kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),kron(Lx,spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));
			
					if fwdpara.ddcompz~=1 && fwdpara_dis.ddcompz_id ==1
						Lz			= spdiags(Lz,-G_np_extra:G_np_extra,nz_loc_u,nz_loc_u)';
					elseif fwdpara.ddcompz~=1 && fwdpara_dis.ddcompz_id ==fwdpara.ddcompz
						Lz			= spdiags(Lz,-G_np_extra:G_np_extra,nz_loc_u,nz_loc_u)';
					else
						Lz			= spdiags(Lz,-G_np_extra:G_np_extra,nz_loc_u,nz_loc_u)';
					end	
			
					Lz 				= kron(kron(spdiags(diag_U_ind_y,0,ny_loc_u,ny_loc_u),spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u)),Lz);
					
					Lxz				= Lx * Lz;
					clear Lx Lz
					
					H	= H  - spdiags(sin(2.*fwdpara_dis.theta_dis(:)),0,n_loc_u,n_loc_u) * Lxz;
					H0	= H0 + spdiags(sin(2.*fwdpara_dis.theta_dis(:)),0,n_loc_u,n_loc_u) * Lxz;
					
		
					
					A2	= M -  spdiags( (1+2.*fwdpara_dis.delta_dis(:)),0,n_loc_u,n_loc_u)*H - H0;
					
					B	=  spdiags( (1+2.*fwdpara_dis.delta_dis(:)),0,n_loc_u,n_loc_u) * H;
					
					B2	= M - spdiags( 2.*(fwdpara_dis.epsilon_dis(:)-fwdpara_dis.delta_dis(:)),0,n_loc_u,n_loc_u) * H;
					
					C	= spdiags( 2.*(fwdpara_dis.epsilon_dis(:)-fwdpara_dis.delta_dis(:)),0,n_loc_u,n_loc_u) * H;
	
					
				otherwise 
					%--------------------- 3D case ---------------------------
					
					dy									= fwdpara.mmy(2) - fwdpara.mmy(1);
					Ly									= zeros(ny_loc_u,length(stencil_ceo_2));
					Ly(fwdpara_dis.u_idy,:)				= ones(ny_loc_m,1)*stencil_ceo_2./(dy^2);



					if fwdpara.ddcompy~=1 && fwdpara_dis.ddcompy_id ==1
					Ly		= spdiags(Ly,-G_np_extra:G_np_extra,ny_loc_u,ny_loc_u)';
					elseif fwdpara.ddcompy~=1 && fwdpara_dis.ddcompy_id ==fwdpara.ddcompy
					Ly		= spdiags(Ly,-G_np_extra:G_np_extra,ny_loc_u,ny_loc_u)';
					else
					Ly		= spdiags(Ly,-G_np_extra:G_np_extra,ny_loc_u,ny_loc_u)';
					end		
					Ly			= kron(Ly,kron(spdiags(diag_U_ind_x,0,nx_loc_u,nx_loc_u),spdiags(diag_U_ind_z,0,nz_loc_u,nz_loc_u)));



				end

				 
		end

end

