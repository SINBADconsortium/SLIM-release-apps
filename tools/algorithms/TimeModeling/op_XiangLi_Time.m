classdef op_XiangLi_Time < opSpot
	%% op_XiangLi_Time(model,src,fwdpara,wavefield_chk_point)
	% 2D/3D time domain forward modeling, 
	%
	%  
	% input arguements
	%	vel					: vectorized velocity model m/s, 
	%	den					: vectorized density mdoel g/cm^3
	%		NOTICE: vel should already include boundaries.
	%	src					: SEE 'fwdpara' setting (source location parameters);
	%	wavefield_chk_point	: wavefield checkpoint, To use this you need to do the modeling first 
	%						  (fwdpara.fmode=1), and save checkpoint. Then you can reconstruct wavefeild
	%						  from Tmax to T1.
	%	input_data			: input data for jacobian of FWI objective function;
	%						  for forward mode (fwdpara.fmode= 3), input data is model perturbation.
	%						  for adjoint mode (fwdpara.fmode=-3), input data is residual wavefield.
	%	fwdpara				: structure of modeling parameters
	%		wave equations
	%			wave_equ:	1, acoustic with only velocity
	%						2, acoustic with variable density 
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
	%		----------------------------------------------------------------------------
	%         parameters for Jacobian and its adjoint of FWI least-squares objective
	%		----------------------------------------------------------------------------
	%		Unit of updates
	%			v_up_type: 'slowness' or 'velocity'
	%
	%
	%
	% put a example here;
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
	% I am not responsiable for any problem that caused by using this code. You may not use it for commercial.
	% If you have any question, concern or sugguestions, please feel free to contact me.
	% ===========================================================================================================
	
	
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = protected)
   		funHandle = []; 
    end % Properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = op_XiangLi_Time(model,src,fwdpara,wavefield_chk_point)
		   
		   nr	  = size(fwdpara.rlx,1);
		   ns	  = size(fwdpara.slx,2);
		   nt	  = length(fwdpara.taxis);	   
		   n_data =  nr * ns * nt;
		   n_g    =  numel(model.vel);
		   
		   % Construct operator
         
  			fun = @(x,mode) op_XiangLi_Time_intrnl(model,src,fwdpara,wavefield_chk_point,x,mode);
            % Construct operator
            op = op@opSpot('XiangLi_Time_domain_Jacobian', n_data, n_g);
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



function y = op_XiangLi_Time_intrnl(model,src,fwdpara,wavefield_chk_point,x,mode)
	if (mode == 1)
		fwdpara.fmode = 3;
		y = XiangLi_Time(model,src,fwdpara,0,x);



	elseif (mode == 2)  
		fwdpara.fmode = -3;
		y = XiangLi_Time(model,src,fwdpara,0,x);

	end
end