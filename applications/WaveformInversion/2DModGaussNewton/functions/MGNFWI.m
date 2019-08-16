function [results,u,misfit,Mig] = MGNFWI(vel_initial,Dobs,wavelet,model,opts)
% This function allows us to process modified gauss-newton full-waveform 
% inversion, in which each gauss-newton subproblem is solved with the L1 solver
% (SPGL1) developed by Ewout van den Berg and Michael P. Friedlander. 
% 
% Example are publiced in the paper Xiang Li, Aleksandr Y. Aravkin, Tristan van
% Leeuwen, and Felix J. Herrmann, “Fast randomized full-waveform inversion with 
% compressive sensing”. 2011.
% 
% useage:
% 	[results,u,misfit,Mig] = MGNFWI(vel_initial,Dobs,wavelet,model,opts)
% 
% Input:
% 	vel_initial: starting model for FWI, 2D matrix; Vertical axis should be
% 	             depth axis; Surface on the top. (unit: meters/second)
% 	Dobs: obeservation data, 3D data cube; Vertical axis should be receiver axis
% 	      lateral axis should be source axis; Third axis should frequency axis
% 	wavelet: data struct with wavelet information. 
% 			wavelet.t0: starting time (second)
% 			wavelet.f0: central frequency (Hz)
% 			wavelet.faxis: frequencies (Hz)
% 	model: data struct with modeling parameters.
%		minf: minimal freqency, Hz
%		maxf: maximal frequency, Hz
%		nf  : number of frequencies in each frequency band
%		ol  : number of overlap frequencies for two adjacent bands
%		snf : number of frequencies for each GN updates 
%		gl  : mesh distance, Meters
%		xbound: x dimension boundary size, percentage
%		zbound: z dimension boundary size, percentage
%		nshot: number of shot positions
%		nsim : number of simultaneous shots used for each GN update
%		sdep : source depth, unit: meters
%		sp   : source positions in x dimension
%		nrec : number of receivers
%		rp   : receiver positions
%		rdep : receiver depth, unit: meters
%		water: estimated water depth, unit: number of grid
%		sourcet: source type. 1, simultaneous source; 2, sequential source.
%		offset: marine data offset mask. 1, with mask; 2, without mask.
%		maxoffst: maximul offset
%		minoffst: minimul offset
%		vmin: estimated maximal velocity For update projection
%		vmax: estimated minimal velocity
%		name: name of results saved on disk
%
% 	opts: data struct with optimization papermeters.
%		fwistep: 1, only do Forward modeling; 2, output misfit and Jacobi;3,do Fwi
%		renewal: raw different set of shots for each GN update. 1,yes; 2,no.
%		miter: number for frequency bands
%		mlinm: number for GN subproblems for each freq band
%		iterations: iteration number for L1 solver
%		verbosity: help spgl1 for more info. 0=quiet, 1=some output, 2=more output.
%		bpTol: help spgl1 for more info. Tolerance for identifying a basis pursuit solution.
%		optTol: help spgl1 for more info. Optimality tolerance (default is 1e-4).
%		decTol: help spgl1 for more info. Larger decTol means more frequent Newton updates
%		quitPareto: help spgl1 for more info. 1, quit; 2 not quit
%       sptrans: sparse transform. 'yes': curvelet domain. 'no': physical domain
%		ny: size of sparsfiy transform
%		nx: size of sparsfiy transform
%		scale: help curvelab, number of scales
%		angle: help curvelab, number of angles
%		finst: help curvelab, finst grid
%		dispresult: display the result interactively, 1 yes, 0 no.
% 			
% Author: Xiang Li
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: 02, 2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

tic 
%========================== initializing =============================

model.t0 = wavelet.t0;
model.f0 = wavelet.f0;
mod_def = Set_Paras( ...
'minf',0, ...                    
'maxf',30 , ...				     
'nf',10, ...                     
'ol',5, ...					
'snf',10, ...			
'xbound',.1, ...			   
'zbound',.1, ... 
'water',1, ...				
'sourcet',2, ...			
'offset',1, ...					
'maxoffst',3000, ...		   
'minoffst',100, ...			  
'vmin',1500, ...			    
'vmax',5000, ...		       
'name','FwiBGmodel');	               

model       = Set_Paras(mod_def,model);clear mod_def;
index       = find(wavelet.faxis>=model.minf & wavelet.faxis<=model.maxf);
wavelet.faxis = wavelet.faxis(index);
[model.nz,model.nx] = size(vel_initial);
model.nx_border   = floor(model.xbound*model.nx);
model.nz_border   = floor(model.zbound*model.nz);

% do some check
if ~isfield(model,'sp'),model.sp = round(linspace(1,model.nx,model.nshot));end
if length(model.rp) ~= model.nrec
	error('number of receivers positions should be the same as number of receivers')	
end

opts_def = Set_Paras( ...
'fwistep',3, ...                      % fwi step. 1, only do Forward modeling; 2, misfit and Jacobi;3,do Fwi
'renewal',1 , ...                     % draw different set of shots for each GN update. 1,yes; 2,no.
'miter',10, ...		                  % iteration number for frequency band.
'mlinm',10, ...			 	  		  % iteration number for relinearization for each freq band
'iterations',20, ...				  % iteration number for L1 solver ...
'totiter_freqband',200, ...           % total L1 solver iteration for each frequency band
'verbosity',1, ...					  % 0=quiet, 1=some output, 2=more output.
'bpTol',1e-3, ...				      % help spgl1. Tolerance for identifying a basis pursuit solution.
'optTol',1e-3, ...                    % help spgl1. Optimality tolerance (default is 1e-4).
'decTol',5e-2, ...                    % help spgl1. Larger decTol means more frequent Newton updates...
'quitPareto',0, ...                   % help spgl1. 1, quit; 2 not quit
'sptrans','yes', ...                  % sparse basis
'ny',model.nz, ...					  % size of the basis
'nx',model.nx, ...                    % size of the basis
'scale',4, ...	                      % help curvelab, number of scales
'angle',16, ...			              % help curvelab, number of angle
'finst',0, ...                        % help curvelab, finst grid ...
'dispresult',0); 
opts  = Set_Paras(opts_def,opts); clear opts_def;


% setup parameters for Frequency domain FD modeling
model.o = [0 0 0];
model.d = [model.gl model.gl 1];
model.n = [model.nz model.nx 1];        
model.nb = [model.nz_border model.nx_border 0];
model.zsrc = model.sdep; 
model.xsrc = (model.sp - 1) * model.gl ;
model.zrec = model.rdep; 
model.xrec = (model.rp - 1) * model.gl ;


% results and updates container
results     = vel_initial;
updates     = [];
resultinfo  = [];


% source type defination, simoultaneous or sequential 
[Q,simp,G] = source_fun(model.gl,model.nx,model.nz,model.nx_border, ...
model.nz_border,model.nshot,model.sp,model.sdep,model.sourcet,model.nsim);

% offset mask defination
if model.offset && model.sourcet == 2
	ofstind = offsetmask(model.nx,model.gl,model.sp(simp),model.minoffst,model.maxoffst);
	ofstR   = opColumnRestrict(1,model.nx*model.nsim,ofstind);
	ofstR   = oppKron2Lo(opDirac(model.snf),ofstR);
	clear ofstind
end


% squared slowness [s^2/km^2]
m_initial   = 1e6./vel_initial(:).^2;


%=========================== FWI loop start here ================================

%---------------------- outer loop over frequency bands -------------------------
iter = 1;
while length(wavelet.faxis)>=model.nf && wavelet.faxis(1)<model.maxf && iter<=opts.miter
% for iter = 1:opts.miter
	dmw_c = [];  % initial guess of update
	disp(['=== Now processing frequency band ',num2str(iter),'==='])
	disp(' ')
    freqband    = wavelet.faxis(1:model.nf);

%--------------------- GN subproblems for each frequency band -------------------
	linm = 1;
	tot_iter = 0;
	GN_iter = 0;
	while tot_iter < opts.totiter_freqband && linm <= opts.mlinm
		disp(['---calculating GN update ',num2str(linm),' for frequency band ',num2str(iter),'---'])
		disp(' ')
		
		idx   = randperm(model.nf);
		idx   = sort(idx(1:model.snf));
		model.freq    = freqband(idx);

		% load observation data
		index  = ((model.nf-model.ol)*(iter-1)+1):((model.nf-model.ol)*(iter-1)+model.nf);
		index  = index(idx);
		uh = vec(distributed(Dobs(:,:,index)));
		RM = oppKron2Lo(opDirac(model.snf), opKron(G,opDirac(model.nrec)));
		uh = RM * uh;

		% calculate jacobian and wavefield difference 
		[u,Mig]  = F(m_initial,Q,model); 
		du_nl    = uh - u;
		disp(['two-norm wavefield difference of current frequency band is: ',num2str(norm(du_nl))])

		% output misfit and Jacobian
		if opts.fwistep == 2
			misfit = norm(du_nl);
			return;
		end
		
		% apply offset mask
		if model.offset && model.sourcet == 2
			du_nl  = ofstR * du_nl;
			Mig    = ofstR * Mig;
		end
		
		% GN optimization
		[dmw_nl,info] = GN_up(Mig,du_nl,opts,dmw_c);
		resultinfo{end + 1} = info;
		dmw_nl      = reshape(real(dmw_nl),model.nz,model.nx);
		dmw_nl      = Smoothedge(dmw_nl,model.water,1);
		
		% do update
		GN_iter = GN_iter + info.iter;
		if GN_iter >= opts.iterations
			m_initial   = m_initial + dmw_nl(:);
			
			% do projection
			index       = find(m_initial < 1e6./model.vmax^2);
			m_initial(index) = 1e6./model.vmax^2;
			index       = find(m_initial > 1e6./model.vmin^2);
			m_initial(index) = 1e6./model.vmin^2;
			vel_initial = reshape(m_initial,model.nz,model.nx);
			dmw_c = [];
			GN_iter = 0;
		end
		
		% save results and updates
		results(:,:,end+1) = sqrt(1e6./vel_initial);
		updates(:,:,end+1) = dmw_nl;
		
		if opts.dispresult
			figure(1);imagesc(model.gl*(1:model.nz),model.gl*(1:model.nx),vel_initial);
			xlabel('x [m]');ylabel('z [m]');colorbar;title(['reverted result from frequency band',num2str(iter)]);
			pause(.1);
		end
		% define saving path
		path1 = which(mfilename);
		path2 = '/functions';
		ind   = strfind(path1,path2);
		path1 = path1(1:ind(end)-1);
		
		% save results to the disk
		eval(['save ',path1,'/results',model.name,' results updates resultinfo']);
		
		% do renewals, draw other shots for the next iteration
		if opts.renewal
			[Q,simp,G] = source_fun(model.gl,model.nx,model.nz,model.nx_border, ...
			model.nz_border,model.nshot,model.sp,model.sdep,model.sourcet,model.nsim);
			if model.offset && model.sourcet == 2
				ofstind = offsetmask(model.nx,model.gl,model.sp(simp),model.minoffst,model.maxoffst);
				ofstR   = opColumnRestrict(1,model.nx*model.nsim,ofstind);
				ofstR   = oppKron2Lo(opDirac(model.snf),ofstR);
				clear ofstind
			end
		end
		tot_iter = tot_iter + info.iter;
		linm = linm + 1;
%---- end of GN subloop ---
	end
%----- end of FWI loop ----
iter              = iter + 1;
ol                = model.nf - model.ol;
wavelet.faxis(1:ol) = [];	
end

time        = toc;
eval(['save ',path1,'/results',model.name,' results updates resultinfo']);