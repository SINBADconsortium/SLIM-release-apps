function [du,r]=Born(v,model1,source,din,mode,dens,ani)
% Born --Born modelling function, Linearized dat from model perturbation.
%
% Function that generates the data cube Time x Rec x Src from a source for
% a given model perturbation and velocity in forward mode and does RTM in
% reverse mode
%
% INPUTS
%
%   'model1' : stucture containing the model parameters
%       'v'  :  square slowness [s^2/m^2]
%     'dens' : density [g/cm^3]
%     'source' : source term  (1 x nt)
%     'anis' : anisotropy parameters (epsilon,delta,theta in a
%     structure)
%     'mode' : 1 for forward, -1 for adjoint
%     'din' : model perturbation [s^2/km^2] in forward mode, data residual
%     or data in reverse mode (needs to put two out put to say it is the data).
%     'source' : source term  (1 x nt)
%
% OUTPUT from multiplication
%
%   'du':   linearized data cube in forward mode, RTM image in reverse mode
%   'r' :   data residual if the input is the data when in reverse mode
%
%
% Author:Mathias Louboutin
%        Philipp Witte
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
% Date: July, 2015
nout=nargout;

%% Extend the model for PML's
if model1.freesurface
    nbz1=model1.space_order/2;
else
    nbz1=40;
end
% Redifine dimensions
model1.nb=[40 40 40];

if nargin>=6 && ~isempty(dens)
    if model1.n(3)~=1
        dens=kron(opExtension(model1.n(3),model1.nb(3)),opExtension(model1.n(2),model1.nb(2)),opExtension(model1.n(1),[nbz1 model1.nb(1)]))*dens(:);
    else
        dens=kron(opExtension(model1.n(2),model1.nb(2)),opExtension(model1.n(1),[nbz1 model1.nb(1)]))*dens(:);
    end
else
    dens=[];
end
if nargin==7 && ~isempty(ani)
    ani.epsilon=kron(opExtension(model1.n(2),model1.nb(2)),opExtension(model1.n(1),[nbz1 model1.nb(1)]))*ani.epsilon(:);
    ani.delta=kron(opExtension(model1.n(2),model1.nb(2)),opExtension(model1.n(1),[nbz1 model1.nb(1)]))*ani.delta(:);
    ani.theta=kron(opExtension(model1.n(2),model1.nb(2)),opExtension(model1.n(1),[nbz1 model1.nb(1)]))*ani.theta(:);
else
    ani=[];
end

if model1.n(3)~=1
    v=kron(opExtension(model1.n(3),model1.nb(3)),opExtension(model1.n(2),model1.nb(2)),opExtension(model1.n(1),[nbz1 model1.nb(1)]))*v(:);
    if mode==1
        din=kron(opExtension(model1.n(3),model1.nb(3)),opExtension(model1.n(2),model1.nb(2)),opExtension(model1.n(1),[nbz1 model1.nb(1)]))*din(:);
    end
    model1.o=[model1.o(1)-nbz1*model1.d(1) model1.o(2)-model1.nb(2)*model1.d(2) model1.o(3)-model1.nb(3)*model1.d(3)];
    model1.n=model1.n+[sum([nbz1 model1.nb(1)]) 2*model1.nb(2) 2*model1.nb(3)];
else
    v=kron(opExtension(model1.n(2),model1.nb(2)),opExtension(model1.n(1),[nbz1 model1.nb(1)]))*v(:);
    if mode==1
        din=kron(opExtension(model1.n(2),model1.nb(2)),opExtension(model1.n(1),[nbz1 model1.nb(1)]))*din(:);
    end
    model1.o=[model1.o(1)-nbz1*model1.d(1) model1.o(2)-model1.nb(2)*model1.d(2) 0];
    model1.n=model1.n+[sum([nbz1 model1.nb(1)]) 2*model1.nb(2) 0];
end


%% Mesh acquisition grid

% Check Acquisiyion Geometry
Check_Acquisition_Geom(model1)
if ~strcmp(model1.type,'marine')
    [model1.xrec,model1.yrec,model1.zrec] = ndgrid(model1.xrec,model1.yrec,model1.zrec);
end
%% Domain decomposition
[v,model,dens,ani] = setup_for_domain_decomp(v,model1,dens,ani);
if mode==1
	[data_dis,~,~,~]= setup_for_domain_decomp(din,model1,[],[]);
end
% Number of workers
Num_of_worker=parpool_size();
if Num_of_worker==0
	Num_of_worker=1;
end% number of labs


% Domain decomposition size
if isfield(model1,'ddcompx')
	Num_of_model_split = model1.ddcompx * model1.ddcompz * model1.ddcompy;
else
	model1.ddcompx=1;
	 model1.ddcompz =1;
	 model1.ddcompy=1;
	Num_of_model_split=1;
end


if  Num_of_model_split >  Num_of_worker
    error(['Number of works (',num2str( Num_of_worker),') can not be smaller ' ...
        'than number of model splits (',num2str( Num_of_worker),'); Need to be N ' ...
        'times larger, then you do not waste works.']);
end

%% Source distribution
% Check if simultenaous shots
if size(model1.xsrc,1)>1
    model1.multi=size(model1.xsrc,1);
else
    model1.multi=0;
end

Num_of_shot_band = floor( Num_of_worker/ Num_of_model_split);
Num_of_all_shots = size(model1.xsrc,2);
shots_parition	  = round(linspace(0,Num_of_all_shots,Num_of_shot_band+1));
src_dis = Composite();
nt=length(0:model1.dt:model1.T);
ntn=length(model1.NyqT);
if mode~=1
    din=reshape(din,length(model1.NyqT),size(model1.xrec,1),Num_of_all_shots);
	data_dis=Composite();
end
for labi = 1:Num_of_worker
    shot_band_index			= ceil(labi/Num_of_model_split);
    shot_id		= shots_parition(shot_band_index)+1:shots_parition(shot_band_index+1);
    slx_dis		= model1.xsrc(:,shot_id);
    sly_dis     = model1.ysrc(:,shot_id);
    slz_dis     = model1.zsrc(:,shot_id);
    src_dis{labi}=struct('slx_band',slx_dis,'sly_band',sly_dis,'slz_band',slz_dis,'shots_id',shot_id,'max_number_of_shot_all_worker', max(diff(shots_parition)));
    if mode~=1
        data_dis{labi}=din(:,:,shot_id);
    end
end
%% Slove the Born modelling or RTM

spmd
	model.dt=model1.dt;
	model.f0=model1.f0;
    % Number of local shots
    Num_of_shots_loc = size(src_dis.slx_band,2);
    if mode==1
        data=cell(Num_of_shots_loc,1);
    else
        data=zeros(length(model.mmz_loc),length(model.mmx_loc),length(model.mmy_loc));
		if nout>1
        	data2=cell(Num_of_shots_loc,1);
		end
    end
    ntrace_loc=0;
    
    for i=1:Num_of_shots_loc
        % Local source geometry
        model.slx	= src_dis.slx_band(:,i);
        model.sly	= src_dis.sly_band(:,i);
        model.slz	= src_dis.slz_band(:,i);
        
		if labindex==1
			if mode==1
        		disp(['J for shot  ' num2str(i) ' over ' num2str(length(src_dis.shots_id)) ' at position ' num2str(model.slx)]);
			else
				disp(['J^T for shot  ' num2str(i) ' over ' num2str(length(src_dis.shots_id)) ' at position ' num2str(model.slx)]);
			end
		end
        model.numsave=src_dis.shots_id(i);
        % Forward Wavefield checkpoints for adjoint (RTM)
        if nt/ntn>1
            model.nsub=sort(unique([1:floor(nt/ntn):nt nt-1 nt]));
        else
            model.nsub=1:nt;
        end

        % Receiver geometry
        if strcmp(model1.type,'marine')
            model.rlx	= unique(model1.xrec(:,src_dis.shots_id(i)));
            model.rly	= unique(model1.yrec(:,src_dis.shots_id(i)));
            model.rlz	= unique(model1.zrec(:,src_dis.shots_id(i)));
            [model.rlx,model.rly,model.rlz]=ndgrid(model.rlx,model.rly,model.rlz);
            model.rlx= model.rlx(:);
            model.rlz= model.rlz(:);
            model.rly= model.rly(:);
        else
            model.rlx	= model1.xrec(:);
            model.rly	= model1.yrec(:);
            model.rlz	= model1.zrec(:);
        end
        
        model = locate_rcv_position(model);
        model = locate_shot_position(model);
        
        % In case of simultenaous sources
        if model1.multi>1
            model.multi=model1.multi;
            srcloc=zeros(size(model.xsrc,1),nt);
            for is=1:length(model.xsrc(:))
                srcloc(is,:)=interpft(model1.R{is}*source,nt);
            end
        else
            srcloc=source';
        end
		
		
		if mode==1
			model.NyqT=model1.NyqT;
			J=opJ(v,model,srcloc,0,dens,ani);
			data{i}=J*(-data_dis(:));
			ntrace = length(model.rec_idx);
			J=[];
		else
			model.nconv=model.nsub;
	        model.save=model1.save;
	        model.NyqT=model1.NyqT;
			% Forward wavefield, we need it for the cross correlations
			F=opF(v,model,dens,ani);
			D=F*srcloc(:);
			F=[];
			if nout>2
				if strcmp(model.save,'RAM')
					data2{i}=D{1}(:)-vec(data_dis(:,model.rec_idx,i));
					J=opJ(v,model,srcloc,D{2},dens,ani);
				else
					data2{i}=D(:)-vec(data_dis(:,model.rec_idx,i));
					J=opJ(v,model,srcloc,[],dens,ani);
				end
			else
				if strcmp(model.save,'RAM')
					J=opJ(v,model,srcloc,D{2},dens,ani);
				else
					J=opJ(v,model,srcloc,[],dens,ani);
				end
			end
			% Check local size of the data to fit the operator size
			if length(model.rec_idx)==0
				dataloc=J'*zeros(size(J,1),1);
			else
				if nout>2
					dataloc=J'* data2{i};
				else
					dataloc=J'* vec(data_dis(:,model.rec_idx,i));
				end
			end
			ntrace = length(model.rec_idx);
			J=[];
		end

		if mode~=1
			data=data+dataloc;
		end
		
		ntrace_loc = ntrace_loc + ntrace;

    end
    
    % Build output size according to mode
    if mode==1
        if Num_of_shots_loc>0
            recidx=model.rec_idx;
        end
        
        codistr	   = codistributor1d(2,codistributor1d.unsetPartition,[1,numlabs]);
        ntrace_loc		= codistributed.build(ntrace_loc,codistr,'noCommunication');
    else
        data=reshape(data,length(model.mmz_loc),length(model.mmx_loc),length(model.mmy_loc));
        data=data(model.u_idz,model.u_idx,model.u_idy);
        
        if nout>1
            if Num_of_shots_loc>0
                recidx=model.rec_idx;
            end
            
            codistr	   = codistributor1d(2,codistributor1d.unsetPartition,[1,numlabs]);
            ntrace_loc		= codistributed.build(ntrace_loc,codistr,'noCommunication');
        end
    end
	labBarrier;
end
%----------------------------------------- Gather result ----------------------------------------%

if mode==1
    if Num_of_model_split>1
        du=zeros(length(model1.NyqT),length(model1.xrec(:)),size(model1.xsrc,2));
        if strcmp(model1.type,'full')
            for i=1:Num_of_worker
                sloc=src_dis{i};
                dloci=data{i};
                Num_of_shots_loc = size(sloc.slx_band,2);
                for s=1:Num_of_shots_loc
                    dlocj=dloci{s};
                    du(:,recidx{i},sloc.shots_id(s))=reshape(dlocj,length(model1.NyqT),length(recidx{i}));
                end
            end
            du=du(:);
        else
            disp('Data is not reshaped correctly due to domain decomposition for this type of acquisition yet...');
        end
    else
        du = gather(output_shot_record(data,model1,ntrace_loc,1));
        du=gather(du);
    end
else
    du=gather_gradient(data,model1);
    du=du(:);
    
    if nargout>1
        if Num_of_model_split>1
            r=zeros(length(model1.NyqT),length(model1.xrec(:)),size(model1.xsrc,2));
            if strcmp(model1.type,'full')
                idx=[];
                for i=1:Num_of_worker
                    sloc=src_dis{i};
                    dloci=data2{i};
                    Num_of_shots_loc = size(sloc.slx_band,2);
                    for s=1:Num_of_shots_loc
                        dlocj=dloci{s};
                        r(:,recidx{i},sloc.shots_id(s))=reshape(dlocj,length(model1.NyqT),length(recidx{i}));
                    end
                end
                r=r(:);
            else
                disp('Data is not reshaped correctly due to domain decomposition for this type of acquisition yet...');
            end
        else
            r = gather(output_shot_record(data2,model1,ntrace_loc,1));
            r=gather(r);
        end
    end
end
end