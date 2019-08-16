function[fun,grad]=GS(v,model1,source,dataM,dens,ani)
% GS -- Gradient and functional computation.
%
% Function that generates the gradient and data misfit
%
% INPUTS
%
%   'model1' : stucture containing the model parameters
%       'v'  :  square slowness [s^2/m^2]
%     'dens' : density [g/cm^3]
%     'source' : source term  (1 x nt)
%     'anis' : anisotropy parameters (epsilon,delta,theta in a
%     structure)
%     'dataM' : True data (does not support files currently)
%     'source' : source term  (1 x nt)
% OUTPUT from multiplication
%
%   'grad':  gradient
%   'fun':  misfit
%
%
% Author:Mathias Louboutin
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
% Date: July, 2015


%% Extend the model for PML's
if model1.freesurface
    nbz1=model1.space_order/2;
else
    nbz1=40;
end
% Redifine dimensions
model1.nb=[40 40 40];
if nargin>=5 && ~isempty(dens)
    if model1.n(3)~=1
        dens=kron(opExtension(model1.n(3),model1.nb(3)),opExtension(model1.n(2),model1.nb(2)),opExtension(model1.n(1),[nbz1 model1.nb(1)]))*dens(:);
    else
        dens=kron(opExtension(model1.n(2),model1.nb(2)),opExtension(model1.n(1),[nbz1 model1.nb(1)]))*dens(:);
    end
else
    dens=[];
end
if nargin==6 && ~isempty(ani)
    ani.epsilon=kron(opExtension(model1.n(2),model1.nb(2)),opExtension(model1.n(1),[nbz1 model1.nb(1)]))*ani.epsilon(:);
    ani.delta=kron(opExtension(model1.n(2),model1.nb(2)),opExtension(model1.n(1),[nbz1 model1.nb(1)]))*ani.delta(:);
    ani.theta=kron(opExtension(model1.n(2),model1.nb(2)),opExtension(model1.n(1),[nbz1 model1.nb(1)]))*ani.theta(:);
else
    ani=[];
end

if model1.n(3)~=1
    v=kron(opExtension(model1.n(3),model1.nb(3)),opExtension(model1.n(2),model1.nb(2)),opExtension(model1.n(1),[nbz1 model1.nb(1)]))*v(:);
    model1.o=[model1.o(1)-nbz1*model1.d(1) model1.o(2)-model1.nb(2)*model1.d(2) model1.o(3)-model1.nb(3)*model1.d(3)];
    model1.n=model1.n+[sum([nbz1 model1.nb(1)]) 2*model1.nb(2) 2*model1.nb(3)];
else
    v=kron(opExtension(model1.n(2),model1.nb(2)),opExtension(model1.n(1),[nbz1 model1.nb(1)]))*v(:);
    model1.o=[model1.o(1)-nbz1*model1.d(1) model1.o(2)-model1.nb(2)*model1.d(2) 0];
    model1.n=model1.n+[sum([nbz1 model1.nb(1)]) 2*model1.nb(2) 0];
end
%% Time setup
ntn=length(model1.NyqT); 
nt=length(0:model1.dt:model1.T);
%% Mesh acquisition grid

% Check Acquisiyion Geometry
Check_Acquisition_Geom(model1)
if ~strcmp(model1.type,'marine')
    [model1.xrec,model1.yrec,model1.zrec] = ndgrid(model1.xrec,model1.yrec,model1.zrec);
end
%% Domain decomposition
[v,model,dens,ani] = setup_for_domain_decomp(v,model1,dens,ani);

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
if size(model1.xsrc,2)>1
    model1.multi=size(model1.xsrc,2);
else
    model1.multi=0;
end

Num_of_shot_band = floor( Num_of_worker/ Num_of_model_split);
Num_of_all_shots = size(model1.xsrc,2);
shots_parition	  = round(linspace(0,Num_of_all_shots,Num_of_shot_band+1));
dataM=reshape(dataM,ntn,length(model1.xrec(:)),Num_of_all_shots);
src_dis = Composite();
data_dis=Composite();

for labi = 1:Num_of_worker
    shot_band_index			= ceil(labi/Num_of_model_split);
    shot_id		= shots_parition(shot_band_index)+1:shots_parition(shot_band_index+1);
    slx_dis		= model1.xsrc(:,shot_id);
    sly_dis     = model1.ysrc(:,shot_id);
    slz_dis     = model1.zsrc(:,shot_id);
    src_dis{labi}=struct('slx_band',slx_dis,'sly_band',sly_dis,'slz_band',slz_dis,'shots_id',shot_id,'max_number_of_shot_all_worker', max(diff(shots_parition)));
    data_dis{labi}=dataM(:,:,shot_id);
end
%% FWI update
fh = @Jls;

spmd
    % Number of local shots
    Num_of_shots_loc = size(src_dis.slx_band,2);
    
    % Init local gradient and misfit
    g   = zeros(length(model.mmz_loc),length(model.mmx_loc),length(model.mmy_loc));
	f0=0;

    
    for i=1:Num_of_shots_loc
        % Local source geometry
        model.slx	= src_dis.slx_band(:,i);
        model.sly	= src_dis.sly_band(:,i);
        model.slz	= src_dis.slz_band(:,i);
        
        % if labindex==1
            % disp(['shot = ',num2str(src_dis.slx_band(i))]);
        % end
		
        if nt/ntn>1
            model.nsub=sort(unique([1:floor(nt/ntn):nt nt-1 nt]));
        else
            model.nsub=1:nt;
        end
		model.nconv=model.nsub;
		
		model.numsave=src_dis.shots_id(i);
		model.save=model1.save;
		model.NyqT=model1.NyqT;


        if strcmp(model1.type,'marine')
            model.rlx	= unique(model1.xrec(:,src_dis.shots_id(i)));
            model.rly	= unique(model1.yrec(:,src_dis.shots_id(i)));
            model.rlz	= unique( model1.zrec(:,src_dis.shots_id(i)));
            [model.rlx,model.rlz,model.rly]=ndgrid(model.rlx,model.rlz,model.rly);
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

		if strcmp(model.save,'disk')
			model.numsave=src_dis.shots_id(i);
		end
		
		% Data for smooth model
		F=opF(v,model,dens,ani);
		Dt= F*source;
		D=data_dis(:,model.rec_idx,i);
		%% Filter input data
		if isfield(model,'filter')
			[D,fk,kx]=fktran(D,1e-3*model.NyqT,unique(model.rlx),1);
			[ff,kkx] = ndgrid(fk,kx);
			F2 = zeros(length(fk),length(kx));
			F2(ff <= 1e3*model1.f0) = 1;
			Sm= opSmooth(size(F2,1),20);
			F2=Sm*F2;
			D=fktran(F2.*D(:,1:size(F2,2)),1e-3*model.NyqT,unique(model.rlx),-1); 
			D=[D zeros(nd(1),nd(2)-size(D,2))];    
		end
		D=D/norm(D(:));
		
		if strcmp(model.save,'disk')
			dloc=reshape(Dt(:),size(D));
		else
			dloc=reshape(Dt{1}(:),size(D));
		end

		%% Filtr data to match synthetic
		if ~isempty(dloc)
			dloc=dloc(:)/norm(dloc(:));
		end
		if strcmp(model.save,'disk')
			floc = .5*norm(Dt(:)-D(:)).^2;
			Jt=opJ(v,model,source,[],dens,ani);
			if isempty(Dt)
				D=zeros(size(Jt,1),1);
			else
				D=(dloc-D(:));
			end
		else
			floc = .5*norm(dloc(:)-D(:)).^2;
			Jt=opJ(v,model,source,Dt{2},dens,ani);
			D=(dloc-D(:));
		end
		f0=f0+floc;
		g=g + Jt'*D(:);
		Jt=[];
		dloc=[];
		D=[];
		Dt=[];
		F=[];
    end

    g=reshape(g,length(model.mmz_loc),length(model.mmx_loc),length(model.mmy_loc));
    g=g(model.u_idz,model.u_idx,model.u_idy);
	
	xloc=[];
	model=[];
	model=[];
	g0=[];
	dm=[];
	labBarrier;
end

grad = gather_gradient(g,model1);
fun=0;
for i=1:Num_of_worker

		fun    = fun  + f0{i};
end

grad=grad(:);
end
