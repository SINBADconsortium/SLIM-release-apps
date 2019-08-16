function s=Bak_propag(v,model1,source,dens,ani)
	
% Bak_propag -- Back propagation function, only used for the dottest
%
% INPUTS to the constructor
%
%   'model1' : stucture containing the model parameters
%       'v'  :  square slowness [s^2/m^2]
%     'dens' : density [g/cm^3]
%     'source' : source term  (1 x nt)
%     'anis' : anisotropy parameters (epsilon,delta,theta in a
%     structure)
%
% OUTPUT from multiplication
%
%   's':    source
%
%
% Author : Mathias Louboutin
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
% Date : July, 2015

%% Extend the model for PML's

if model1.freesurface
    nbz1=model1.space_order/2;
else
    nbz1=40;
end
% Redifine dimensions
model1.nb=[40 40 40];
if nargin>=4 && ~isempty(dens)
	if model1.n(3)~=1
	    dens=kron(opExtension(model1.n(3),model1.nb(3)),opExtension(model1.n(2),model1.nb(2)),opExtension(model1.n(1),[nbz1 model1.nb(1)]))*dens(:);
	else
	    dens=kron(opExtension(model1.n(2),model1.nb(2)),opExtension(model1.n(1),[nbz1 model1.nb(1)]))*dens(:);
    end
else
    dens=[];
end
if nargin==5 && ~isempty(ani)
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
end

Num_of_shot_band = floor( Num_of_worker/ Num_of_model_split);
Num_of_all_shots = size(model1.xsrc,2);
shots_parition	  = round(linspace(0,Num_of_all_shots,Num_of_shot_band+1));
source=reshape(source,length(model1.NyqT),size(model1.xrec,1),Num_of_all_shots);
nt=length(0:model1.dt:model1.T);

src_dis = Composite();
data_dis=Composite();

for labi = 1:Num_of_worker
    shot_band_index			= ceil(labi/Num_of_model_split);
    shot_id		= shots_parition(shot_band_index)+1:shots_parition(shot_band_index+1);
    slx_dis		= model1.xsrc(:,shot_id);
	sly_dis     = model1.ysrc(:,shot_id);
	slz_dis     = model1.zsrc(:,shot_id);
    src_dis{labi}=struct('slx_band',slx_dis,'sly_band',sly_dis,'slz_band',slz_dis,'shots_id',shot_id,'max_number_of_shot_all_worker', max(diff(shots_parition)));
    data_dis{labi}=source(:,:,shot_id);
end

%% Slove the Adjoint model

spmd

    % Number of local sources
	Num_of_shots_loc = size(src_dis.slx_band,2);
	ntrace_loc=0;
	data=zeros(Num_of_shots_loc,nt);


	for i=1:Num_of_shots_loc

        model.slx	= src_dis.slx_band(:,i);
	    model.sly	= src_dis.sly_band(:,i);
	    model.slz	= src_dis.slz_band(:,i);
	
		if labindex==1
	   		disp(['shot = ',num2str(src_dis.slx_band(i))]);
		end
    
		model.nsub=[];
		model.save='RAM';
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

        F=opF(v,model,dens,ani);
		
	    data(i,:)=F'*vec(data_dis(:,:,i));

	    ntrace_loc = ntrace_loc + length(model.rec_idx);

	end
end

s=zeros(Num_of_all_shots,nt);
for i=1:Num_of_worker
	dloc=data{i};
	sloc=src_dis{i};
	Num_of_shots_loc = size(sloc.slx_band,2);
	for j=1:Num_of_shots_loc
		jloc=sloc.shots_id(j);
		s(jloc,:)=dloc(j,:);
	end
end