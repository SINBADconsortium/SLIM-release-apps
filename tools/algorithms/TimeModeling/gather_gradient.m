function grad = gather_gradient(g,model1)

%% Compute the gradient on the full grid
% Uses the gradient computed for each source and domain and put everything
% back together at the correct place 
% Author : Xiangli edited by Mathias Louboutin
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
% Date : July, 2015
nx		=	model1.n(2);
ny		=	model1.n(3);
nz		=	model1.n(1);
ngx		=   round(linspace(0,nx,model1.ddcompx+1));
ngz		=	round(linspace(0,nz,model1.ddcompz+1));
ngy		=	round(linspace(0,ny,model1.ddcompy+1));
    
Num_of_worker=length(g);    
Num_of_model_split=model1.ddcompx*model1.ddcompy*model1.ddcompz;
grad=zeros(prod(model1.n),1);
for i=1:Num_of_worker
        model_parition_index	= mod(i-1,Num_of_model_split)+1;
		[ddcompz_id,ddcompx_id,ddcompy_id]	= ind2sub([model1.ddcompz,model1.ddcompx,model1.ddcompy],model_parition_index);


		% find model partition
		x_id = ngx(ddcompx_id)+1:ngx(ddcompx_id+1);
		z_id = ngz(ddcompz_id)+1:ngz(ddcompz_id+1);
		y_id = ngy(ddcompy_id)+1:ngy(ddcompy_id+1);
	
		[ZZ,XX,YY] = ndgrid(z_id,x_id,y_id);
		
		Model_idx  = sub2ind([nz,nx,ny],ZZ(:),XX(:),YY(:));
		
	
		grad(Model_idx)    = grad(Model_idx)  + vec(g{i});
end

grad=reshape(grad,model1.n);
if ny==1
	if model1.freesurface
		grad=grad(model1.space_order/2+1:(end-model1.nb(2)),model1.nb(2)+1:(end-model1.nb(2)));
	else
		grad=grad(1+model1.nb(1):(end-model1.nb(1)),model1.nb(2)+1:(end-model1.nb(2)));
	end
else
	if model1.freesurface
		grad=grad(model1.space_order/2+1:(end-model1.nb(2)),model1.nb(2)+1:(end-model1.nb(2)),model1.nb(3)+1:(end-model1.nb(3)));
	else
		grad=grad(1+model1.nb(1):(end-model1.nb(2)),model1.nb(2)+1:(end-model1.nb(2)),model1.nb(3)+1:(end-model1.nb(3)));
	end
end

if isfield(model1,'water')
	if length(model1.water)==1
		grad(1:model1.water,:)=0;
	else
		grad(model1.water)=0;
	end
end