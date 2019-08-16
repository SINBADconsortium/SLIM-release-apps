function [f,g] = Jls(m,Q,D,model,dens,ani)
% least-squares misfit for time domain FWI
%
% use:
%   [f,g] = Jls(m,Q,D,model,deani)
%
% input:
%   m - gridded squared slowness [km^2/s^2], see also F.m
%   Q - gridded source functions, see also F.m
%   D - `observed' data consistent with output of F(m,Q,model)
%   model - Local struct with model parameters for opFM3d.m and opJ3D.m
%   ani (optional) - Structure with Thomson parameters and TTI tilt angle
%
% output
%   f - misfit
%   g - gradient
%

% Author: Mathias louboutin
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: July, 2015


F=opF(m,model,dens,ani);
Dt= F*Q;

%% Filter input data
if isfield(model,'filter')
	[D,fk,kx]=fktran(D,1e-3*model.NyqT,unique(model.rlx),1);
	[ff,kkx] = ndgrid(fk,kx);
	F2 = zeros(length(fk),length(kx));
	F2(ff <= 1e3*model.f0) = 1;
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
    f = .5*norm(Dt(:)-D(:)).^2;
    Jt=opJ(m,model,Q,[],dens,ani);
    if isempty(Dt)
        D=zeros(size(Jt,1),1);
    else
        D=(dloc-D(:));
    end
else
    f = .5*norm(dloc(:)-D(:)).^2;
    Jt=opJ(m,model,Q,Dt{2},dens,ani);
    D=(dloc-D(:));
end

if nargout > 1
    g = Jt'*D(:);
end
end