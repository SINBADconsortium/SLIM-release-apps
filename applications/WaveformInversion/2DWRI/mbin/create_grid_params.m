function [model_PDE,int_model2PDE, int_PDE2model,Q_scale_factor,int_model2PDE_z]=create_grid_params(m,model,freqs,params,units)
%
% This function takes a model + corresponding grid and resamples it
% according to the highest frequency that needs to be modeled. Output
% includes model information structure and interpolation/resampling
% operators
%
% input:
% m       - model to be resampled to a new grid
% model   - structure with modeling and geometry info
% freqs   - frequencies in the current frequency batch
% params  - parameters for the points per wavelength and PML size
% units   - indicates if the input model is in 'm/s' or in '(s/m)^2'.
%
% output:
% model_PDE      - structure with all modeling settings and information for
% the PDE solves 
% int_model2PDE  - SPOT operator which resamples from model grid to PDE
% grid
% int_PDE2model  - SPOT operator which resamples from PDE grid to model
% grid
% Q_scale_factor - 
%
% Author: Bas Peters
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%         
% Date: June, 2015

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.



model_PDE = model;

%get parameters
nppw    = params.nppw; %get desired minimum number of points per wavelength
nwl_pml = params.nwl_pml; %get the desired number of wavelengths in the PML

%convert m from slowness squared to m/s
if strcmp(units,'(s/m)^2')==1
    v=1e3./sqrt(m);
elseif strcmp(units,'m/s')==1
    v = m;
end

min_vel     = min(v(:));
max_vel     = max(v(:));
min_lambda  = min_vel/max(freqs);
max_lambda  = max_vel/min(freqs);

pml_width = nwl_pml*max_lambda;         %physical size of the PML
model_PDE.d = [min_lambda/nppw min_lambda/nppw];  %spacing between grid points
%now we have the desired d, based on the nppw, but we also want this d to
%generate a new grid which has the same endpoint as the original grid, in
%order to use extrapolation when going from a PDE grid to a model grid
% %which is finer.
c1 = model.z(end)/model_PDE.d(1);
c2 = model.x(end)/model_PDE.d(2);
c1 = ceil(c1);
c2 = ceil(c2);
model_PDE.d = [model.z(end) model.x(end)]./[c1 c2];

model_PDE.nb  = ceil([pml_width pml_width]./model_PDE.d);   %number of grid points in the PML

%new axis vectors
x = [model.o(2):model_PDE.d(2):(model.n(2)-1)*model.d(2)];
z = [model.o(1):model_PDE.d(1):(model.n(1)-1)*model.d(1)];

int_model2PDE = opKron(opLInterp1D_mod(model.x,x),opLInterp1D_mod(model.z,z));
int_model2PDE_z =opLInterp1D_mod(model.z,z);
int_PDE2model = opKron(opLInterp1D_mod(x,model.x),opLInterp1D_mod(z,model.z));
%int_model2PDE = opKron(opLInterp1D(model.x,x),opLInterp1D(model.z,z));
%int_PDE2model = opKron(opLInterp1D(x,model.x),opLInterp1D(z,model.z));

model_PDE.x = x;
model_PDE.z = z;
[~, ~, model_PDE.n] = grid2odn(z,x);

model_PDE.freq  = freqs;
model_PDE.nfreq = length(freqs);

Q_scale_factor = 100*(1/(model_PDE.d(1)*model_PDE.d(2)));