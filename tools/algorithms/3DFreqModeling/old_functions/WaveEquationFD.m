%-------------------------------------------------------------------------------
% Computes the value of the predicted data, gradient and value of the misfit
% function using a finite difference frequency domain modeling of the wave 
% equation.
% 
% This is the core function of the 3D frequency domain modeling code.
% 
% Use: 
%   [misfit, grad, d_pred] = WaveEquationFD(model, acq, opts, flog)
%
% Input:
%  model      - Contains information about the current model to which the 
%               forward modeling is to be applied. Read further for more details
%               on its content.
%  acq        - Contains information about the acquisition and its geometry. 
%               Read further for more details on its content.
%  
%  opts       - [OPTIONAL] argument containing some tweaking parameters. Read 
%               further for more details.
%  flog       - [OPTIONAL] file identifier for the (already open) log file. 
%               Set to 0 (or simply do not provide it) if no verbosity is 
%               desired.
%
% Output:
%   misfit    - value of the misfit function for the given observed data
%   grad      - gradient of the misfit function evaluated at given model for the
%               given data
%   d_pred    - data cube as an array whose dimension is (xrec,yrec,zrec,xsrc,
%               ysrc,zsrc,nfreq). This is not vectorized.
%
%%%%%%%%%%%%%%%%%%%
%
% IMPORTANT: Be careful with the output arguments. If you do not need the 
% predicted data, for instance, use
%      [misfit, grad] = WaveEquationFD(model, acq, opts, flog)
% rather than 
%      [misfit, grad, ~] = WaveEquationFD(model, acq, opts, flog)
% because in the later, Matlab thinks that there are 3 output arguments. 
% To circumvent this, you can use the opts structure to force this function to
% only output specific variables, regardless of how you call this function.
%
%%%%%%%%%%%%%%%%%%%%
%
% Structure of "model" (input):
% ---------------------------------------
% *Must* contain:
%  model.v    - velocity model in meters per second or slowness squared; this
%               is a *3D object* with dimension (nx,ny,nz), NOT A VECTOR!
%  model.ov   - origin of the velocity model
%  model.dv   - Array with the distance between each point in each direction
%               for this velocity model. This is used to compute the physical
%               size of the model - which should not change under any 
%               circunstance.
%  model.unit - String set either to "m/s" (meters per second) or "s2/m2" 
%               (slowness squared)
% 
% 
% Structure of "acq" (input):
% ---------------------------------------
% *Must* contain:
%  acq.freq   - array with frequencies
%  acq.{zsrc,xsrc,ysrc} - vectors describing source array
%  acq.{zrec,xrec,yrec} - vectors describing receiver array.
%  acq.sources- array with dimensions (xsrc*ysrc*xsrc, nsrc) containing the 
%               source position, where nsrc is the number of sources.
% Optional:
%  acq.f0     - peak frequency of Ricker wavelet, 0 for no wavelet.
%  acq.t0     - phase shift [s] of wavelet. Default to 0.
%  data       - observed data, used to compute the value of the value of the
%               misfit function and the gradient. Can be neglected if the only
%               output desired it d_pred
%
%
% Structure of "opts" (input, optional):
% ---------------------------------------
% 
% par_solver  - Can be passed as a structure or a (cell) array of structures.
%               If this is a single structure, this will be used as parameter
%               for the solver for all frequencies. If this is an array, it 
%               must have nfreq entries, where nfreq is the number of frequency
%               slices in your data. In this case, par_solver{k} will be used
%               for the k-th frequency. Notice that some of the elements of the
%               cell can be empty e.g. par_solver{m}=[]. In this case, the 
%               parameters for the m-th frequency are set to default, but not 
%               for the other frequency slices. This can be totally ignored, in
%               which case the default parameters are used for all frequency 
%               slices.
% par_helm    - Same as par_solver, but this is passed to helmholtz_solver
%               function.
% rng_freq    - An array containing the index of the frequencies of interest.
%               For instance, if this value is set to [1 2 5], this function 
%               will take into account only acq.freq(1), acq.freq(2) and 
%               acq.freq(5), and ignore all other possible frequencies.
% rng_src     - Same as rng_freq, but for acq.sources.
%
% weight_fun  - Function handler to the weighting to be used in the loop over
%               sources for generating the computed data and misfit function.
%               The interface should look like 
%               weight_fun = @(l)(...)
%               where "l" is source index. The default function is called
%               "offset_mask". 
% src_est     - Function handler to the source estimage function, to be used 
%               in the loop over sources. Only relevant if data is being 
%               output. The interface should look like
%               @(u,l,k)(...);
%               where "u" is the unscaled computed data, "l" is the source
%               index and "k" is the frequency index. The default function
%               is called "src_estimate".
% 
% IMPORTANT: If outputing d_pred, the outputed data will have the same dimension
%            of acq.nd, regardless of how many sources/frequencies were 
%            actually used. The data for the unused sources/frequencies will 
%            be set to zero. You have been warned!! 
%       
% output_grad   - Force this function to ignore the computation of the gradient
% output_d_pred - Forces this function to throw away the predicted data after
%                 using it (to save memory).
% output_misfit - Force this function to ignore the computation of the value
%                 the misfit function. 
%
% Author: Rafael Lago
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: September, 2014
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%-------------------------------------------------------------------------------
function [misfit, grad, d_pred] = WaveEquationFD(model, acq, opts, flog)

if ~exist('opts','var')
  opts = [];
end

if ~exist('flog','var')
  flog = 0;
end

% Acquisition dimensions
%-------------------------------
nrec    = length(acq.xrec)*length(acq.yrec)*length(acq.zrec);
rec_dom = acq.nd(1:3);

rng_freq = check_field(opts,'rng_freq',1:length(acq.freq));
  nfreq = length(acq.freq);
rng_src = check_field(opts,'rng_src', 1:size(acq.sources,2));
    nsrc = size(acq.sources,2);

% Parse parameters
%-------------------------------
% Check if parameters for the linear solver were given
par_solver = check_field(opts,'par_solver',cell(nfreq,1));
if length(par_solver)==1
  tmp(1:nfreq) = {par_solver};
  par_solver = tmp;
end
% Check if parameters for the discretization was passed
par_helm = check_field(opts,'par_helm',cell(nfreq,1));
if length(par_helm)==1
  tmp(1:nfreq) = {par_helm};
  par_helm = tmp;
end

% Define wavelet
%-------------------------------
t0 = check_field(acq,'t0',0);
f0 = check_field(acq,'f0',0);

wlet = exp(1i*2*pi*acq.freq*t0);
if f0>0
  % Ricker wavelet with peak-frequency acq.f0
  wlet = (acq.freq).^2.*exp(-(acq.freq/acq.f0).^2).*wlet;
end

% Check what needs to/can be output by this function
output_misfit = check_field(opts,'output_misfit',nargout>0);
output_grad   = check_field(opts,'output_grad'  ,nargout>1);
output_d_pred = check_field(opts,'output_d_pred',nargout>2);

if output_misfit || output_grad
  % This makes easier to compute the data residual
  % Notice that this change remains inside this function. Outside this function
  % the acq.data is still a 7D object.
  acq.data = reshape(acq.data,nrec,nsrc,nfreq);
  if ~isfield(acq,'data')
      error(['Can only compute misfit or the gradient if observed data is'...
            ' passed as acq.data\n Aborting...']);
  end
  % Check if a weighting function was passed; if not, set it to identity
  weight_fun = check_field(opts,'weight_fun',@(l)(offset_mask(acq,0,l)));
  src_est    = check_field(opts,'src_est',   @(u,l,k)(src_estimate(weight_fun(l),u,acq.data(:,l,k))));
end


% Loop over frequencies
% and sources
% Well... "main loop"
%------------------------------
misfit = 0;
d_pred = [];
grad   = zeros(model.nv);
if output_d_pred
  d_pred = zeros(nrec,nsrc,nfreq);
end
freq_counter = 1;
plog(flog,'\n* Start modeling...:\n');
for k = rng_freq

  plog(flog,'* Frequency ', freq_counter ,' of ', length(rng_freq) ,': ', acq.freq(k),'Hz \n');
  plog(flog,'* Setting up the PDE solver...\n');
  
  %-------------
  if output_grad
      [H, A] = helmholtz_solver(model,acq.freq(k),par_helm{k},par_solver{k},flog);
      
      plog(flog,'* Obtaining the derivative of the discrete operator...\n');
      [dH, ~ ] = helmholtz_3d_derivative(H.pg2cg(model.v),H.d,acq.freq(k),model.unit);
      dH = conj(Htransp(dH,H.idx));
  else
      H = helmholtz_solver(model,acq.freq(k),par_helm{k},par_solver{k},flog);
  end
  plog(flog,'* Chosen solver: ', H.solver_name ,'\n');
  Nt = prod(H.nt);
  
  % Abbreviations:
  % cg - computational grid
  % rg - receivers grid
  % sg - sources grid
  % pg - physical grid
  %----------------------------------------------------------------------------
  [xcg,ycg,zcg] = odn2grid(model.ov,H.d,H.n);   %No PML here
  [xpg,ypg,zpg] = odn2grid(model.ov,H.dv,H.nv); %No PML here
  
  Lpad   = opKron(opExtension(H.n(3),[H.pml.top H.pml.bottom],0),...
                  opExtension(H.n(2),H.pml.y,0),...
                  opExtension(H.n(1),H.pml.x,0));
  Lcg2rg = opKron(opLInterp1D(zcg,acq.zrec),opLInterp1D(ycg,acq.yrec),opLInterp1D(xcg,acq.xrec));
  Lcg2sg = opKron(opLInterp1D(zcg,acq.zsrc),opLInterp1D(ycg,acq.ysrc),opLInterp1D(xcg,acq.xsrc));
  Lcg2pg = opKron(opLInterp1D(zcg,zpg),opLInterp1D(ycg,ypg),opLInterp1D(xcg,xpg));
  
  sg2cg = @(x)(reshape(Lpad*(Lcg2sg'*x(:)),H.nt));
  rg2cg = @(x)(reshape(Lpad*(Lcg2rg'*x(:)),H.nt));
  cg2rg = @(x)(Lcg2rg*(Lpad'*x(:)));
  
  % Volume of the cubic cell; needed to scale the source more or less properly.
  cub_vol = prod(H.d)/prod(H.dv);
  
  plog(flog,'* System size: ', Nt,' -  Number of shots: ',size(acq.sources,2),'\n');
  plog(flog,'* ---------------------------------------------------\n');
  for l = rng_src
      
      %-------------------------------------------------
      % Forward Problem
      %       uᵢ = H⁻¹qᵢ
      %-------------------------------------------------
      Uk = H.solve(cub_vol*wlet(k)*sg2cg(acq.sources(:,l)),zeros(Nt,1));
        
      %---------------
      if output_d_pred
        d_pred(:,l,k) = cg2rg(Uk);
      end

      %-------------------------------------------------
      % Define data residual
      %   d_res =   dᵢ + H⁻¹qᵢ
      %-------------------------------------------------
      if output_misfit || output_grad
        wsrc  = src_est(cg2rg(Uk),l,k);
        Uk = Uk*wsrc;
        d_res = weight_fun(l)*(cg2rg(Uk)-acq.data(:,l,k));
      end
      
      %-------------------------------------------------
      % Misfit evaluation
      %      f =   f    + || dᵢ + H⁻¹qᵢ||
      %-------------------------------------------------
      if output_misfit
        misfit = misfit + .5*norm(d_res)^2;
      end
      
      %-------------------------------------------------
      % Adjoint/Backward Problem
      %       vᵢ = H⁻*(dᵢ - H⁻¹qᵢ)
      %-------------------------------------------------
      if output_grad
        Vk = A.solve(rg2cg(cub_vol*(weight_fun(l)*d_res)),zeros(Nt,1));
        
      %-------------------------------------------------
      % ∇f = ∇f -             (uᵢ* (∂H/∂m)* vᵢ )
      %-------------------------------------------------
        grad = grad - H.cg2pg( real(conj(Uk).*Hmvp(dH,H.idx,Vk)) );
      end
  end
  
  freq_counter = freq_counter + 1;
end

if output_d_pred
  d_pred = reshape(d_pred,acq.nd);
end

clear A, H;
end