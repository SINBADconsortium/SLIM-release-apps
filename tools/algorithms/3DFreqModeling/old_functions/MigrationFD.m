%-------------------------------------------------------------------------------
%
%                (PRE-ALPHA VERSION!!! MIGHT NOT WORK AT ALL!!)
%                
%------------------------------------------------------------------------------- 
% 3D Frequency domain finite differences modeling operator. 
% Given a model and a smooth background model, this operator computes the 
% model perturbation as in
% 
% δm = ∇F(m0)⁺(F(m) - F(m0))
% 
% where m0 is the background velocity model, and m is the velocity model.
% 
% This is nothing but a wrapper to WaveEquationFD.
%
% Use: 
%   dmodel = MigrationFD(model, b_model, acq, opts, flog)
%
% Input:
%  model      - Contains information about the current model to which the 
%               forward modeling is to be applied. Read further for more details
%               on its content.
%  b_model    - Like model, but for the background velocity model.
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
%   d_pred    - data cube as an array whose dimension is (xrec,yrec,zrec,xsrc,
%               ysrc,zsrc,nfreq). This is not vectorized.
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
% 
%  acq.freq   - array with frequencies
%  acq.{zsrc,xsrc,ysrc} - vectors describing source array
%  acq.{zrec,xrec,yrec} - vectors describing receiver array.
%  acq.sources- array with dimensions (xsrc*ysrc*xsrc, nsrc) containing the 
%               source position, where nsrc is the number of sources.
% Optional:
%  acq.f0     - peak frequency of Ricker wavelet, 0 for no wavelet.
%  acq.t0     - phase shift [s] of wavelet. Default to 0.
%  
% *Optional*:
%  data       - observed data, used to compute the value of the value of the
%               misfit function and the gradient. Can be neglected if the only
%               output desired it d_pred.
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
% par_helm    - Same as par_solver, but this is passed to discrete_helmholtz
%               function.
% rng_freq    - An array containing the index of the frequencies of interest.
%               For instance, if this value is set to [1 2 5], this function 
%               will take into account only acq.freq(1), acq.freq(2) and 
%               acq.freq(5), and ignore all other possible frequencies.
% rng_src     - Same as rng_freq, but for acq.sources. 
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
%
%                (PRE-ALPHA VERSION!!! MIGHT NOT WORK AT ALL!!)
%
%-------------------------------------------------------------------------------
function dmodel = MigrationFD(model, b_model, acq, opts, flog)
   
   if ~exist('flog','var')
      flog = 0;
   end
   
   % Gets F(m)
   acq.data = ForwardFD(model,acq,opts,flog);
   
   opts.output_misfit = 0;
   opts.output_grad   = 1;
   opts.output_d_pred = 0;
   
   % Effectively computes
   %  -δm      = ∇F(m0)⁺(F(m0) - F(m))
   [~, dmodel] = WaveEquationFD(b_model, acq, opts, flog);
   
   % There you go
   dmodel = -dmodel;
end


