classdef oppDF < oppSpot
% pSPOT wrapper for DF.m
%
% use:
%   J = oppDF(m,Q,model)
%
% see DF.m for further documentation
%

% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2012
% 
% Modified by: Curt Da Silva
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        mt,Q,model,nt,params;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = oppDF(mt,Q,model,params)
           nsrc  = size(Q,2);
           nrec  = length(model.xrec)*length(model.zrec);
           nfreq = length(model.freq);
           m = nsrc*nrec*nfreq;
           n = length(mt);
           op = op@oppSpot('oppDF', m, n);
           op.cflag     = 1;
           op.linear    = 1;
           op.children  = [];
           op.sweepflag = 0;
           op.mt        = mt;
           op.Q         = Q;
           model = fwi2d_model_compatibility(model);
           op.model     = model;
           op.nt        = nsrc*nrec;
           if exist('params','var')==0 || isempty(params)
               params = struct;
           end
           params       = default_fwi_params2d(params);
           params.wri   = false;
           op.params    = params;
       end 
       
    end
    
    
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = multiply(op,x,mode)
           if mode==1
               y = PDEfunc_dist(PDEopts.JACOB_FORW,op.mt,op.Q,x,[],op.model,op.params);
               y = pSPOT.utils.distVectorize(y);
           else
               y = PDEfunc_dist(PDEopts.JACOB_ADJ,op.mt,op.Q,x,[],op.model,op.params);
               y = vec(gather(y));
           end                      
       end %multiply       
    end %protected methods
    
end %classdef

    
    
    
