classdef opDF < opSpot
% SPOT wrapper for the Jacobian of F.m
%
% use:
%   J = opDF(m,Q,model)
%
% see PDEfunc.m for further documentation
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
        mt,Q,model,nfreq,nt,params;
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opDF(mt,Q,model,params)
           nsrc  = size(Q,2);
           nrec  = length(model.xrec)*length(model.zrec);
           nfreq = length(model.freq);
           m = nsrc*nrec*nfreq;
           n = length(mt);           
           
           op = op@opSpot('opDF', m, n);
           op.cflag     = 1;
           op.linear    = 1;
           op.children  = [];
           op.sweepflag = 0;
           op.mt        = mt;
           op.Q         = Q;
           op.model     = model;
           op.nfreq     = nfreq;
           op.nt        = nsrc*nrec;
           if exist('params','var')==0 || isempty(params)
               params = struct;
           end
           params.wri = false;
           params = default_fwi_params2d(params);
           op.params = params;
       end        
    end
        
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = multiply(op,x,mode)
           mt = op.mt; Q = op.Q; model = op.model; params = op.params;            
           nsrc = size(Q,2); nrec = length(model.xrec)*length(model.zrec); nfreq = length(model.freq);
           if mode == 1
               y = PDEfunc(PDEopts.JACOB_FORW,mt,Q,x,[],model,params);
           else
               y = PDEfunc(PDEopts.JACOB_ADJ,mt,Q,x,[],model,params);
           end   
           y = vec(y);
       end %multiply
       
    end %protected methods
    
end %classdef

    
    
    
