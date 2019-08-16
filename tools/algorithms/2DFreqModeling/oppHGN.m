classdef oppHGN < oppSpot
% pSPOT wrapper for the Gauss Newton Hessian of F.m
%
% Usage:
%   H = oppHGN(m,Q,D,model,{params})
%
% see PDEfunc.m for further documentation
%

% Author: Curt Da Silva
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: March, 2015
% 
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        mt,Q,model,params;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = oppHGN(mt,Q,model,params)
           m = numel(mt);
           n = numel(mt);
           
           op = op@oppSpot('FWI GN Hessian', m, n);
           op.cflag     = 0;  
           op.linear    = 1;
           op.children  = []; 
           op.sweepflag = 0;
           op.mt        = mt;
           op.Q         = Q;           
           
           op.model     = model;
           if exist('params','var')==0||isempty(params)           
               params = struct; 
           end
           params = default_fwi_params2d(params);
           params.wri = false;
           op.params = params;
           
       end        
    end
    
    
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = multiply(op,x,~)           
           y = PDEfunc_dist(PDEopts.HESS_GN,op.mt,op.Q,x,[],op.model,op.params);
       end %multiply
       
    end %protected methods
    
end %classdef

    
    
    
