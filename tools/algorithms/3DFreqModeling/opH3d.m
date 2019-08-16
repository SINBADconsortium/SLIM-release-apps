classdef opH3d < opSpot
% SPOT wrapper for the full Hessian of F3d.m
%
% use:
%   H = opH3d(m,Q,D,model,params)
%
% see PDEfunc3D.m for further documentation
%

% Author: Curt Da Silva, 2015
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
        mt,Q,D,model,params;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opH3d(mt,Q,D,model,params)            
           n = numel(mt);            
                           
           op = op@opSpot('opH3d', n,n);
           op.cflag     = 0;  
           op.linear    = 1;
           op.sweepflag = 0;
           op.mt        = mt;
           op.Q         = Q;
           op.D         = D;
           op.model     = model;
           op.params = params; 
       end                      
       
    end
        
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = multiply(op,x,~)
           y = PDEfunc(PDEopts.HESS,op.mt,op.Q,x,op.D,op.model,op.params);       
       end %multiply
       
    end %protected methods
    
end %classdef

    
    
    
