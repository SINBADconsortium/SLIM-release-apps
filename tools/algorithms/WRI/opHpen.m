classdef opHpen < opSpot
% SPOT wrapper for the WRI Hessian
%
% Usage:
%   H = opHGNpen(m,Q,D,model,params);
%
% Author: Curt Da Silva
% 
% March 2015
%
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
    
    properties ( SetAccess = protected )
        mt,D,Q,model,params;
    end
    
    methods
        function op = opHpen(m,Q,D,model,params)
           op = op@opSpot('WRI Hessian',length(m),length(m));
           op.cflag = 0;
           op.linear = 1;
           op.children = [];
           op.sweepflag = 0;
           op.mt = m;
           op.Q = Q;
           nsrc = size(Q,2); nrec = length(model.xrec)*length(model.zrec); nfreq = length(model.freq);           
           op.D = reshape(D,nrec,nsrc,nfreq);
           op.model = model;
           params.wri = true;
           op.params = params;
           
        end
    end
    
    methods ( Access = protected )
        function y = multiply(op,x,~)            
            y = PDEfunc(PDEopts.HESS,op.mt,op.Q,x,op.D,op.model,op.params);            
        end
    end
end