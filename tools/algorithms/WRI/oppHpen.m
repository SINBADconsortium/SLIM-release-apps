classdef oppHpen < oppSpot
    % pSPOT wrapper for the WRI Hessian.
    %
    % Usage:
    %   H = oppHpen(m,Q,D,model,params);
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
    % software.
    %
    
    properties ( SetAccess = protected )
        mt,D,Q,model,params;
    end
    
    methods
        function op = oppHpen(m,Q,D,model,params)
            op = op@oppSpot('WRI Hessian distributed',length(m),length(m));
            op.cflag = 0;
            op.linear = 1;
            op.children = [];
            op.sweepflag = 0;
            op.mt = m;
            op.Q = Q;          
            op.D = D;
            op.model = model;
            params.wri_mode = true;
            op.params = params;            
        end
    end
    
    methods ( Access = protected )
        function y = multiply(op,x,~)
            y = PDEfunc_dist(PDEopts.HESS,op.mt,op.Q,x,op.D,op.model,op.params);            
        end
    end
end