classdef oppHGNpen < oppSpot
    % pSPOT wrapper for WRI Gauss-Newton Hessian.
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
        function op = oppHGNpen(m,Q,D,model,params)
            op = op@oppSpot('WRI GN Hessian',length(m),length(m));
            op.cflag = 0;
            op.linear = 1;
            op.children = [];
            op.sweepflag = 0;
            op.mt = m;
            op.Q = Q;          
            nsrc = size(Q,2); nrec = length(model.xrec); nfreq = length(model.freq);
            if size(D,2) ~= nfreq && nfreq > 1
                D = pSPOT.utils.distVec2distArray(pSPOT.utils.distVectorize(D),[nrec*nsrc,nfreq]);
            end
            op.D = D;
            op.model = model;
            op.params = params;            
        end
    end
    
    methods ( Access = protected )
        function y = multiply(op,x,~)
            y = PDEfunc_dist(PDEopts.HESS_GN,op.mt,op.Q,x,op.D,op.model,op.params);            
        end
    end
end