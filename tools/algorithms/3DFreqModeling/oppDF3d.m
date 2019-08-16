classdef oppDF3d < oppSpot
% pSPOT wrapper for 3D migration/demigration. See PDEfunc_dist.m
% for more details.
%
% Curt Da Silva, 2015
% 
% Usage:
%  J = oppDF3d(m,Q,model,params);
% 
% Input:
%   m      - 3D model vector
%   Q      - source weight matrix
%   model  - model parameters struct
%   params - performance parameters struct
%
% Output:
%   J      - nrec*nsrc*nfreq x prod(model.n) SPOT operator

    properties ( SetAccess = protected )
        mt,Q,model,params;
    end
    
    methods
        function op = oppDF3d(mt,Q,model,params)
            nsrc = size(Q,2);
            nrec = numel(model.xrec)*numel(model.yrec)*numel(model.zrec);
            nfreq = length(model.freq);
            op = op@oppSpot('oppDF3d', nsrc*nrec*nfreq,numel(mt));
            op.mt = mt; 
            op.cflag = true;
            op.sweepflag = false;
            op.Q = Q;
            op.model = model;
            op.params = params;            
        end
    end
    methods ( Access = protected )
        function y = multiply(op,x,mode)
           if mode==1
               y = PDEfunc_dist(PDEopts.JACOB_FORW,op.mt,op.Q,x,[],op.model,op.params);
               y = pSPOT.utils.distVectorize(y);
           else              
               nsrc = size(op.Q,2);
               nrec = numel(op.model.xrec)*numel(op.model.yrec)*numel(op.model.zrec);
               nfreq = length(op.model.freq);
               x = pSPOT.utils.distVec2distArray(x,[nrec,nsrc*nfreq]);
               y = PDEfunc_dist(PDEopts.JACOB_ADJ,op.mt,op.Q,x,[],op.model,op.params);
           end
        end
    end
end