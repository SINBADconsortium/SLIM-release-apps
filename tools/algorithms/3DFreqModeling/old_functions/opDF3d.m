classdef opDF3d < opSpot
% pSPOT wrapper for DF3d.m
%
% use:
%   J = opDF3d(m,Q,model)
%
% see DF3d.m for further documentation
%
%
% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2013
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        mt,Q,model,nfreq,nt;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opDF3d(mt,Q,model)
            nsrc  = size(Q,2);
            nrec  = length(model.xrec)*length(model.zrec)*length(model.yrec);
            nfreq = length(model.freq);
            m = nsrc*nrec*nfreq;
            n = length(mt);
           
           op = op@opSpot('opDF3d', m, n);
           op.cflag     = 1;  
           op.linear    = 1;
           op.children  = []; 
           op.sweepflag = 0;
           op.mt        = mt;
           op.Q         = Q;
           op.model     = model;
           op.nfreq     = nfreq;
           op.nt        = nsrc*nrec;
       end 
       
    end
    
    
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = multiply(op,x,mode)
           if mode == 1
                y = DF3d(op.mt,op.Q,x, 1,op.model);
           else %adjoint
                y = DF3d(op.mt,op.Q,x,-1,op.model);  
           end
       end %multiply
       
    end %protected methods
    
end %classdef

