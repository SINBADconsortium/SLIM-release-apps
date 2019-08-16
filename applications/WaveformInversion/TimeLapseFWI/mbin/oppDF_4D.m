classdef oppDF_4D < oppSpot
% pSPOT wrapper for DF.m
%
% use:
%   J = oppDF_4D(m,model,{gather})
%
% see DF_4D.m for further documentation
%

% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2012
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
       function op = oppDF_4D(mt,Q,model,dogather)
            nsrc  = size(Q,2);
            nrec  = length(model.xrec)*length(model.zrec);
            nfreq = length(model.freq);
            m = nsrc*nrec*nfreq;
            n = length(mt);
            if nargin < 4
                dogather = 0;
            end
           
           op = op@oppSpot('oppDF_4D', m, n);
           op.cflag     = 1;  
           op.linear    = 1;
           op.children  = []; 
           op.sweepflag = 1;
           op.mt        = mt;
           op.Q         = Q;
           op.model     = model;
           op.nfreq     = nfreq;
           op.nt        = nsrc*nrec;
           op.gather    = dogather;
       end 
       
    end
    
    
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = multiply(op,x,mode)
           if mode == 1
                y = DF_4D(op.mt,op.Q,x, 1,op.model,op.gather);
           else %adjoint
                y = DF_4D(op.mt,op.Q,x,-1,op.model,op.gather);  
           end
       end %multiply
       
    end %protected methods
    
end %classdef

    
    
    
