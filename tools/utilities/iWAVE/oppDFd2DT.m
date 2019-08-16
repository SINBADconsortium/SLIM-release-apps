classdef oppDFd2DT < opSpot
    % SPOT wrapper for DFD2DT
    %
    % use:
    %    J = oppDFd2DT(m,model,options)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        mt,model,nsr,nt,options;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppDFd2DT(mt,model,options)
            if isfield(model,'simultsrc')
                ismulti = model.simultsrc;
            else
                ismulti = 0;
            end
            
            if isfield(model,'romb')
                isromb = model.romb;
            else
                isromb = 0;
            end
            
            if isfield(model,'encod')
                isphase = model.encod;
            else
                isphase = 0;
            end
            
            if isfield(model,'ysrc')
                nsrc = length(model.ysrc)*length(model.xsrc)*length(model.zsrc);
            else
                nsrc = length(model.xsrc)*length(model.zsrc);
            end
            
            if isfield(model,'yrec')
                nrec = length(model.yrec)*length(model.xrec)*length(model.zrec);
            else
                nrec = length(model.xrec)*length(model.zrec);
            end
            
            nt = length(model.t);
            if ismulti == 0 && isphase == 0 && isromb == 0;
                m = nsrc*nrec*length(model.t);
            else
%                m = options.supershots*nrec*length(model.t);
                m = size(model.srcdata,3)*nrec*length(model.t);
            end
            if options.bulkonly == 1
                n = length(mt)/2;
            else
                n = length(mt);
            end
            nsr = nsrc*nrec;
            if nargin < 4
                dogather = 0;
            end
            
            op = op@opSpot('oppDFd2DT', m, n);
            op.cflag     = 0;
            op.linear    = 1;
            op.children  = [];
            op.sweepflag = 1;
            op.mt        = mt;
            op.model     = model;
            op.nsr       = nsrc*nrec;
            op.nt        = nt;
            op.options   = options;
        end
        
    end
    
    
    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            if mode == 1
                y = DFd2DT(op.mt,x, 1,op.model,op.options);
            else %adjoint
                y = DFd2DT(op.mt,x,-1,op.model,op.options);
                model = op.model;
                dt          = 1000*(model.t(2) - model.t(1));
                y            = y / dt * prod(model.d);
            end
        end %multiply
        
    end %protected methods
    
end %classdef
