classdef opBorn < opSpot
    %opBorn SPOT wrapper for the linearized born modeling and migration
    % function "Born.m"
    %
    %   opJ = opBorn(m0,model,q,dens,ani)
    %
    %
    % INPUT for constructor:
    %
    %   m0: squared slowness [s^2/km^2]
    %   model: structure containing model parameters
    %   q: source wavelet function
    %   dens: density (optional)
    %   ani: structure with anisotropy parameters (optional)
    %
    %
    % FORWARD mode: 
    %
    %   dD = opJ*dm
    %
    % ADJOINT mode:
    %
    %   dm = opJ^T*dD
    %
    %
    %   dD: linearized (single scattered) seismic data
    %   dm: model perturbation/migrated image
    %
    %
    % Author: Philipp Witte
    %
    %
    
    properties
        
        v=[];
        model=[];
        q=[];
        thomsen=[];
        dens=[];
        
    end
    
    methods
        %% Constructor
        function op=opBorn(v,model,q,dens,thomsen)
            
            nrec = length(model.xrec);
            nsrc = length(model.xsrc);
            ns = length(model.NyqT);
            
            % Operator
            op=op@opSpot('opj',ns*nrec*nsrc,length(v(:)));
            
            % Model,velocity,density, anisotropy and data
            op.model=model;
            op.v=v;
            op.q=q;
            
            if nargin==4
            	op.dens=dens;
            elseif nargin == 5
                op.thomsen = thomsen;
            end
            
        end
    end
    
    methods ( Access = protected )
        %Multiplication function
        function y=multiply(op,x,mode)
            
            if mode==2  
                % adjoint mode: RTM
                y = Born(op.v,op.model,op.q,x,-1);

            else
                % forward mode: linearized modeling
                y = Born(op.v,op.model,op.q,x,1);

            end
        end
    end
    
end
        
         