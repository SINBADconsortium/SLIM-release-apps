classdef opTopMute < opSpot
%opTopMute Automatic top mute of marine seismic shot record(s). Mutes the 
% direct wave and everything outside of the wedge defined by the direct
% wave.
%
% T = opTopMute(model)
%
% INPUT:
%   model:    Structure with modeling parameters
%
% USAGE:
% The operator is self-adjoint, i.e. forward and adjoint mode are the same.
%
%   Dout = T*Din
%
% Din: Vectorized input shot record(s) of length ns x nrec x nsrc
% Dout: Muted shot record(s)
%
%
% Author: Philipp Witte
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%         
% Date: February, 2016

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

properties
        model=[];
end
    
methods
    %% Constructor
    function op=opTopMute(model)

        nrec = length(model.xrec);
        nsrc = length(model.xsrc);
        ns = length(model.NyqT);

        % Operator
        op=op@opSpot('opj',ns*nrec*nsrc,ns*nrec*nsrc);

        % Model,velocity,density, anisotropy and data
        op.model=model;

    end
end
    
methods ( Access = protected )
        %Multiplication function
	function y=multiply(op,x,~)    
        
        % Get model paramters
        nShot = length(op.model.xsrc);
        nt = length(op.model.NyqT);
        nrec = length(op.model.xrec);
        dRec = op.model.xrec(2)-op.model.xrec(1);

        % Reshape inpute data
        Din = reshape(x,nt,nrec,nShot);

        % Smoothing function
        S  = opKron(opSmooth(nrec,100),opSmooth(nt,100));

        % Width of muting window
        offsetDirectWave = 1.5*op.model.T;
        idxOffset = offsetDirectWave/dRec;
        dx = round(idxOffset-idxOffset/10);

        % Loop over shots
        for i=1:nShot

            % find index of current source
            x0 = round(op.model.xsrc(i)/dRec);
            if (x0 > nrec); x0=nrec; elseif(x0<1); x0=1; end

            % use first shot to find time sample where mute starts + slope
            if i==1
                zmax = find(abs(Din(:,x0,i))==max(abs(Din(:,x0,i))));
                z0 = round(zmax-zmax/3);
                slope = (nt-zmax)/dx;
            end

            mask = ones(nt,nrec);
            mask(1:z0,:)=0;

            % left side muting window
            if x0 < dx
                x = 1;
                zIntercept = round(z0+slope*(x0-x));
                zax = z0+1:1:zIntercept;
            else
                x = x0-dx;
                zax = z0+1:1:nt;
            end
            xax = round(linspace(x0,x,length(zax)));

            for j=1:length(zax)
               mask(zax(j),1:xax(j)) = 0;
            end

            % right side muting window
            if (nrec-x0 < dx)
                x = nrec;
                zIntercept = round(z0+slope*(x-x0));
                zax = z0+1:1:zIntercept;
            else
                x = x0+dx;
                zax = z0+1:1:nt;
            end
            xax = round(linspace(x0,x,length(zax)));

            for j=1:length(zax)
               mask(zax(j),xax(j):end) = 0;
            end

            % Smooth and apply mask
            mask = reshape(S*mask(:),nt,nrec);
            Din(:,:,i) = Din(:,:,i).*mask;
        end
        y = vec(Din);
    end
end
end
