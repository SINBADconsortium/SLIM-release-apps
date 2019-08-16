classdef opFFTsym_conv_mask < opSpot
% Syntax:
% op = opFFTsym_mask(n, mask)
% "mask" is optional
%
% Description:
%  One-dimensional symmetric real fast Fourier transform (FFT) that
%  does not include self-adjoint scaling and does not double the
%  magnitude of the positive frequencies. WILL NOT PASS
%  DOT-TEST. ONLY meant for computing convolutions where the FFT
%  appears with a corresponding inverse. In this case the adjoint
%  (mode = 2) actually implements the inverse.
%
% Input list:
% n: length of the vector
% mask: (optional input) frequency mask, MUST be logical if present
%
% Output list:
% op: the operator
%
% Author: Ning Tu
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: Oct/17/2013
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( SetAccess = private )
        mask;
        l;
    end % Properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - Public
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        % Constructor
        function op = opFFTsym_conv_mask(n, mask)
            if mod(n,2) == 1
                l = ceil(n/2);
            else
                l = (n/2) + 1;
            end
            
            if not(exist('mask','var'))
                mask = true(l,1);
                m = l;
            else
                m = sum(double(mask));
            end
            if not(islogical(mask))
                mask = logical(mask);
            end

            % Construct operator
            op = op@opSpot('FFTsym',m,n);
            op.mask = mask;
            op.l = l;
            op.cflag     = 1;
            op.linear    = 1;
            op.sweepflag = true;
        end % Constructor
        
    end % Methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - protected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = protected )
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            if mode == 1
                y = fft(x);
                y = y(1:op.l,:);
                y = y(op.mask,:);
            else
                y = zeros(op.l,size(x,2));
                y(op.mask,:) = x;
                y = ifft(y,op.n,'symmetric');
            end
        end % Multiply
       
    end % Methods
        
end % Classdef