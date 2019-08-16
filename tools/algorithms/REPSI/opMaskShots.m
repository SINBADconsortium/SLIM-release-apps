classdef opMaskShots < opSpot
%OPMASK  Selection mask.
%   
%   opMaskShots(NR,NS,IDX) creates an operator that takes the input x,
%   reshapes it into a NR-by-NS 2D array, then masks (set to zero) the
%   columns specified by IDX, which is a vector of column indices. 
%   When x represents a time/frequency slice of a seismic line prestack
%   datacube in the 'canonical' representaiton, this operator represents
%   masking of shots specified by the index list IDX
%   
%

%   Author      : Tim Lin
%                 Seismic Laboratory for Imaging and Modelling
%                 Department of Earth & Ocean Sciences
%                 The University of British Columbia
%         
%   Date        : Faburary, 2011



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
        maskShotIdx; % Binary mask
        nr;
        ns;
    end % Properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - Public
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opMaskShots(varargin)
          
            if nargin == 3
                nr = varargin{1};
                ns = varargin{2};
                idx = varargin{3};
            else
                error('Invalid number of parameters specified.')
            end

 
            if spot.utils.isposintmat(idx) || isempty(idx)
                if ~isempty(idx) && (max(idx) > ns)
                    error('Index exceeds operator dimensions.');
                end
            else
                error('Subscript indices must either be real positive integers or logicals.');
            end
           
           
           % Construct operator
           op = op@opSpot('MaskShots',nr*ns,nr*ns);
           op.maskShotIdx = idx;
           op.nr = nr;
           op.ns = ns;
        end % Constructor
        
    end % Methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - protected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = protected )
   
        % Multiplication
        function y = multiply(op,x,mode)
            x = reshape(x,[op.nr op.ns]);
            x(:,op.maskShotIdx) = 0;
            x = x(:);
            y = x;
        end % Multiply
  
    end % Methods
        
end % Classdef
