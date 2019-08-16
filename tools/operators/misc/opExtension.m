classdef opExtension < opSpot
% Extension operator. Pads input with constant values or zeros.
%
% use:
%   op = opExtension(n,nb,flag)
%
% input:
%   n    - length of input
%   nb   - number of points to add left nb(1) and right nb(2). If
%          length(nb)=1, use the same number of each side.
%   flag - 0: padd with zeros, 
%          1: padd with boundary value, 
%          2: periodic padding (not implemented)
%          3: symmetric padding 
%

% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2012
%
% Updated by Curt Da Silva, 2017
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
       nb;
       extension_mode;
    end % Properties


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opExtension(n,nb,extension_mode)
          if length(nb)==1
              nb = [nb nb];
          end
          nb = nb(1:2);
          m  = n + sum(nb);
          if nargin<3
              extension_mode = 1;
          end
          
          % Construct operator
          op = op@opSpot('Extension', m, n);
          op.nb = nb;
          op.sweepflag = 1;
          op.extension_mode = extension_mode;
       end % Constructor             
       

    end % Methods
       
 
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = multiply(op,x,mode)
          y = opExtension_intrnl(op.n,op.nb,op.extension_mode,x,mode);
       end % Multiply          

    end % Methods
   
end % Classdef


%=======================================================================


function y = opExtension_intrnl(n,nb,extension_mode,x,mode)
    nx = size(x,2);
    switch extension_mode
        case 0
            if mode==1
                y = [zeros(nb(1),nx);x;zeros(nb(2),nx)];
            else
                 y = x(nb(1)+1:end-nb(2),:);
            end
        case 1
            if mode==1
                y = [repmat(x(1,:),nb(1),1);x;repmat(x(end,:),nb(2),1)];
            else
                y = x(nb(1)+1:end-nb(2),:);
                y(1,:)   = y(1,:) + sum(x(1:nb(1),:),1);
                y(end,:) = y(end,:) + sum(x(end-nb(2)+1:end,:),1);
            end
        case 2
            if mode == 1
                
            else
                
            end
            
        case 3            
            if mode==1
                y = [flip(x(1:nb(1),:)); x; flip(x(end-nb(2)+1:end,:))];                                
            else
                y = x(nb(1)+1:end-nb(2),:);
                y(1:nb(1),:)   = y(1:nb(1),:) + x(nb(1):-1:1,:);
                y(end-nb(1)+1:end,:) = y(end-nb(1)+1:end,:) + x(end:-1:end-nb(1)+1,:);
            end            
    end
            
   
end

