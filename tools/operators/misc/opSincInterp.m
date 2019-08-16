classdef opSincInterp < opSpot
    %OPSINCINTERP - Uses sinc interpolation to interpolate a function
    %  between two grids, xin and xout. Normalizes the sinc function to the
    %  output grid, unless the output grid is just a single point. Utilizes
    %  optional Kaiser windowing.
    %
    %  Curt Da Silva, 2015
    %
    %  Usage:
    %    op = opSincInterp(xin,xout,r);
    %
    %  Input:
    %    xin    - input grid
    %    xout   - output grid
    %    r      - optional Kaiser window half length (default: 0, no Kaiser
    %             window)
    % 
    %  Output:
    %    op     - sinc interpolation SPOT operator
    
    properties
        sinc_func, window;
    end
    
    methods
        function op = opSincInterp(xin,xout,r)
            if exist('r','var')==0
                r = 0;
            end
            if length(xout)==1
                % Convert xin and xout grids normalized to xin units
                
            else                
                % Convert xin and xout grids normalized to xout units
                %min_xout = min(xout);
                %xout = xout - min_xout;
                %xin  = xin - min_xout;
                dx = xout(2)-xout(1);
                xin = xin/dx; xout = xout/dx;
            end
            [XOUT,XIN] = ndgrid(xout,xin);
            op = op@opSpot('Sinc interpolation',length(xout),length(xin));
            
            op.sinc_func = sinc(XOUT-XIN);
            if r > 0
                r_b =[1.24,2.94, 4.53, 6.31, 7.91, 9.42, 10.95, 12.53, 14.09, 14.18];
                op.window = kaiser_window(XOUT-XIN,r,r_b(r));
            else
                op.window = [];
            end
            op.sweepflag = true;
        end
    end
    
    methods (Access = protected)
        function y = multiply(op,x,mode)
            if mode==1
                if ~isempty(op.window)
                    y = (op.window .* op.sinc_func)*x;
                else
                    y = op.sinc_func*x;
                end
            else
                if ~isempty(op.window)
                    y = (op.window .* op.sinc_func)'*x;
                else
                    y = op.sinc_func'*x;
                end
            end
        end
    end
    
end

