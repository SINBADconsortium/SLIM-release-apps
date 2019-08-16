classdef opDecimate < opSpot
% Signal decimation with Kaiser-window low-pass filtering
%
% Curt Da Silva, 2015
%
% Usage:
%   A = opDecimate(xin,xout);
% 
% Input:
%   xin   - input grid (n x 1)
%   xout  - output grid (must be a subset of xin) (nc x 1)
% 
% Usage:
%   A = opDecimate(n,r);
%   
% Input:
%   n     - integer input grid length
%   r     - decimation factor < 1, will be truncated to k/p for k=1,2,3,4, p an integer
% 
% Output:
%   A     - decimation operator
%
    properties (SetAccess = protected)
        k,p;
    end
    
    methods
        function op = opDecimate(n,r)
            [~,den] = rat(r);
            
            maxk = 4;
            [K,P] = ndgrid(1:maxk,1:den);
            [~,I] = min(vec(abs(r-K./P)));
            [k,p] = ind2sub([maxk,den],I);
            nr = floor(k*n/p); nc = n;
            
            op = op@opSpot('Decimation operator', nr,nc);                
                        
            op.k = k; op.p = p;
            op.sweepflag = true;
        end
    end
    
    methods ( Access = protected)
        function y = multiply(op,x,mode)
            n = op.n; k = op.k; p = op.p;
            [~,fout] = spacefreq_grid(n*k,1/k);
            fnyq = max(fout);
            cutoff_freq = 0.8*fnyq/(p);
            opts = struct; opts.trans_width = 0.3*fnyq/(p);
            opts.ripple = 1e-3; opts.full = true;
            w = filter_bank(n*k,1/k,'kaiser',cutoff_freq,opts);
            if mode==1
                % Interpolate k-fold by zero-padding in space                                
                if k~= 1                                        
                    y = zeros(n*k,size(x,2));
                    y(1:k:n*k,:) = k*x;
                else
                    y = x;
                end
                % Lowpass 
                y = conv2(y,w,'same');
                
                % Decimate                
                y = y(1:p:k*n,:);
            else
                y = zeros(k*n,size(x,2));
                y(1:p:k*n,:) = x;
                y = conv(y,w([end:-1:1]),'same');
                if k~=1
                    y = y(1:k:n*k,:);
                else
                    y = x; 
                end
            end
              
            
        end
    end
end