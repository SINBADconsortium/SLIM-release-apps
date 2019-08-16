% PML_FUNC - Quadratic pml function 
%             C ( (x-eta)/eta )^2           if x in bottom pml layer
% sigma(x) =  0                             if x in interior layer
%             C ( (x-eta+1)/eta )^2         if x in top pml layer
% 
% func(x) = 1/(1+1i*sigma(x)/omega);
% 
% func(x) should be evaluated on the grid 1:nx
%
% Curt Da Silva, 2015
% 
% Usage:
%   func = pml_func(nx,omega,np_bot,np_top,C);
%
% Input:
%   nx       - number of points (including pml) on grid
%   omega    - temporal frequency (Hz)
%   np_bot   - number of pml points on bottom of the domain (can be 0 - no pml)
%   np_top   - number of pml points on top of the domain (can be 0 - no pml)
%   C        - pml constant (default: 10)
% 
% Output:
%   pml_func - pml function handle
% 

function func = pml_func(nx,omega,np_bot,np_top,C)  
    if exist('C','var')==0
        C = 10;
    end
    assert( np_bot > 0 || np_top > 0);
    if np_bot ~= np_top && min(np_bot,np_top) > 0
        M = max(np_bot,np_top);
        np_bot = M; np_top = M;        
    end
    
    % number of points in the non-pml domain
    nx = nx - np_bot - np_top;
    if nx==0 %only pml
        eta = 1/2;
    elseif np_bot==0 && np_top > 0
        eta = 1/(2 + nx/(np_top-1));
    else       
        eta = 1/(2 + nx/(np_bot-1));
    end
    
    % mapping [1:nx] -> [0,1]
    x_norm = @(x) (x-1) / (nx+np_bot+np_top-1) ;
    
    if np_bot > 0 , lower = @(x) C/eta * ( (x_norm(x)- eta)/eta  ).^2; else lower = @(x) 0*x; end
    if np_top > 0 , top = @(x) C/eta * ( ( x_norm(x)- 1 + eta )/eta ).^2; else top = @(x) 0*x; end
    
    sigma = @(x) lower(x) .* (x_norm(x) <= eta ) + top(x) .* (x_norm(x) >= 1-eta);
   
    func = @(x) (1 + 1i * sigma(x)/omega).^(-1);
    
end 
