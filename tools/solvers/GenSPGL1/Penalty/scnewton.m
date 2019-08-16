function [x] = scnewton(fh, x0, tol, para)
% Newton's method for problems with a strictly diagonal Hessian.
% 
% function [x] = newton(func, x0)
%
% Input:  func  - function handle returning [value, gradient, hessian]
%                 where hessian is a vector containing the diagonal of the
%                 Hessian.
%         x0    - initial value of parameter
%         tol   - tolerance

% initialization
maxIter = 20;
x = x0;
iter = 0;
[f, g, h] = fh(x, para);

%fprintf('%1.5e, %1.5e\n', f, norm(g));  
while(norm(g) > tol)&&(iter < maxIter)
    
    % search direction
    d = -g./h;
    
    % Armijo linesearch
    step = 1;
    mult = 0.8;
    c = 0.01;
    lsiter = 0;
    while( fh(x + step*d, para) > f - c*step*abs(g'*d) || x + step*d < 0)
        step = step*mult;
        lsiter = lsiter + 1;
        if lsiter > 20
            return;
        end
    end
    
    % update
    x = x + step*d;
    iter = iter + 1;
    [f, g, h] = fh(x);

%    fprintf(1,'%d: %1.5e, %1.5e, %1.5e\n', iter, f, norm(g), step);  

 
end  % while