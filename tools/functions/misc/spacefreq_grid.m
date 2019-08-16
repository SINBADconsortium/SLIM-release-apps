function [x,f] = spacefreq_grid(n,dt)
    
    if exist('dt','var')==0, dt = 1; end
    
    if mod(n,2)==0
        x = (-n/2+1:n/2)*dt;
        f = (-1/2+1/n:1/n:1/2)*(1/dt);
    else
        x = (-(n-1)/2:(n-1)/2)*dt;
        f = (-1/2+1/n:1/n:1/2)*(1/dt);
    end
end