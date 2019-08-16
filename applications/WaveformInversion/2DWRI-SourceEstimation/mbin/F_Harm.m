function D = F_Harm(m,Pr,nt,dt,omega)
    %% function to perform the forward modeling of Harmoic oscillator
    % Usage:
    % D = F_Harm(m,Pr,nt,dt,omega)
    %
    % Input:
    % m       - model parameter
    % Pr      - Projector operator for receiver
    % nt      - number of time 
    % dt      - time interval
    % omega   - frequency
    %
    % Output
    % D       - observed data
    %
    %  Author:
    %  Zhilong Fang, SLIM, UBC
    %  2016/01 
    
    [H dH] = GenH_Harmonic_Oscillator(m,nt,dt);
    s = sin(omega*([0:nt-1]*dt));
    s = s(:);
    ss = zeros(length(s)*2,1);
    for i = 1:length(s)
        ss(2*i-1)=s(i);
    end
    
    u      = H\ss;
    P      = kron(speye(nt), Pr);
    D      = P * u;