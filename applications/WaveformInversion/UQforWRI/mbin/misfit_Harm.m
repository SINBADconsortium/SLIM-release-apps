function [f g H ] = misfit_Harm(m,D,Pr,nt,dt,omega)
    %% function to calculate the misfit function, gradient and hessian for harmonic oscillator
    % Usage:
    % [f g H ] = misfit_Harm(m,D,Pr,nt,dt,omega)
    %
    % Input:
    % m       - model parameter
    % D       - Data
    % Pr      - Projector operator for receiver
    % nt      - number of time 
    % dt      - time interval
    % omega   - frequency
    %
    % Output
    % f       - misfit function
    % g       - gradient
    % H       - Hessian
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
    Dcal   = P * u;
    % dHu    = zeros(2*nt,2);
    % for i = 1:nt-1
    %     dHu(2+2*i-1,:) = u(2*i-1:2*i)';
    % end
    
    J      = Test_GenDHu(m,nt,dt,u);
    
    f      = 0.5 * (norm(Dcal - D))^2;
%    g      = dH' * (dHu' * (H'\(P'*(Dcal-D))));
    g      = -J' * (H'\(P'*(Dcal-D)));
    
    A      = (J' * (H' \ P'));
    H      = A*A';