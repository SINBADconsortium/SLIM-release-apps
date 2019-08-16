function [f g Hout fu] = misfit_Harm_pen(m,D,Pr,nt,dt,omega,lambda)
    %% function to calculate the misfit function, gradient and hessian for penalty method of harmonic oscillator
    % Usage:
    % [f g H ] = misfit_Harm(m,D,Pr,nt,dt,omega,lambda)
    %
    % Input:
    % m       - model parameter
    % D       - Data
    % Pr      - Projector operator for receiver
    % nt      - number of time 
    % dt      - time interval
    % omega   - frequency
    % lambda  - penalty parameter
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
    
    P      = kron(speye(nt), Pr);
    A      = [lambda*H; P];
    S      = [lambda*ss; D];
%    u      = H\ss;
    u      = A\S;
    
    Hu     = full(A'*A);
    [Uu Su Vu] = svd(Hu);
    log_Su     = log(1./diag(Su))/2;
    flog_Su    = sum(log_Su);
    
    
    Dcal   = P * u;
    % dHu    = zeros(2*nt,2);
    % for i = 1:nt-1
    %     dHu(2+2*i-1,:) = u(2*i-1:2*i)';
    % end
    
    J      = Test_GenDHu(m,nt,dt,u);
    
    v      = H*u-ss;
    f      = 0.5 * (norm(Dcal - D))^2 + 0.5*lambda^2*norm(v)^2;
    g      = J' * v * lambda^2;
    
    B      = H' \ (P');
    C      = B' * B;
    C      = eye(size(C,1)) + 1/lambda^2*C;
    C      = inv(C);
    
    
    Hout      = J' * (B * C * B') * J;
    fu        = f - flog_Su;
    
    