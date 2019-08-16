function [stdm smp fsmp] = Sample_Harm_pen(m,D,Pr,nt,dt,omega,lambda,nsmp, mp, Hp, alphaD)
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
    m = m(:); 
    
    for ismp = 1:nsmp
    
        [H dH] = GenH_Harmonic_Oscillator(m,nt,dt);
        s = sin(omega*([0:nt-1]*dt));
        s = s(:);
        ss = zeros(length(s)*2,1);
        for i = 1:length(s)
            ss(2*i-1)=s(i);
        end
    
        P      = kron(speye(nt), Pr);
        A      = [lambda*H; alphaD*P];
        S      = [lambda*ss; alphaD*D];
        %    u      = H\ss;
        u      = A\S;
    
        Hu     = full(A'*A);
        L      = chol(Hu);
        u      = u + (L\randn(size(L,1),1));
    
    
        Dcal   = P * u;
        % dHu    = zeros(2*nt,2);
        % for i = 1:nt-1
        %     dHu(2+2*i-1,:) = u(2*i-1:2*i)';
        % end
    
        J      = Test_GenDHu(m,nt,dt,u);
    
        v      = H*u-ss;
        f      = 0.5 * alphaD^2 * (norm(Dcal - D))^2 + 0.5*lambda^2*norm(v)^2 + .5 * (m-mp)'*(Hp*(m-mp));
        g      = J' * v * lambda^2;
        g      = g + Hp *(m-mp);
    
        Hout      = lambda^2 * J' * J;
        Hout      = Hout + Hp;
        Lm        = chol(Hout);
        dm        = Hout \ g;
        m         = m - dm;
        m         = m + Lm\randn(size(Lm,1),1);
        smp(:,ismp)  = m;
        fsmp(ismp)   = f;
        if mod(ismp,1000) == 1
            fprintf('%2d,',ismp);
        end
    end
    
    for i = 1:size(smp,1)
        stdm(i) = std(smp(i,:));
    end
    