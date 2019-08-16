function [H dH] = GenH_Harmonic_Oscillator(m,nt,dt);
    %% function to use generate the stencil of harmonic oscillator system
    % Usage:
    % [H, dH] = GenH_Harmonic_Oscillatort(m,nt,dt,omega);
    %
    % d[u1](t) = [-2m1m2 -m1^2][u1](t) + [sin(wt)]. 
    % d[u2]      [1      0    ][u2]      [0]
    %
    % Input:
    % m       - model parameter
    % nt      - number of time 
    % dt      - time interval
    % omega   - frequency
    %
    % Output
    % H       - Full stencil
    % dH      - Derivative of stencil
    %
    %  Author:
    %  Zhilong Fang, SLIM, UBC
    %  2016/01 
    
    B  = [-2*m(1)*m(2) -m(1)^2;1 0];
    II = eye(2)/dt;
    C  = -II - B;
    
    HI  = kron(speye(nt), II);
    HII = kron(spdiags(ones(nt,1),-1,nt,nt), C);
    
    H   = HI + HII;
    dH  = [2*m(2) 2*m(1); 2*m(1) 0];