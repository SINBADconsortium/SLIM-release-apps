function [y dy ddy] = slowsquare2vel(m)
    %% function to convert slowness square to velocity
    % Usage:
    % [y J H] = slowsquare2vel(m)
    % 
    % Input:
    % m - slowness square [s^2/km^2]
    % Output
    % y   - velocity
    % dy  - first order derivative
    % ddy - second order derivative
    %
    %  Author:
    %  Zhilong Fang, SLIM, UBC
    %  2016/01
    
    y   = m.^(-.5);
    dy  = -0.5 * m.^(-1.5);
    ddy = 1.5*0.5 * m.^(-2.5);