function [f1, f2] = NLfunForward_parallel(x, g, params)

%------------------------------------------------------------------
% Use:
%   NLfunForward(x, g, params)

% Author: Haneet Wason and Rajiv Kumar
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmospheric Sciences
%         The University of British Columbia
%         
% Date: March, 2015

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%------------------------------------------------------------------

% separate L1 R1 L2 R2
e1 = params.mhnumr*params.rk;
e2 = params.mhnumc*params.rk;

L1 = x(1:e1);
R1 = x(e1+1:e1+e2);
L1 = reshape(L1,params.mhnumr,params.rk);
R1 = reshape(R1,params.mhnumc,params.rk);

L2 = x(e1+e2+1:2*e1+e2);
R2 = x(2*e1+e2+1:end);
L2 = reshape(L2,params.mhnumr,params.rk);
R2 = reshape(R2,params.mhnumc,params.rk);

if isempty(g)
    f1 = shiftzs_HSS_MH([vec(L1*R1');vec(L2*R2')],params,1);
    f2 = 0;
else
    fp1 = shiftzs_HSS_MH(g,params,-1);
    fps1 = reshape(fp1(1:params.mhnumr*params.mhnumc),params.mhnumr,params.mhnumc);
    fps2 = reshape(fp1(params.mhnumr*params.mhnumc+1:end),params.mhnumr,params.mhnumc);
    f1 = [vec(fps1*R1); vec(fps1'*L1);vec(fps2*R2); vec(fps2'*L2)];  % gradient w.r.t. L1, R1, L2, R2
    f2 = [vec(fp1)];
end

end % function

