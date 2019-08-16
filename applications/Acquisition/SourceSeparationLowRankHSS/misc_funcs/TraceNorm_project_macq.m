function [x] = TraceNorm_project_macq(x,weights,B,params)

% Force the rows of L and R to have norm at most B.
%
% (LOut, ROut) = arg min_(U,V) || [L;R]-[U;V] ||_F^2 s.t. || [U;V] ||_mr < B,
%
% where || A ||_mr is the maximum Euclidean norm of a row of A.

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

c = sqrt(B/(0.5*norm(x)^2));
L1 = min(1,c)*L1(:);
R1 = min(1,c)*R1(:);
L2 = min(1,c)*L2(:);
R2 = min(1,c)*R2(:);
x = [vec(L1);vec(R1);vec(L2);vec(R2)];

end % function end

