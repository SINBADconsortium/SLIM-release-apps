function [U, Alpha] = Solve_WRI_SrcEst2(A, P, Q, D, lambda);

nsrc = size(Q,2);

for i = 1:nsrc
        b = zeros(size(P,1),1);
        B = [lambda*A  -lambda * Q(:,i); ...
                P b];
        d = [zeros(size(A,1),1); D(:,i)];
        x = (B'*B) \ (B'*d);
        U(:,i) = x(1:end-1);
        Alpha(i) = x(end); 
end
