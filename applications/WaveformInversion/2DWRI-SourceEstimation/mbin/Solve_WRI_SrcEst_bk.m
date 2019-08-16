function [U, Alpha] = Solve_WRI_SrcEst(A, P, Q, D, lambda);
%%  function [U, Alpha] = Solve_WRI_SrcEst(A, P, Q, D, lambda);
% Function to solve the WRI augmented system with source estimcation
% Input : 
% A : Helmholtz Matrix
% P: Receiver projection matrix
% Q: Source
% D: Data
% lambda: Penalty parameter
%
% Output
% U : wavefield
% Alpha: source wavelet
%
% Author: Zhilong Fang
% Date: Aug. 2016

LAMBDA = lambda^2;
B = LAMBDA * A'*A + P'*P;

%R = chol(B);

[LL,UU,Pp,Qp,R] = lu(B);
Hinv = @(x) Qp*(UU\(LL\(Pp*(R\(x))))); 
Hinvadj = @(x) R'\(Pp'*(LL'\(UU'\(Qp'*x))));

AQ = LAMBDA*A'*Q;
Pd  = P'*D;

Y     = Hinv([AQ Pd]);

%Y     = R\(R'\([AQ Pd]));
 
 YAQ = Y(:,1:size(AQ,2));
 YPd  = Y(:,size(AQ,2)+1:end);
 
 for i = 1:size(AQ,2)
        q               = Q(:,i);
        d               = D(:,i);
        e               = AQ(:,i);
        a               = LAMBDA*q'*q - e'*YAQ(:,i);
        b               = -e'*YPd(:,i);
        Alpha(i) = b/a;
        S(:,i)         = Pd(:,i) - e*Alpha(i);
 end
 
 %U = R\(R'\S);

 U = Hinv(S);

 Alpha = - Alpha(:);

