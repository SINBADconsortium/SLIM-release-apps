function [U, Alpha] = Solve_WRI_SrcEst(A, P, Q, D, lambda, flagseq);
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

if nargin < 6
  flagseq = 0;
end

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

 YAQ  = Y(:,1:size(AQ,2));
 YPd  = Y(:,size(AQ,2)+1:end);

 if flagseq > 0
    for i = 1:size(AQ,2)
        q               = Q(:,i);
        d               = D(:,i);
        e               = AQ(:,i);
        a               = LAMBDA*q'*q - e'*YAQ(:,i);
        b               = -e'*YPd(:,i);
        Alpha(i) = b/a;
        S(:,i)         = Pd(:,i) - e*Alpha(i);
        U(:,i)        = YPd(:,i) - Alpha(i) * YAQ(:,i);
    end
  else
    a = 0;
    b = 0;
    for i = 1:size(AQ,2)
        q               = Q(:,i);
        d               = D(:,i);
        e               = AQ(:,i);
        a               = a + LAMBDA*q'*q - e'*YAQ(:,i);
        b               = -e'*YPd(:,i) + b;
    end
    Alpha = b / a;
    S     = Pd - Alpha * AQ;
    U    = YPd - Alpha * YAQ;
    Alpha = Alpha * ones(size(Q,2),1);
  end

 %U = R\(R'\S);

 %U = Hinv(S);

 Alpha = - Alpha(:);
