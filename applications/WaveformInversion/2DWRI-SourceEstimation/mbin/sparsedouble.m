function [X] = sparsedouble(A)

n = size(A);

if n(1) > n(2)
        X = sparse(n(1),n(2));
        x = sparse(n(2),1);
  
        for i = 1:n(2)
                a      = x;
                a(i)   = 1;
                b      = A * a;
                X(:,i) = b;
        end
else
        X = sparse(n(2),n(1));
        x = sparse(n(1),1);
  
        for i = 1:n(1)
                a      = x;
                a(i)   = 1;
                b      = A' * a;
                X(:,i) = b;
        end
        X = X';
end
    
