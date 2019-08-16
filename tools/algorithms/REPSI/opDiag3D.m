function op = opDiag3D(n1,n2,n3,d)
% OPDIAG3D          An optimized version of opDiag acting on the first axis of a datacube
%                   (Multiply each 2D slice with corresponding element of a vector d)

if (n1 ~= length(d))
    error('OPDIAG3D: vector length mismatch to n1');
end
    
d = d(:);

subfunc_handle = @(x,mode) opDiag3D_intrnl(n1,n2,n3,d,x,mode);
m = n1*n2*n3;
n = n1*n2*n3;

op = opFunction(m,n,subfunc_handle); % return a SPOT operator using constructor opFunction

function y = opDiag3D_intrnl(n1,n2,n3,d,x,mode)

switch mode
    case 0
        y = {m,n,[0,1,0,1],{'opDiag3D'}};
   
    case 1
        x = reshape(x,n1,n2,n3);
        x = permute(x, [2 3 1]);
        y = x .* repmat(reshape(d, [1 1 n1]), [n2 n3]);
        y = permute(y, [3 1 2]);
        y = y(:);

    case 2
        x = reshape(x,n1,n2,n3);
        x = permute(x, [2 3 1]);
        y = x .* repmat(reshape(conj(d), [1 1 n1]), [n2 n3]);
        y = permute(y, [3 1 2]);
        y = y(:);        
end
