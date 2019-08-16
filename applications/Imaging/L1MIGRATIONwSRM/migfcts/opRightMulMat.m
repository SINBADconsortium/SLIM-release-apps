function op = opRightMulMat(sizeA,W)
% Syntax:
% op = opRightMulMat(sizeA,W)
%
% Description:
% perform matrix multiplication on the right hand side, i.e., op*A(:) = vec(A*W)
% 
% Note:
% It is equivalent to opKron(W.',I) but faster.
%
% Input list:
% sizeA: size of A, a 2x1 vector
% W: the matrix that you multiply
% 
% Output list:
% op: the operator
% 
% Author: Ning Tu
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: Feb/14/2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

sizeW = size(W);
if sizeW(1) ~= sizeA(2)
	error('Fatal: matrix dimension mismatch.')
end
n = sizeA(1)*sizeA(2);
m = sizeA(1)*sizeW(2);

funchandle = @(x,mode) opRightMulMat_intrnl(x,mode);
op=opFunction(m,n,funchandle);

    function y = opRightMulMat_intrnl(x,mode)
        if (mode == 0)
            y = {m,n,[0,1,0,1],{'Dirac'}};
        elseif mode == 1
        	x = reshape(x, sizeA);
            y = x*W;
            y = y(:);
        elseif mode ==2
        	x = reshape(x, [sizeA(1) sizeW(2)]);
        	y = x*W';
        	y = y(:);
        end
    end
end