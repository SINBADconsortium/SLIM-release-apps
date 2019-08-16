function test_suite = test_transpose
test_suite=buildFunctionHandleTestSuite(localfunctions);
end

function test_transpose_elementary_ops
    m = 2; n = 2;
    check_2Input('opBernoulli',m,n);
    check_2Input('opBinary',m,n);
    check_1Input('opDCT',m);
    check_2Input('opDCT2',m*n,m*n);
    check_2Input('opDFT',m,n);
    check_2Input('opDFT2',m*n,m*n);
    check_1Input('opDirac',m);
    check_2Input('opEmpty',m,0);
    check_4Input('opExtend',m,n,2*m,2*n);
    check_2Input('opEye',m,m);
    check_2Input('opEye',m,n);
    check_2Input('opGaussian',m,n)
    check_1Input('opHaar',m*2^5);
    check_2Input('opHaar2',64,32);
    check_2Input('opHadamard',n,n);
    check_2Input('opHeaviside',n,n);
    check_2Input('opOnes',m,n);
    check_2Input('opZeros',m,n);
end

function check_1Input(operatorName,m)
    A = eval(strcat(operatorName,'(',num2str(m),')'));
    checkEqual(A);
end
function check_2Input(operatorName,m,n)
    A = eval(strcat(operatorName,'(',num2str(m),',',num2str(n),')'));
    checkEqual(A);
end
function check_4Input(operatorName,v1,v2,v3,v4)
    A = eval(strcat(operatorName,'(',num2str(v1),',',num2str(v2),',',num2str(v3),',',num2str(v4),')'));
    checkEqual(A);
end
function checkEqual(A)
    assertEqual(sparse(real(double(A.'))),sparse(real(double(A'))));
    assertEqual(sparse(imag(double(A.'))),sparse(imag(-double(A'))));
end