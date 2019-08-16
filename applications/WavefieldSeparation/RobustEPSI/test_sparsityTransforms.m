function test_sparsityTransforms()
% Tests for existence oand operation of RWT and CurveLab operators

n = 32;

Curv = opCurvelet(n,n,2,16,1,'WRAP',0);
W = opTranspose(opWavelet(n,1));
S = opKron(Curv,W);

retcode = dottest(S,2);

if retcode == 0
    disp(' ')
    disp('TEST_SPARSITYTRANSFORMS successfully completed without issue')
end
