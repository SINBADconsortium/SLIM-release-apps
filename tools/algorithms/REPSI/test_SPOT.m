function test_sparco()
% Tests for existence of basic SLIM-SPOT operators

n = 128;
x = randn(n,1);
s = randn(n,1);

I = opDirac(n);
F = opDFT(n);
Is = opDiag(F*s);
CONV = I * F' * Is * F;


retcode = dottest(CONV);

if retcode == 0
    disp(' ')
    disp('TEST_SPOT successfully completed without issue')
end
