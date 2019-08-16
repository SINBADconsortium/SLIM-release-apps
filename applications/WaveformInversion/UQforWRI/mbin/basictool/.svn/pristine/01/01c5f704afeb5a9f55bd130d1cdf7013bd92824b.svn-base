function [Y] = LowPassFilter(X, CutFreq, CutWidth,  dFreq)


n 	= 	size(X,1); 
F	= 	opDFT(n);
XX	=	F * X;

k1      =  floor(CutFreq/dFreq);
k2      =  floor((CutFreq+CutWidth) / dFreq);

b       = ones(n,1);
b(1:k1)         = 1;
b(k2:end)   = 0;
z                     = [0:k2-k1] / (k2-k1) * pi/2;
b(k1:k2)      = cos(z);
FF                  = opDiag(b);

%Y                      = FF * X;
Y                      = real(F' * (FF*XX));


