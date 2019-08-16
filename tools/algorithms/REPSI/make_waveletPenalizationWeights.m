function y = make_waveletPenalizationWeights(n, n_wind_pos, n_wind_neg)
%% Makes a linear penalization weight starting from the middle of an n-length sample.
%  Anyting after (n/2 + n_wind_pos) gets penalized by a linear taper, likewise for anything before (n/2 - n_wind_neg - 1)
%  ASSUMES EVEN VALUES OF N
%
%  Tim Lin 2010

validateattributes(n,{'numeric'},{'positive','integer','finite','real','nonempty','nonnan','even'})
checkFor_validArrayIndex(n_wind_pos, n/2)
checkFor_validArrayIndex(n_wind_neg, n/2)

y = zeros(n,1);
y(1:n/2-n_wind_neg) = [n/2-n_wind_neg:-1:1];
y(n/2+n_wind_pos+1:n) = [1:n/2-n_wind_pos];
