function [out]=MVP_data_aug_syst(in,mode,P,H,lambda)

%   AFUN(X,'notransp') accepts a vector input X and returns the
%   matrix-vector product A*X while AFUN(X,'transp') returns A'*X.

if strcmp(mode,'notransp')==1
    out=lambda*(H*in);
    out=[out;P*in];
elseif strcmp(mode,'transp')==1
    out = lambda*(H'\in(1:size(H,1))) + P'*in(1+size(H,1):end);
    %out=Hmvp(lambda*conj(R),idx,in(1:max(size(R))))+P'*in(max(size(R))+1:end);
end
