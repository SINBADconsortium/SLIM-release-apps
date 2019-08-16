function [out]=proj_card(in,sigma)

%find value of smallest element not to threshold
list=sort(abs(in(:)),'descend');
target=list(sigma);
out=in;
out(abs(out)<target)=0;
