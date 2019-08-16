function [s,taxis]= ifftreal(fspec,faxis,outtaxis)
%  [s,taxis]= ifftreal(fspec,faxis)
%
% Input
% fspec    : input frequency spectrum, iFFT along the vertical axis;
% faxis    : frequency axis
% outtaxis : output time axis
%
% Output
% s     :  time domain signal
% taxis : time axis
%
% Author: Xiang Li


[m,n]=size(fspec);
nf = length(faxis);
df = faxis(2) - faxis(1);
if m==1 | n==1
	m = max([m n]);
	fspec = fspec(:);
end
if m < nf
	warning('Signal length is less than frequency axis, will cut taxis')
	faxis(m+1:end) = [];
	nf = m;
elseif m > nf
	warning('Signal length is large than frequency axis, will extend time axis')	
	tmp = m - nf;
	faxis(end:end+m) = faxis(end):df:faxis(end)+df*tmp;
	nf  = m;
end

% padding conjuage part
idx1 = 1:nf;
if isreal(sum(abs(fspec(end,:))))
	idx2 = nf-1:-1:2;
else
	idx2 = nf:-1:2;
end
fspec = [fspec;conj(fspec(idx2,:))];
% iFFt transform
s     = real(ifft(fspec));

% Time axis
nt     = size(s,1);
dt     = 1/(nt*df);
taxis  = dt*(0:nt-1)';


