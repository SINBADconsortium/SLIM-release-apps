function [fspec,faxis]= fftreal(s,taxis,fmin,fmax,outfaxis)
% [fspec,faxis]= fftreal(s,taxis,fmin,fmax,outfaxis)
%
% Input
% s    : input signal, FFT along the vertical axis;
% taxis: time axis
% fmin : minimal output frequecies
% fmax : maximal output frequecies
% outfaxis: output frequencies
%
% Output
% fspec:  frequency spectrum
% faxis: frequency axis
%
% Author: Xiang Li 

intrerp_method = 'v5cubic';
if nargin < 5, outfaxis = 0; end
if nargin < 4, fmax=Inf;end
if nargin < 3, fmin=[]; end

% time interval
dt = taxis(2) - taxis(1);
nt = length(taxis);
[m,n] = size(s);
if m==1 | n==1
	m = max([m n]);
	s = s(:);
end

if m < nt
	warning('Signal length is less than time axis, will cut taxis')
	taxis(m+1:end) = [];
	nt = m;
elseif m > nt
	warning('Signal length is large than time axis, will extend time axis')	
	tmp = m - nt;
	taxis(end:end+m) = taxis(end):dt:taxis(end)+dt*tmp;
	nt  = m;
end	
	
% FFT transfrom
fspec=fft(s);
nf = floor(nt/2+1);
fspec=fspec(1:nf,:);% save only the positive frequencies
clear s;

% Frequency axis
fnyq=1. / (2*(dt));
df=2*fnyq/nt;
faxis=df*(0:nf-1)';


% output certain frequency range
if fmin | fmax < Inf
	idx = find(faxis>=fmin & faxis<=fmax);
	fspec = fspec(idx,:);
	faxis = faxis(idx);
end

% output certain frequencies
if outfaxis ~= 0
	if min(outfaxis) < min(faxis) &  max(outfaxis) > max (faxis)
		error('output frequency range exceed exist frequency range')
	end
	fspec = interp1(faxis,fspec,outfaxis,intrerp_method);
	faxis = outfaxis;
end

