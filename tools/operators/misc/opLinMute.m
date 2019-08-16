function op = opLinMute(sig_length, mute_start, mute_end)
% function op = opLinMute(length, mute_start, mute_end)
% Apply a window that mutes the shallow water column with a taper
% that starts from mute_start to mute_end.


% Compute the mask

mask = zeros(sig_length,1);
mask(1:mute_start-1) = 0;
mask(mute_end+1:end) = 1;
taper_length = mute_end-mute_start+1;
taper = (1+sin((pi*(0:taper_length-1))/(taper_length-1)-pi/2))/2;
mask(mute_start:mute_end) = taper(:);

% create function handle
subfunc_handle = @(x,mode) opLinMute_intrnl(x,mode);
m = sig_length;
n = m;

op = opFunction(m,n,subfunc_handle); % return a SPOT operator using constructor opFunction


	function y = opLinMute_intrnl(x,mode)
		if mode == 0
		   y = {m,n,[0,1,0,1],{'opLinMute',mask}};

		else
		  y = x.*mask;
		  
		end
	end
end
