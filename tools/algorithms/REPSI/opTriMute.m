function op = opTriMute(data,top,t_offset,slope,value)
% OPTRIMUTE     Masks the early times of the traces of a set of 2D semismic data.
%               The mask grows linearly with the offset, making the mask symmetric triangular.
%               Right now the mask assumes that the shots are colocated with the reciever.
%
%               opTriMute(data,top,t_offset,slope) will mask the seismic data cube DATA from t=0 to TOP on all
%               shot records (TOP in number of time samples), as well as mask all samples on a mask growing
%               with offset by SLOPE (defined in absolute number of sample coordinate) with intercept from t=0
%               defined as T_OFFSET. In other words, mask will mute samples with (nt <= top), and also with 
%               (nt < noffset * slope + t_offset), all defined with absolute number of samples.
%                   NOTE: Seismic datacube must be organized as: Coord1 = Time, Coord2 = reciever, Coord3 = shot#



% Compute the mask

Nt = size(data,1);
Nr = size(data,2);
Ns = size(data,3);
clear data

if nargin < 5
    value = 0;
end

Mask = make_mask(Nt,Nr,top,t_offset,slope);

% create function handle
subfunc_handle = @(x,mode) opTriMute_intrnl(x,mode);
m = Nt*Nr*Ns;
n = Nt*Nr*Ns;

op = opFunction(m,n,subfunc_handle); % return a SPOT operator using constructor opFunction


	function y = opTriMute_intrnl(x,mode)
		if mode == 0
		   y = {m,n,[0,1,0,1],{'opTriMute',Mask}};

		else
		  y = reshape(x,Nt*Nr,Ns);
		  for k=1:Ns   % loop over shots
		    MaskTemp  = Mask(:,Nr+1-k:2*Nr-k);
		    y(MaskTemp(:),k)  = value;
		  end
		  y = y(:);
		  
		end
	end
end