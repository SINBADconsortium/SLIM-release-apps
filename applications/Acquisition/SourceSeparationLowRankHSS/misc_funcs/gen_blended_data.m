function [delay_s2, D2_shift, D_blend] = gen_blended_data(rseed, dt, nt, ns, nr, D1, D2) 

%--------------------------------------------------------------------
% gen_blended_data generates blended data for over/under acquisition
%
% Use:
%   gen_blended_data(rseed, dt, nt, ns, nr, D1, D2)
%
% Input:   
%     rseed - random seed
%        dt - time sampling interval
%        nt - number of time samples
%        ns - number of shots
%        nr - number of receivers
%        D1 - original 3-D data matrix from source at shallower depth
%        D2 - original 3-D data matrix from source at deeper depth

% Author: Haneet Wason
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmospheric Sciences
%         The University of British Columbia
%         
% Date: April, 2014

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%--------------------------------------------------------------------

% Random time delays for the second source
% less than 1 second
rng(rseed);
t_s2 = 0:dt:(1-dt);
rperm_s2 = randperm(length(t_s2));
delay_s2 = t_s2(rperm_s2(1:ns));

% Time samples
t = 0:dt:(nt-1)*dt;

% Angular frequency range
wn = [0:nt/2,-nt/2+1:-1]*2*pi/(nt*dt);

% Shifting operator
Ft = opKron(opDirac(nr),opDFT(nt));

% Add delays to second source (apply time-shift)
D2_shift = zeros(nt,nr,ns);
for s = 1:ns
    shot = D2(:,:,s);
    Ftshot = reshape(Ft*shot(:), size(shot));
    for k = 1:size(Ftshot,2)
        Ftshot(:,k) = Ftshot(:,k).*exp(-1i*wn'*delay_s2(s));
    end 
    D2_shift(:,:,s) = reshape(Ft'*Ftshot(:),size(shot));
end

clear shot Ftshot

% Blend shots from both sources
D_blend = D1 + D2_shift;

end

