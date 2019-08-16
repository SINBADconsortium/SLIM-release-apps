function y = timeweighting(data,gamma,dt,invflag)
% TIMEWEIGHTING  Apply a time weighting to seismic line datacube
%
%   TIMEWEIGHTING(data,gamma,dt) will apply a time weighting according
%   to exp(gamma*t), where gamma is a positive constant. The function
%   assumes that the first axis is the time axis.
%
%   TIMEWEIGHTING(data,gamma,dt,'inv') will undo the weighting 
%

% Check border case
if (gamma == 0)
    y = data;
    return
end

nt = size(data,1);
time = [0:dt:dt*(nt-1)];
weights = exp(gamma.*time);

if not(exist('invflag','var'))
    for k = 1:nt
        data(k,:,:) = weights(k) .* data(k,:,:);
    end
elseif strcmpi(invflag,'inv')
    for k = 1:nt
        data(k,:,:) = data(k,:,:) ./ weights(k);
    end
else
    error('variable invflag must be a string ''inv''')
end
    
y = data;