function disableUselessWarnings(status)
%DISABLEUSELESSWARNINGS Disable useless warnings in pSpot
%
%   disableUselessWarnings(true) turns off the unecessary warnings forced
%   upon hapless innocent users by Matlab Parallel Computing Toolbox
%
%   disableUselessWarnings(false) turns it back on

if status
    spmd,warning('off','distcomp:codistributed:mTimes:changeOutputCodistr');end
else
    spmd,warning('on','distcomp:codistributed:mTimes:changeOutputCodistr');end
end
