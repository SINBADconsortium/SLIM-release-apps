
function [f g] = funHUB(r, params)
hub = params.hub;
[f g] = hubers(r,hub,hub);
end
